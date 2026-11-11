"""Geometry repair and explicit fallback publication for batch model refreshes."""
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

from data.disulfides import pdb_disulfides
from reec_backend import CandidateRejected
from reec_reports import report_path
from reec_review import inspect, inspect_sterics, dihedral


CODES = dict(zip('ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR'.split(),
                 'ACDEFGHIKLMNPQRSTVWY'))


def review_source(source, target, sequence_length, work):
    candidate = Path(work) / 'source-review.pdb'
    shutil.copyfile(source, candidate)
    Path(str(candidate)+'.json').write_text(json.dumps({'contacts': [], 'target_sequence_length': sequence_length}))
    try:
        review = inspect(candidate)
        if review['valid']:
            sterics = inspect_sterics(candidate)
            review['serialized_sterics'] = sterics
            if sterics['severe_clashes']:
                review['valid'] = False
                review['failures'].append(f'Final serialized model has {sterics["severe_clashes"]} severe clashes')
        return review
    except (ValueError, KeyError) as exc:
        return {'valid': False, 'failures': ['Source geometry could not be inspected'],
                'inspection_error_type': type(exc).__name__}


def run_native(command, root, work, label):
    from reec_thermal import run_guarded, read_temperature
    environment = os.environ.copy()
    environment.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    try:
        read_temperature()
    except (OSError, ValueError, RuntimeError, subprocess.TimeoutExpired):
        return subprocess.run([str(v) for v in command], cwd=root, env=environment).returncode
    with (work / f'{label}.log').open('w') as stream:
        result = run_guarded([str(v) for v in command], cwd=root, environment=environment,
                             log=stream, temperature_log=work / f'{label}.temperature.jsonl',
                             timeout=1800, pause_at=67, resume_below=61)
    (work / f'{label}.execution.json').write_text(json.dumps(result, indent=2)+'\n')
    return result['exit_code']


def repair_existing(target, source, options, root, *, contact_report=None, orient=True):
    """Repair full geometry while softly retaining the source side-chain torsions."""
    import data.protutils as pu
    from reec_backend import publish, metadata
    from modeling_backend import publication_guard, finalize
    root, source = Path(root), Path(source).resolve()
    (root / 'tmp').mkdir(exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix=f'{target}-geometry-repair-', dir=root / 'tmp'))
    snapshot = work / 'source.pdb'
    shutil.copyfile(source, snapshot)
    residues, atoms = {}, {}
    for line in snapshot.read_text().splitlines():
        if line.startswith('ATOM  ') and line[21] == 'A':
            number, name = int(line[22:26]), line[12:16].strip()
            if line[17:20] not in CODES:
                raise CandidateRejected('Geometry repair source has an unsupported protein residue')
            residues[number] = CODES[line[17:20]]
            atoms[number, name] = tuple(float(line[i:i+8]) for i in (30, 38, 46))
    sequence = pu.prots[target]['sequence']
    if not residues or max(residues) > len(sequence):
        raise CandidateRejected('Repair source does not match the target residue numbering')
    row = ''.join(residues.get(i, '-') for i in range(1, len(sequence)+1))
    alignment = work / 'repair.ali'
    alignment.write_text(f'>P1;repair_source\nstructureX:repair_source:FIRST:A:LAST:A::::\n{row}*\n'
                         f'>P1;{target}\nsequence:{target}:1:A::::::\n{sequence}*\n')
    contacts = [f'disulfide {a} {b}' for a, b in pdb_disulfides(snapshot)]
    for contact in contact_report or []:
        if contact['kind'] == 'disulfide':
            continue
        mode, weight = ('hard', 100) if contact['hard'] else ('soft', 10)
        contacts.append(f'contact {mode} {contact["res1"]} {contact["atom1"]} '
                        f'{contact["res2"]} {contact["atom2"]} '
                        f'{contact["min"]} {contact["max"]} {weight}')
    # General source-conformation priors; no receptor-specific pocket exception.
    for number in residues:
        names = {a for n, a in atoms if n == number}
        branch = next((a for a in ('CG', 'CG1', 'OG', 'OG1', 'SG') if a in names), None)
        quads = [('N', 'CA', 'CB', branch)] if branch else []
        if branch == 'CG':
            tip = next((a for a in ('CD', 'CD1', 'SD', 'OD1', 'ND1') if a in names), None)
            if tip:
                quads.append(('CA', 'CB', 'CG', tip))
        for quad in quads:
            if all((number, a) in atoms for a in quad):
                mean = math.radians(dihedral(*(atoms[number, a] for a in quad)))
                contacts.append('dihedral '+ ' '.join(f'{number} {a}' for a in quad)+f' {mean:.10f} 0.5')
    contact_file = work / 'repair.contacts'
    contact_file.write_text('\n'.join(contacts)+'\n')
    last_report = None
    for attempt in range(1, 4):
        candidate = work / f'attempt-{attempt}.pdb'
        command = [root / 'bin/reec', '--alignment', alignment, '--target', target,
                   '--template', f'repair_source={snapshot}', '--contacts', contact_file,
                   '--output', candidate, '--models', '1', '--iterations',
                   str(max(6000, options.iterations or 0)*attempt), '--seed', str(options.seed+attempt-1),
                   '--anneal-steps', '0', '--nofail']
        if any(line.startswith('HETATM') for line in snapshot.read_text().splitlines()):
            command.extend(['--environment', snapshot, '--environment-hydrogens'])
        print(f'Repairing {target} geometry: attempt {attempt}/3; work {work}', file=sys.stderr, flush=True)
        code = run_native(command, root, work, f'attempt-{attempt}')
        last_report = Path(str(candidate)+'.json')
        if candidate.exists():
            options.last_rejected_model = candidate
        if code != 0:
            continue
        if not inspect(candidate)['valid']:
            continue
        oriented = work / f'attempt-{attempt}.oriented.pdb'
        publish(candidate, oriented, metadata(target, source, f'geometry-repaired {source.name}'),
                orient=target if orient else None)
        hydrated = work / 'repaired.pdb'
        script = work / 'hydrate.phew'
        script.write_text(f'LOAD "{oriented}" A A\nSTRAND A\nHYDRO\nSAVE "{hydrated}"\n')
        with publication_guard(hydrated):
            code = run_native([root / 'bin/phew', script], root, work, f'hydrogenate-{attempt}')
            if code != 0:
                continue
            try:
                finalize(None, hydrated, candidate, target=target)
            except CandidateRejected:
                continue
        status_path = Path(str(hydrated)+'.json')
        status = json.loads(status_path.read_text())
        status['validated_ready'] = True
        status['validated_sha256'] = hashlib.sha256(hydrated.read_bytes()).hexdigest()
        status['geometry_repair_work'] = str(work)
        status_path.write_text(json.dumps(status, indent=2)+'\n')
        return hydrated, work
    raise CandidateRejected(f'Geometry repair did not pass in three attempts; work {work}', last_report)


def copy_pair(source, destination):
    """Atomically replace each file; outer publication guards protect the pair."""
    source, destination = Path(source), Path(destination)
    if not source.is_file():
        raise FileNotFoundError(source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    for suffix in ('', '.json'):
        original = report_path(source) if suffix else source
        final = Path(str(destination)+suffix)
        if not original.is_file():
            if suffix:
                final.unlink(missing_ok=True)
            continue
        temporary = final.with_name(final.name+f'.copy-{os.getpid()}')
        try:
            shutil.copyfile(original, temporary)
            os.replace(temporary, final)
        finally:
            temporary.unlink(missing_ok=True)


def fallback_model(target, destination, options, rejection, root):
    """Emit a finite protein model and clearly retain its quality status."""
    import data.protutils as pu
    from modeling_backend import publication_guard
    root, destination = Path(root), Path(destination)
    family = pu.family_from_protid(target)
    candidates = [getattr(options, 'last_rejected_model', None), destination,
                  root / f'pdbs/{family}/{target}.active.pdb',
                  root / f'pdbs/_validated/{target}.active.pdb',
                  root / f'pdbs/{family}/{target}.inactive.pdb']
    for source in candidates:
        if source is None or not Path(source).is_file():
            continue
        try:
            text = Path(source).read_text()
            points = [tuple(float(line[i:i+8]) for i in (30, 38, 46))
                      for line in text.splitlines() if line.startswith('ATOM  ')]
        except (OSError, ValueError):
            continue
        if not points or not all(math.isfinite(v) for p in points for v in p):
            continue
        report = {'valid': False, 'status': 'fallback', 'refresh_accepted': False,
                  'reason': str(rejection), 'source': str(source), 'contacts': [],
                  'target_sequence_length': len(pu.prots[target]['sequence'])}
        log_root = root / 'tmp'
        log_root.mkdir(exist_ok=True)
        work = Path(tempfile.mkdtemp(prefix=f'{target}-fallback-', dir=log_root))
        backup = work / 'source.pdb'
        shutil.copyfile(source, backup)
        try:
            review = review_source(backup, target, len(pu.prots[target]['sequence']), work)
        except (ValueError, KeyError, OSError):
            review = {'valid': False, 'failures': ['Independent geometry review could not complete']}
        source_report = report_path(source)
        try:
            prior_status = json.loads(source_report.read_text()) if source_report.is_file() else {}
        except (ValueError, OSError):
            prior_status = {}
        report['valid'] = review['valid'] and prior_status.get('valid', False)
        report['serialized_geometry'] = {k: v for k, v in review.items() if k not in ('model', 'optimizer_report')}
        from reec_backend import metadata, receptor_frame
        try:
            text = receptor_frame(target, text)
            report['receptor_frame'] = True
        except (ValueError, KeyError):
            report['receptor_frame'] = False
        text = '\n'.join(metadata(target, source, 'batch fallback'))+'\n'+text
        with publication_guard(destination):
            destination.write_text('REMARK 250 DOHM FALLBACK: refresh did not pass; inspect sidecar JSON.\n'+text)
            Path(str(destination)+'.json').write_text(json.dumps(report, indent=2)+'\n')
        print(f'Fallback protein model written: {destination}; quality status: {report_path(destination)}.',
              file=sys.stderr, flush=True)
        return
    raise CandidateRejected('No finite protein candidate or prior receptor is available for fallback') from rejection
