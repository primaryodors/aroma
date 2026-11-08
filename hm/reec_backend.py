"""REEC integration shared by dohm and fixfail; MODELLER stays selectable."""
import argparse
import contextlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
HM = ROOT / 'hm'
sys.path.insert(0, str(ROOT))
from data import protutils as pu
from data.disulfides import pdb_disulfides


class CandidateRejected(ValueError):
    """A modeled candidate failed validation and can be retried."""
    def __init__(self, message, report=None):
        super().__init__(message)
        self.report = Path(report) if report is not None else None


def requested(args, default='modeller'):
    backend = os.environ.get('AROMA_HM_BACKEND', default)
    for i, arg in enumerate(args):
        if arg == '--backend' and i + 1 < len(args):
            backend = args[i + 1]
        elif arg.startswith('--backend='):
            backend = arg.split('=', 1)[1]
    if backend not in ('reec', 'modeller'):
        raise ValueError('Backend must be reec or modeller')
    return backend == 'reec'


def legacy_backend_args(args):
    result = []
    index = 0
    while index < len(args):
        if args[index] == '--backend':
            index += 2
        elif args[index].startswith('--backend='):
            index += 1
        else:
            result.append(args[index])
            index += 1
    return result


@contextlib.contextmanager
def directory(path):
    previous = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


def command(args):
    print(' '.join(str(x) for x in args), flush=True)
    subprocess.run([str(x) for x in args], cwd=ROOT, check=True)


@contextlib.contextmanager
def workspace(prefix, keep=False):
    path = Path(tempfile.mkdtemp(prefix=prefix, dir=ROOT / 'tmp'))
    try:
        yield path
    except BaseException:
        print(f'REEC work files retained: {path}', file=sys.stderr)
        raise
    else:
        if keep:
            print(f'REEC work files retained: {path}')
        else:
            shutil.rmtree(path)


def parser(repair=False):
    p = argparse.ArgumentParser(description=('Receptor repair' if repair else 'Receptor homology modeling')
                                + ' (REEC or MODELLER)')
    p.add_argument('receptor')
    p.add_argument('legacy', nargs='*')
    p.add_argument('--backend', choices=['reec', 'modeller'], default='reec')
    p.add_argument('--template', action='append', default=[], help='PDB template; repeat for an ensemble')
    p.add_argument('--chain', default='A', help='Template protein chain')
    p.add_argument('--alignment', type=Path, help='PIR alignment; target ID is the receptor ID')
    p.add_argument('--contacts', type=Path, action='append', default=[])
    p.add_argument('--output', type=Path, help='Publish here instead of the receptor model')
    p.add_argument('--iterations', type=int, help='Initial modeling iteration budget')
    p.add_argument('--models', type=int, help='Override the workflow ensemble size')
    p.add_argument('--seed', type=int, default=0, help='Reproducible REEC starting-conformation seed')
    p.add_argument('--anneal-steps', type=int, default=0, help='Targeted Monte Carlo rounds; disabled by default, use 64 to enable')
    p.add_argument('--anneal-temperature', type=float, default=12, help='Initial annealing temperature in REEC score units')
    p.add_argument('--anneal-residues', help='Comma-separated residue numbers; defaults to receptor binding-site residues')
    p.add_argument('--max-attempts', type=int, default=100 if repair else 0,
                   help='Maximum attempts; dohm defaults to retrying until success (0)')
    p.add_argument('--keep-work', action='store_true', help='Retain PIR, template, and contact files')
    p.add_argument('--skip-reshape', action='store_true', help='Repair the supplied conformation directly')
    return p


def receptor(value):
    state = 'inactive' if value.endswith('i') else 'active'
    if value.endswith(('a', 'i')):
        value = value[:-1]
    value = value.upper()
    pu.load_prots()
    if value not in pu.prots:
        raise ValueError(f'Unknown receptor: {value}')
    return value, state


def template_argument(value):
    if '=' in value:
        name, path = value.split('=', 1)
        return f'{name}={Path(path).resolve()}'
    return str(Path(value).resolve())


def target_alignment(rcp, path, template=None, chain='A'):
    sequence = pu.prots[rcp]['sequence']
    prefix = ''
    if template:
        codes = dict(zip(
            'ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR'.split(),
            'ACDEFGHIKLMNPQRSTVWY'))
        present = set()
        numbered = True
        for line in Path(template).read_text().splitlines():
            if not line.startswith('ATOM  ') or line[12:16].strip() != 'CA' or line[21] not in (' ', chain):
                continue
            number = int(line[22:26])
            if not 1 <= number <= len(sequence) or codes.get(line[17:20]) != sequence[number-1] or number in present:
                numbered = False
                break
            present.add(number)
        if numbered and present:
            row = ''.join(amino if i in present else '-' for i, amino in enumerate(sequence, 1))
            name = Path(template).stem
            prefix = f'>P1;{name}\nstructureX:{name}:FIRST:{chain}::::::\n{row}*\n'
    path.write_text(prefix + f'>P1;{rcp}\nsequence:{rcp}:1:A::::::\n{sequence}*\n')


def validated_alignment(rcp, work, chain='A'):
    if not re.fullmatch(r'OR5K\d+', rcp):
        return None
    source = ROOT / 'pdbs' / '_validated' / 'OR5K1.active.pdb'
    if not source.is_file():
        raise ValueError(f'Validated OR5K1 template is missing: {source}')
    rows = {}
    for line in (ROOT / 'data' / 'sequences_aligned.txt').read_text().splitlines():
        parts = line.split()
        if parts and parts[0] in ('OR5K1', rcp):
            rows[parts[0]] = parts[1]
    for name in ('OR5K1', rcp):
        if name not in rows or rows[name].replace('-', '') != pu.prots[name]['sequence']:
            raise ValueError(f'GPCR alignment does not match the sequence for {name}')
    codes = dict(zip(
        'ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR'.split(),
        'ACDEFGHIKLMNPQRSTVWY'))
    present = {}
    for line in source.read_text().splitlines():
        if not line.startswith('ATOM  ') or line[12:16].strip() != 'CA' or line[21] not in (' ', chain):
            continue
        number = int(line[22:26])
        sequence = pu.prots['OR5K1']['sequence']
        if (line[26].strip() or number in present or not 1 <= number <= len(sequence)
                or codes.get(line[17:20]) != sequence[number-1]):
            raise ValueError('Validated OR5K1 residue numbering does not match its sequence')
        present[number] = True
    if not present:
        raise ValueError(f'Validated OR5K1 has no CA atoms in chain {chain}')
    position = 0
    template_row = ''
    for amino in rows['OR5K1']:
        if amino != '-':
            position += 1
        template_row += amino if amino == '-' or position in present else '-'
    width = max(len(template_row), len(rows[rcp]))
    template = work / f'{rcp}_tpl.pdb'
    shutil.copyfile(source, template)
    ali = work / f'{rcp}.validated.ali'
    ali.write_text(f'>P1;{rcp}_tpl\nstructureX:{rcp}_tpl:FIRST:{chain}:LAST:{chain}::::\n'
                   f'{template_row.ljust(width, "-")}*\n>P1;{rcp}\nsequence:{rcp}:1:A::::::\n'
                   f'{rows[rcp].ljust(width, "-")}*\n')
    return ali, [f'{rcp}_tpl={template}'], 'pdbs/_validated/OR5K1.active.pdb'


def metadata(rcp, source=None, provenance=''):
    lines = []
    if source and Path(source).exists():
        lines = [s for s in Path(source).read_text().splitlines()
                 if s.startswith(('REMARK   1 ', 'REMARK 265 HM_TEMPLATES'))]
    if provenance:
        lines = [s for s in lines if not s.startswith('REMARK 265 HM_TEMPLATES')]
        lines.append(f'REMARK 265 HM_TEMPLATES: {provenance}')
    for name, region in pu.prots[rcp].get('region', {}).items():
        if name.startswith(('TMR', 'HXR')):
            lines.append(f'REMARK 650 HELIX {name} {region["start"]} {region["end"]}')
    for tag, residue in pu.prots[rcp].get('bw', {}).items():
        lines.append(f'REMARK 800 SITE BW {tag} {residue}')
    return lines


def receptor_frame(rcp, text):
    """Use the same membrane and BW anchors as UPRIGHT/BWCENTER, by rigid rotation."""
    from reec_review import sub, dot, cross

    def unit(v):
        length = math.sqrt(dot(v, v))
        if length < 1e-8:
            raise ValueError('Degenerate receptor orientation anchors')
        return tuple(x / length for x in v)

    def mean(points):
        return tuple(sum(p[i] for p in points) / len(points) for i in range(3))

    lines = text.splitlines()
    ca = {int(s[22:26]): tuple(float(s[i:i+8]) for i in (30, 38, 46))
          for s in lines if s.startswith('ATOM  ') and s[12:16].strip() == 'CA'}
    extracellular, cytoplasmic = [], []
    regions = pu.prots[rcp]['region']
    if any(f'TMR{h}' not in regions for h in range(1, 8)):
        raise ValueError('Receptor orientation requires TMR1 through TMR7 annotations')
    for helix in range(1, 8):
        region = regions[f'TMR{helix}']
        start = [ca[n] for n in range(region['start'], region['start'] + 4)]
        end = [ca[n] for n in range(region['end'] - 3, region['end'] + 1)]
        extracellular.extend(start if helix % 2 else end)
        cytoplasmic.extend(end if helix % 2 else start)
    y = unit(sub(mean(extracellular), mean(cytoplasmic)))
    centers = []
    for name in ('TMR1', 'TMR4'):
        region = regions[name]
        centers.append(mean([ca[n] for n in range(region['start'], region['end'] + 1)]))
    direction = sub(centers[1], centers[0])
    z = unit(sub(direction, tuple(c * dot(direction, y) for c in y)))
    x = cross(y, z)
    center = mean([ca[pu.resno_from_bw(rcp, f'{h}.50')] for h in range(1, 8)])
    for i, line in enumerate(lines):
        if line.startswith(('ATOM  ', 'HETATM')):
            point = sub(tuple(float(line[j:j+8]) for j in (30, 38, 46)), center)
            values = [dot(point, axis) for axis in (x, y, z)]
            if any(v < -999.999 or v > 9999.999 for v in values):
                raise ValueError('Oriented coordinates exceed PDB capacity')
            lines[i] = line[:30] + ''.join(f'{v:8.3f}' for v in values) + line[54:]
    return '\n'.join(lines) + '\n'


def publish(model, destination, remarks, orient=None):
    from reec_review import inspect
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    body = model.read_text()
    if orient:
        body = receptor_frame(orient, body)
    text = '\n'.join(remarks) + '\n' + body
    temp = destination.with_name(destination.name + f'.reec-{os.getpid()}')
    temp_report = Path(str(temp) + '.json')
    report = json.loads(Path(str(model) + '.json').read_text())
    if not report.get('valid'):
        raise ValueError('A failed model cannot be published')
    # Re-evaluate selected contact distances after PDB coordinate rounding.
    atoms = {(int(s[22:26]), s[12:16].strip()): tuple(float(s[i:i+8]) for i in (30, 38, 46))
             for s in body.splitlines() if s.startswith(('ATOM  ', 'HETATM'))}
    for contact in report['contacts']:
        distance = math.dist(atoms[contact['res1'], contact['atom1']],
                             atoms[contact['res2'], contact['atom2']])
        contact['distance'] = distance
        contact['violation'] = max(0, contact['min'] - distance, distance - contact['max'])
    report['mandatory_violations'] = sum(c['hard'] and c['violation'] > .05 for c in report['contacts'])
    if report['mandatory_violations']:
        raise ValueError('Mandatory contacts failed after PDB coordinate rounding')
    try:
        temp.write_text(text)
        temp_report.write_text(json.dumps(report, indent=2) + '\n')
        review = inspect(temp)
        if not review['valid']:
            raise ValueError(f'Serialized PDB failed independent review: {review["failures"]}')
        report['serialized_geometry'] = {k: v for k, v in review.items()
                                         if k not in ('model', 'optimizer_report')}
        report['receptor_frame'] = bool(orient)
        temp_report.write_text(json.dumps(report, indent=2) + '\n')
        os.replace(temp, destination)
        os.replace(temp_report, str(destination) + '.json')
    finally:
        temp.unlink(missing_ok=True)
        temp_report.unlink(missing_ok=True)


def model(rcp, args, work, templates, alignment, contacts, environment=None, fit_environment=False,
          environment_hydrogens=False):
    contact_file = work / 'receptor.contacts'
    contact_file.write_text('\n'.join(contacts) + '\n')
    result = work / 'result.pdb'
    cmd = [ROOT / 'bin/reec', '--alignment', alignment, '--target', rcp,
           '--output', result, '--iterations', args.iterations, '--models', args.models,
           '--contacts', contact_file, '--chain', args.chain, '--seed', getattr(args, 'seed', 0)]
    cmd.extend(['--anneal-steps', getattr(args, 'anneal_steps', 0),
                '--anneal-temperature', getattr(args, 'anneal_temperature', 12),
                '--clash-scale', getattr(args, 'clash_scale', 1.)])
    if getattr(args, 'anneal_residues', None):
        cmd.extend(['--anneal-residues', args.anneal_residues])
    for template in templates:
        cmd.extend(['--template', template])
    if getattr(args, 'candidate_prefix', None):
        cmd.extend(['--candidate-prefix', args.candidate_prefix])
    for contact in args.contacts:
        cmd.extend(['--contacts', contact.resolve()])
    if environment:
        cmd.extend(['--environment', environment])
        if fit_environment:
            cmd.append('--fit-environment')
        if environment_hydrogens:
            cmd.append('--environment-hydrogens')
    try:
        command(cmd)
    except subprocess.CalledProcessError as exc:
        # Retain the numerical evidence even when the candidate cannot be published.
        failure = ROOT / 'tmp' / f'{rcp}.reec-failure.json'
        report = Path(str(result) + '.json')
        if report.exists():
            shutil.copyfile(report, failure)
            print(f'Validation report: {failure}', file=sys.stderr)
        if exc.returncode == 2:
            raise CandidateRejected('REEC rejected the candidate ensemble', report if report.exists() else None) from exc
        raise
    return result


def run_dohm(argv):
    from dohm import main
    return main(argv)


def run_fixfail(argv):
    from fixfail import main
    return main(argv)


def entrypoint(function, argv):
    try:
        (ROOT / 'tmp').mkdir(exist_ok=True)
        return function(argv)
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        print(f'Modeling workflow: {error}', file=sys.stderr)
        return 1
