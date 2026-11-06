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


def requested(args):
    backend = os.environ.get('AROMA_HM_BACKEND', 'modeller')
    for i, arg in enumerate(args):
        if arg == '--backend' and i + 1 < len(args):
            backend = args[i + 1]
        elif arg.startswith('--backend='):
            backend = arg.split('=', 1)[1]
    if backend not in ('reec', 'modeller'):
        raise ValueError('Backend must be reec or modeller')
    return backend == 'reec'


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
    p = argparse.ArgumentParser(description='REEC receptor repair' if repair else 'REEC homology modeling')
    p.add_argument('receptor')
    p.add_argument('legacy', nargs='*')
    p.add_argument('--backend', choices=['reec', 'modeller'], default='reec')
    p.add_argument('--template', action='append', default=[], help='PDB template; repeat for an ensemble')
    p.add_argument('--chain', default='A', help='Template protein chain')
    p.add_argument('--alignment', type=Path, help='PIR alignment; target ID is the receptor ID')
    p.add_argument('--contacts', type=Path, action='append', default=[])
    p.add_argument('--output', type=Path, help='Publish here instead of the receptor model')
    p.add_argument('--iterations', type=int, default=6000)
    p.add_argument('--models', type=int, default=2)
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


def known_disulfides(rcp, sources):
    """Preserve declared source links and the project's receptor-specific topology."""
    pairs = set()
    for source in sources:
        if Path(source).exists():
            pairs.update(pdb_disulfides(source))
    bw = lambda value: pu.resno_from_bw(rcp, value)
    first = [bw('3.25'), bw('45.40'), bw('3.40')]
    second = [bw('45.50'), bw('45.60'), bw('5.50')]
    if rcp == 'OR5W2':
        first.append(6)
        second.append(bw('45.38'))
    if rcp in {'OR4D10': 6, 'OR10A7': 3, 'OR51E2': 4}:
        second[1] = {'OR4D10': 6, 'OR10A7': 3, 'OR51E2': 4}[rcp]
    if rcp in {'OR52M1': 8, 'OR5L1': 6, 'OR5L2': 6, 'OR2AT4': 6}:
        first[1] = {'OR52M1': 8, 'OR5L1': 6, 'OR5L2': 6, 'OR2AT4': 6}[rcp]
    if rcp in {'OR56A1': 23, 'OR56A3': 20, 'OR56A4': 19, 'OR56A5': 19}:
        first.append({'OR56A1': 23, 'OR56A3': 20, 'OR56A4': 19, 'OR56A5': 19}[rcp])
        second.append(bw('45.48'))
    if rcp == 'OR1B1':
        first[1] = bw('45.34')
    extra = {'OR1C1': ('3.44', '5.53'), 'OR1N2': ('3.41', '4.49'),
             'OR2T10': ('3.41', '4.49'), 'OR1S1': ('3.55', '5.60'),
             'OR1S2': ('3.55', '5.60'), 'OR2A4': ('3.55', '5.60')}
    if rcp in extra:
        first.append(bw(extra[rcp][0]))
        second.append(bw(extra[rcp][1]))
    if rcp == 'OR11A1':
        second[2] = bw('5.46')
    for a, b in zip(first, second):
        if a and b and pu.aalet_at_resno(rcp, a) == 'C' and pu.aalet_at_resno(rcp, b) == 'C':
            pairs.add(tuple(sorted((a, b))))
    return sorted(pairs)


def receptor_contacts(rcp, sources=()):
    lines = [f'disulfide {a} {b}' for a, b in known_disulfides(rcp, sources)]

    def add(bw1, bw2, atom1, atom2, lo, hi, hard=True):
        a, b = pu.resno_from_bw(rcp, bw1), pu.resno_from_bw(rcp, bw2)
        if not a or not b:
            return
        aa, bb = pu.aalet_at_resno(rcp, a), pu.aalet_at_resno(rcp, b)
        if aa in atom1 and bb in atom2:
            lines.append(f'contact {"hard" if hard else "soft"} {a} {atom1[aa]} '
                         f'{b} {atom2[bb]} {lo} {hi} {100 if hard else 10}')

    add('6.55', '45.51', {'Y': 'OH', 'H': 'ND1|NE2'},
        {'D': 'OD1|OD2', 'E': 'OE1|OE2'}, 2.5, 3.2)
    family = int(re.sub(r'\D', '', pu.family_from_protid(rcp)) or 0)
    if rcp.startswith('OR') and family >= 50:
        add('4.52', '5.50', {'R': 'NE|NH1|NH2', 'K': 'NZ'},
            {'D': 'OD1|OD2', 'E': 'OE1|OE2'}, 2.6, 3.6)
    add('5.58', '7.53', {'Y': 'OH'}, {'Y': 'OH'}, 3.0, 5.0, False)
    return lines


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


def composite_alignment(rcp, work):
    command([sys.executable, HM / 'build_alignment_file.py'])
    text = (HM / 'allgpcr.ali').read_text()
    match = re.search(r'^>P1;' + re.escape(rcp) + r'\s*\n([^\n]*)\n([^*]*)\*', text, re.M)
    if not match:
        raise ValueError(f'Alignment builder omitted {rcp}')
    target = match.group(2) + '*'
    template = work / f'{rcp}_tpl.pdb'
    with directory(HM):
        result = pu.custom_pdb_template(target, rcp, str(template))
    provenance, template_row = result.split('\n', 1)
    # Experimental PIR rows omit trailing gap-only columns. Pad both rows explicitly.
    template_row = ''.join(template_row.split()).rstrip('*')
    target = ''.join(target.split()).rstrip('*')
    width = max(len(template_row), len(target))
    template_row = template_row.ljust(width, '-') + '*'
    target = target.ljust(width, '-') + '*'
    ali = work / 'composite.ali'
    ali.write_text(f'>P1;{rcp}_tpl\nstructure:{rcp}_tpl:FIRST:A:LAST:A::::\n'
                   f'{template_row}\n>P1;{rcp}\n{match.group(1)}\n{target}\n')
    return ali, [f'{rcp}_tpl={template}'], provenance


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


def model(rcp, args, work, templates, alignment, contacts, environment=None, fit_environment=False):
    contact_file = work / 'receptor.contacts'
    contact_file.write_text('\n'.join(contacts) + '\n')
    result = work / 'result.pdb'
    cmd = [ROOT / 'bin/reec', '--alignment', alignment, '--target', rcp,
           '--output', result, '--iterations', args.iterations, '--models', args.models,
           '--contacts', contact_file, '--chain', args.chain]
    for template in templates:
        cmd.extend(['--template', template])
    for contact in args.contacts:
        cmd.extend(['--contacts', contact.resolve()])
    if environment:
        cmd.extend(['--environment', environment])
        if fit_environment:
            cmd.append('--fit-environment')
    try:
        command(cmd)
    except subprocess.CalledProcessError:
        # Retain the numerical evidence even when the candidate cannot be published.
        failure = ROOT / 'tmp' / f'{rcp}.reec-failure.json'
        report = Path(str(result) + '.json')
        if report.exists():
            shutil.copyfile(report, failure)
            print(f'Validation report: {failure}', file=sys.stderr)
        raise
    return result


def run_dohm(argv):
    args = parser().parse_intermixed_args(argv)
    rcp, state = receptor(args.receptor)
    if state != 'active':
        raise ValueError('dohm produces active models; use fixfail for an inactive model')
    unknown = set(args.legacy) - {'nodel', 'riglig'}
    if unknown:
        raise ValueError(f'Unknown options: {sorted(unknown)}')
    family = pu.family_from_protid(rcp)
    source = ROOT / 'pdbs' / family / f'{rcp}.inactive.pdb'
    destination = args.output or ROOT / 'pdbs' / family / f'{rcp}.active.pdb'
    command(['make', 'bin/reec'])
    with workspace(f'{rcp}-reec-', args.keep_work or 'nodel' in args.legacy) as work:
        provenance = ''
        if args.template:
            templates = [template_argument(t) for t in args.template]
            alignment = args.alignment.resolve() if args.alignment else work / 'target.ali'
            if not args.alignment:
                target_alignment(rcp, alignment)
        else:
            if args.alignment:
                raise ValueError('--alignment also requires --template')
            alignment, templates, provenance = composite_alignment(rcp, work)
        environment = (source if 'riglig' in args.legacy and source.exists()
                       and any(s.startswith('HETATM') for s in source.read_text().splitlines()) else None)
        contacts = receptor_contacts(rcp, [source])
        output = model(rcp, args, work, templates, alignment, contacts, environment, bool(environment))
        publish(output, destination, metadata(rcp, source, provenance), orient=rcp)
    print(f'REEC homology model: {destination}')
    return 0


def run_fixfail(argv):
    args = parser(True).parse_intermixed_args(argv)
    rcp, state = receptor(args.receptor)
    family = pu.family_from_protid(rcp)
    original = ROOT / 'pdbs' / family / f'{rcp}.{state}.pdb'
    destination = args.output or original
    options = list(args.legacy)
    input_path = original
    odor = None
    actions = {'reset', 'fit', 'predock', 'loop', 'dock'}
    if options and options[0] not in actions:
        item = options.pop(0)
        candidate = Path(item).resolve()
        if candidate.exists():
            input_path = candidate
        else:
            input_path = ROOT / 'out' / family / rcp / f'{rcp}~{item}.{state}.model1.pdb'
        if not input_path.exists():
            raise ValueError(f'Input PDB not found: {item}')
    if '~' in input_path.name:
        from data import odorutils
        odorutils.load_odors()
        odor = odorutils.find_odorant(input_path.name.split('~', 1)[1].split(f'.{state}.')[0])
    requested_actions = set(options) & actions
    if args.output and requested_actions & {'dock', 'loop'}:
        raise ValueError('Dock actions use the receptor model path; omit --output')
    ligand_tokens = [x for x in options if x not in actions]
    if len(ligand_tokens) % 2:
        raise ValueError('Ligand contacts use pairs: LIGAND_ATOM RES:ATOM')
    if any(x in requested_actions for x in ('fit', 'predock', 'dock', 'loop')) and not odor:
        raise ValueError('Docking and fitting actions require a ligand or docked PDB')
    if 'reset' in requested_actions and state == 'inactive':
        raise ValueError('Reset builds an active model; omit reset for inactive repair')
    if 'reset' in requested_actions:
        reset_argv = [rcp, '--backend', 'reec', '--iterations', str(args.iterations),
                      '--models', str(args.models)]
        run_dohm(reset_argv)
    command(['make', 'bin/reec'] + ([] if args.skip_reshape else ['bin/ic']))
    attempts = 100 if 'loop' in requested_actions else 1
    for attempt in range(attempts):
        print(f'REEC fixfail attempt {attempt + 1}/{attempts}', flush=True)
        with workspace(f'{rcp}-reec-fix-', args.keep_work) as work:
            template = work / 'repair.pdb'
            if 'predock' in requested_actions:
                command(['/bin/bash', './dock.sh', rcp, odor['full_name'], 'noi'])
            if 'fit' in requested_actions:
                from data import odorutils
                sdf = odorutils.ensure_sdf_exists(odor['full_name'])
                command(['make', 'bin/cavity_fit'])
                command([ROOT / 'bin/cavity_fit', original, original.with_suffix('.cvty'),
                         sdf, '-a', '-c', '65536', '-n', '1', '-o', template])
            else:
                shutil.copyfile(input_path, template)
            if not args.skip_reshape:
                manifest = ('OR_ClassI' if rcp.startswith('OR') and int(re.sub(r'\D', '', family)) >= 50
                            else 'OR_ClassII' if rcp.startswith('OR') else 'TAAR' if rcp.startswith('TAAR')
                            else 'VN1R' if rcp.startswith('VN1R') else None)
                if not manifest:
                    raise ValueError('No receptor reshape manifest; use --skip-reshape')
                shape = ROOT / 'data' / f'{manifest}_{"i" if state == "inactive" else "a"}.ic'
                command([ROOT / 'bin/ic', template, shape, 'save'])
            alignment = args.alignment.resolve() if args.alignment else work / 'target.ali'
            if not args.alignment:
                target_alignment(rcp, alignment, template, args.chain)
            contacts = receptor_contacts(rcp, [original, input_path])
            ligand_lines = [s for s in template.read_text().splitlines() if s.startswith('HETATM')]
            environment = template if ligand_lines else None
            if ligand_tokens:
                if not environment:
                    raise ValueError('Ligand contacts require HETATM records')
                for ligand_atom, receptor_atom in zip(ligand_tokens[::2], ligand_tokens[1::2]):
                    matches = {int(s[22:26]) for s in ligand_lines if s[12:16].strip() == ligand_atom}
                    if len(matches) != 1:
                        raise ValueError(f'Ambiguous or absent ligand atom: {ligand_atom}')
                    residue, atom = receptor_atom.split(':', 1)
                    contacts.append(f'contact hard {residue} {atom} {matches.pop()} {ligand_atom} 2.7 3.8 100')
            templates = [str(template)] + [template_argument(t) for t in args.template]
            if original.exists():
                templates.append(str(original))
            output = model(rcp, args, work, templates, alignment, contacts, environment)
            publish(output, destination, metadata(rcp, original))
        if 'dock' not in requested_actions and 'loop' not in requested_actions:
            break
        if args.output:
            raise ValueError('Dock actions use the receptor model path; omit --output')
        command(['/bin/bash', './dock.sh', rcp, odor['full_name'], 'noi'])
        dock = ROOT / 'out' / family / rcp / f'{rcp}~{odor["full_name"]}.{state}.dock'
        if dock.exists():
            text = dock.read_text()
            poses = re.search(r'(\d+) pose\(s\) found', text)
            score = re.search(r'^Total:\s+(-?[\d.]+)', text, re.M)
            if poses and int(poses[1]) >= 4 and score and float(score[1]) < 0 and 'Disqualified because:' not in text:
                print('SUCCESS')
                break
    print(f'REEC repaired model: {destination}')
    return 0


def entrypoint(function, argv):
    try:
        (ROOT / 'tmp').mkdir(exist_ok=True)
        return function(argv)
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        print(f'REEC workflow: {error}', file=sys.stderr)
        return 1
