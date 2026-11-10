"""Model-engine adapters; workflow orchestration stays in dohm.py and fixfail.py."""
from dataclasses import dataclass
from contextlib import contextmanager
import importlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
from types import SimpleNamespace

import reec_backend as reec


def options(argv, repair=False):
    args = reec.parser(repair).parse_intermixed_args(
        ['--help' if value == 'help' else value for value in argv])
    args.backend = 'reec' if reec.requested(argv, default='reec') else 'modeller'
    if args.iterations is not None and args.iterations < 1:
        raise ValueError('Iterations must be positive')
    if args.models is not None and args.models < 1:
        raise ValueError('Models must be positive')
    if not 0 <= args.seed <= 2147483647:
        raise ValueError('Seed must be between 0 and 2147483647')
    if args.max_attempts < (1 if repair else 0):
        raise ValueError('Invalid maximum attempt count')
    if not 0 <= args.anneal_steps <= 100000:
        raise ValueError('Annealing steps must be between 0 and 100000')
    if not math.isfinite(args.anneal_temperature) or args.anneal_temperature <= 0:
        raise ValueError('Annealing temperature must be positive and finite')
    if args.anneal_residues is not None and not re.fullmatch(r'[1-9][0-9]*(,[1-9][0-9]*)*', args.anneal_residues):
        raise ValueError('Annealing residues must be comma-separated positive integers')
    for name in ('alignment', 'output'):
        value = getattr(args, name)
        if value is not None:
            setattr(args, name, value.resolve())
    args.contacts = [value.resolve() for value in args.contacts]
    args.template = [reec.template_argument(value) for value in args.template]
    if not repair and args.alignment and not args.template:
        raise ValueError('--alignment also requires --template')
    return args


def forwarded_options(args):
    result = ['--backend', args.backend, '--chain', args.chain]
    for name in ('iterations', 'models', 'alignment', 'seed', 'anneal_steps', 'anneal_temperature', 'anneal_residues'):
        value = getattr(args, name)
        if value is not None:
            result.extend(['--'+name.replace('_', '-'), str(value)])
    for name in ('template', 'contacts'):
        for value in getattr(args, name):
            result.extend(['--'+name, str(value)])
    if args.keep_work:
        result.append('--keep-work')
    if getattr(args, 'coordinate_metals', False):
        result.append('--coordinate-metals')
    return result


def pir_rows(path):
    return {name.strip(): (description.strip(), ''.join(row.split()))
            for name, description, row in re.findall(
                r'^>P1;([^\n]+)\n([^\n]+)\n([^*]*)\*', Path(path).read_text(), re.M)}


def custom_alignment(args, receptor, work):
    """Prepare the same explicit templates and PIR rows for either engine."""
    from Bio.Align import PairwiseAligner, substitution_matrices
    codes = dict(zip('ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR'.split(),
                     'ACDEFGHIKLMNPQRSTVWY'))
    templates, sequences = [], {}
    for item in args.template:
        if '=' in item:
            name, filename = item.split('=', 1)
        else:
            filename, name = item, Path(item).stem
        if not re.fullmatch(r'[A-Za-z0-9_.-]+', name) or name == receptor:
            raise ValueError(f'Invalid template ID: {name}')
        source = Path(filename)
        target = work / (name+'.pdb')
        if source.resolve() != target.resolve():
            shutil.copyfile(source, target)
        sequences[name] = ''.join(codes[line[17:20]] for line in target.read_text().splitlines()
                                 if line.startswith('ATOM  ') and line[12:16].strip() == 'CA'
                                 and line[21] in (' ', args.chain))
        if not sequences[name]:
            raise ValueError(f'Template has no protein CA atoms: {source}')
        templates.append(f'{name}={target}')
    if args.alignment:
        rows = pir_rows(args.alignment)
        if receptor in rows and all(name in rows for name in sequences):
            return args.alignment, templates, ' '.join(args.template)
    target_sequence = reec.pu.prots[receptor]['sequence']
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
    aligner.open_gap_score, aligner.extend_gap_score = -10, -.5
    aligned, widths = {}, [0]*(len(target_sequence)+1)
    for name, sequence in sequences.items():
        pair = aligner.align(target_sequence, sequence)[0]
        insertions, residues, position = ['']*len(widths), [], 0
        for target, template in zip(pair[0], pair[1]):
            if target == '-':
                insertions[position] += template
            else:
                residues.append(template)
                position += 1
        widths = [max(width, len(chars)) for width, chars in zip(widths, insertions)]
        aligned[name] = insertions, residues
    output = work / f'{receptor}.custom.ali'
    text = ''
    for name, (insertions, residues) in aligned.items():
        row = ''.join(insertions[i].ljust(widths[i], '-') + residues[i] for i in range(len(target_sequence)))
        row += insertions[-1].ljust(widths[-1], '-')
        text += f'>P1;{name}\nstructure:{name}:FIRST:{args.chain}:LAST:{args.chain}::::\n{row}*\n'
    row = ''.join('-'*widths[i]+amino for i, amino in enumerate(target_sequence))+'-'*widths[-1]
    output.write_text(text+f'>P1;{receptor}\nsequence:{receptor}:1:A::::::\n{row}*\n')
    return output, templates, ' '.join(args.template)


def protein_alignment(source, knowns, target, destination):
    """REEC receives protein columns; its fixed environment carries ligand atoms."""
    rows = pir_rows(source)
    selected = {name: rows[name] for name in (*knowns, target) if name in rows}
    if target not in selected:
        raise ValueError(f'Target is absent from alignment: {target}')
    width = max(len(row) for _, row in selected.values())
    padded = {name: (description, row.ljust(width, '-')) for name, (description, row) in selected.items()}
    target_row = padded[target][1]
    columns = [i for i, amino in enumerate(target_row) if amino != '.']
    text = ''
    for name, (description, row) in padded.items():
        protein = ''.join(row[i] for i in columns)
        if '.' in protein or '/' in protein:
            raise ValueError('REEC accepts a single protein chain plus a separate ligand environment')
        text += f'>P1;{name}\n{description}\n{protein}*\n'
    destination.write_text(text)
    return destination


def api(args):
    if args.backend == 'modeller':
        modeller = importlib.import_module('modeller')
        automodel = importlib.import_module('modeller.automodel')
        return SimpleNamespace(**{name: getattr(modeller, name) for name in
            ('Environ', 'Selection', 'RigidBody', 'secondary_structure', 'forms', 'features', 'physical')},
            AutoModel=automodel.AutoModel, DOPEHRLoopModel=automodel.DOPEHRLoopModel,
            autosched=automodel.autosched, refine=automodel.refine)
    return SimpleNamespace(Environ=lambda: Environment(args), AutoModel=ReecModel,
        DOPEHRLoopModel=ReecModel, Selection=lambda *ranges: ranges, RigidBody=lambda value: value,
        secondary_structure=SimpleNamespace(Alpha=lambda value: ('alpha', value)),
        forms=SimpleNamespace(Gaussian=Gaussian),
        features=SimpleNamespace(Distance=lambda *atoms: ('distance', atoms),
                                 Dihedral=lambda *atoms: ('dihedral', atoms)),
        physical=SimpleNamespace(xy_distance='distance', chi1_dihedral='dihedral', chi2_dihedral='dihedral'),
        autosched=SimpleNamespace(slow='slow'), refine=SimpleNamespace(very_slow='very_slow'))


def prepare(args):
    """Build the selected native engine before workflow preprocessing starts."""
    if args.backend == 'reec' and not getattr(args, '_reec_built', False):
        reec.command(['make', 'bin/reec', 'bin/reec_steric_audit'])
        args._reec_built = True


@dataclass(frozen=True)
class Atom:
    name: str
    residue: int


class Atoms:
    def __getitem__(self, value):
        name, residue, *_ = value.split(':')
        return Atom(name, int(residue))


class Residues:
    def __init__(self, model):
        self.model = model

    def __getitem__(self, value):
        if value == -1:
            records = self.model.environment_records()
            residue = int(records[-1][22:26]) if records else len(reec.pu.prots[self.model.sequence]['sequence'])
        else:
            residue = int(str(value).split(':')[0])
        return SimpleNamespace(intnum=residue)


@dataclass(frozen=True)
class Gaussian:
    group: object
    feature: tuple
    mean: float
    stdev: float


class Restraints:
    def __init__(self):
        self.items = []
        self.rigid_bodies = []

    def add(self, value):
        self.items.append(value)


class Environment:
    def __init__(self, args):
        self.options = args
        self.io = SimpleNamespace(atom_files_directory=['.'], hetatm=False)


class ReecModel:
    """Execute the scripts' model callbacks against the REEC numerical engine."""
    def __init__(self, env, alnfile, knowns, sequence):
        self.env, self.alnfile, self.sequence = env, Path(alnfile).resolve(), sequence
        self.knowns = (knowns,) if isinstance(knowns, str) else tuple(knowns)
        self.starting_model, self.ending_model, self.max_var_iterations = 0, 9, 1000
        self.atoms, self.residues, self.restraints = Atoms(), Residues(self), Restraints()
        self.outputs = []
        self.loop = SimpleNamespace(outputs=self.outputs)
        self.disulfides = []
        self.reference_templates = []

    def residue_range(self, first, last):
        return int(first.split(':')[0]), int(last.split(':')[0])

    def patch(self, residue_type, residues):
        if residue_type != 'DISU':
            raise ValueError(f'Unsupported REEC patch: {residue_type}')
        self.disulfides.append(tuple(residue.intnum for residue in residues))

    def template_paths(self):
        result = []
        for name in self.knowns:
            matches = [Path(directory) / (str(name) + '.pdb') for directory in self.env.io.atom_files_directory]
            path = next((path.resolve() for path in matches if path.exists()), None)
            if path is None:
                raise ValueError(f'Template PDB not found: {name}')
            result.append((name, path))
        return result

    def environment_records(self):
        if not self.env.io.hetatm:
            return []
        return [line for line in self.template_paths()[0][1].read_text().splitlines()
                if line.startswith('HETATM')]

    def contact_lines(self):
        lines = [f'disulfide {a} {b}' for a, b in sorted(set(self.disulfides))]
        for item in self.restraints.items:
            if isinstance(item, Gaussian):
                kind, atoms = item.feature
                fields = ' '.join(f'{atom.residue} {atom.name}' for atom in atoms)
                if item.stdev <= 0:
                    raise ValueError('Restraint standard deviation must be positive')
                if kind == 'distance':
                    low = max(.01, item.mean-item.stdev)
                    lines.append(f'contact soft {fields} {low} {item.mean+item.stdev} {1/item.stdev**2}')
                elif kind == 'dihedral':
                    lines.append(f'dihedral {fields} {item.mean} {item.stdev}')
                else:
                    raise ValueError(f'Unsupported restraint feature: {kind}')
            elif item[0] == 'alpha':
                start, end = item[1]
                for residue in range(start, end-3):
                    lines.append(f'contact soft {residue} CA {residue+4} CA 5.8 6.4 25')
                    lines.append(f'contact soft {residue} O {residue+4} N 2.7 3.2 25')
            else:
                raise ValueError(f'Unsupported REEC restraint: {item}')
        if hasattr(self, 'select_loop_atoms'):
            for start, end in self.select_loop_atoms():
                if start <= end:
                    lines.append(f'mobile {start} {end}')
        environment = {int(line[22:26]) for line in self.environment_records()}
        for start, end in self.restraints.rigid_bodies:
            if any(residue not in environment for residue in range(start, end+1)):
                raise ValueError('REEC rigid bodies currently support ligand residues only')
        return lines

    def make(self):
        self.restraints = Restraints()
        self.disulfides = []
        if hasattr(self, 'special_patches'):
            self.special_patches(None)
        if hasattr(self, 'special_restraints'):
            self.special_restraints(None)
        args = self.env.options
        contacts = self.contact_lines()
        prepare(args)
        with reec.workspace(f'{self.sequence}-engine-', args.keep_work) as work:
            rows = pir_rows(self.alnfile)
            templates = [f'{name}={path}' if name in rows else str(path) for name, path in self.template_paths()]
            templates.extend(str(Path(path).resolve()) for path in self.reference_templates)
            records = self.environment_records()
            environment = None
            if records:
                environment = work / 'environment.pdb'
                # Retain the protein frame for ligand alignment; only HETATM records become fixed atoms.
                environment.write_text('\n'.join(line for line in self.template_paths()[0][1].read_text().splitlines()
                                                  if line.startswith(('ATOM  ', 'HETATM', 'CONECT', 'ENDMDL'))) + '\n')
            anneal_residues = args.anneal_residues
            if anneal_residues is None and self.sequence in reec.pu.prots:
                target_size = len(rows[self.sequence][1].replace('-', '').replace('.', ''))
                pocket = {reec.pu.resno_from_bw(self.sequence, bw) for bw in reec.pu.bsrs}
                anneal_residues = ','.join(str(number) for number in sorted(
                    number for number in pocket if number and 1 <= number <= target_size))
            engine_args = SimpleNamespace(iterations=self.max_var_iterations,
                models=self.ending_model-self.starting_model+1, contacts=args.contacts, chain=args.chain,
                candidate_prefix=work / 'candidate', seed=args.seed,
                anneal_steps=args.anneal_steps, anneal_temperature=args.anneal_temperature,
                anneal_residues=anneal_residues, clash_scale=getattr(args, 'clash_scale', 1.),
                retain_rejected=getattr(args, 'retain_rejected', False))
            alignment = protein_alignment(self.alnfile, self.knowns, self.sequence, work / 'protein.ali')
            try:
                output = reec.model(self.sequence, engine_args, work, templates, alignment,
                                    contacts, environment, bool(environment), bool(environment))
            except reec.CandidateRejected:
                if getattr(engine_args, 'last_rejected_model', None):
                    args.last_rejected_model = engine_args.last_rejected_model
                raise
            report = json.loads(Path(str(output)+'.json').read_text())
            if not report.get('valid'):
                raise ValueError('REEC returned an invalid model')
            candidates = []
            for index in range(engine_args.models):
                source = Path(f'{engine_args.candidate_prefix}.{index}.pdb')
                candidate_report = json.loads(Path(str(source)+'.json').read_text())
                if not candidate_report.get('valid'):
                    continue
                destination = Path(f'{self.sequence}.reec.{self.starting_model+index}.pdb').resolve()
                shutil.copyfile(source, destination)
                shutil.copyfile(str(source)+'.json', str(destination)+'.json')
                candidates.append({'name': destination.name, 'molpdf': candidate_report['final_energy'], 'failure': None})
        self.outputs[:] = candidates


def construct(model_class, env, *, references=(), **kwargs):
    source = Path(kwargs['alnfile']).resolve()
    rows = pir_rows(source)
    knowns = kwargs['knowns']
    names = ((knowns,) if isinstance(knowns, str) else tuple(knowns)) + (kwargs['sequence'],)
    if all(name in rows for name in names):
        width = max(len(rows[name][1]) for name in names)
        normalized = source.with_name(source.stem+'.model.ali')
        normalized.write_text(''.join(f'>P1;{name}\n{rows[name][0]}\n{rows[name][1].ljust(width, "-")}*\n'
                                      for name in names))
        kwargs['alnfile'] = str(normalized)
    model = model_class(env, **kwargs)
    if isinstance(model, ReecModel):
        model.reference_templates = [str(Path(path).resolve()) for path in references]
    return model


@contextmanager
def publication_guard(destination):
    """Restore the previous model if shared postprocessing or validation fails."""
    paths = [Path(destination), Path(str(destination)+'.json')]
    previous = {path: path.read_bytes() if path.exists() else None for path in paths}
    try:
        yield
    except BaseException:
        for path, content in previous.items():
            if content is None:
                path.unlink(missing_ok=True)
            else:
                path.write_bytes(content)
        raise


def finalize(model, destination, candidate=None, *, target=None):
    from reec_review import inspect, inspect_sterics
    destination = Path(destination)
    if isinstance(model, ReecModel):
        report = json.loads(Path(str(Path(candidate or model.outputs[0]['name']).resolve())+'.json').read_text())
    elif model is None:
        source_report = Path(str(candidate)+'.json')
        report = (json.loads(source_report.read_text()) if source_report.is_file()
                  else {'engine': 'validated', 'contacts': []})
    else:
        report = {'engine': 'MODELLER', 'contacts': []}
    if target is not None:
        report['target_sequence_length'] = len(reec.pu.prots[target]['sequence'])
    atoms = {(int(line[22:26]), line[12:16].strip()): tuple(float(line[i:i+8]) for i in (30,38,46))
             for line in destination.read_text().splitlines() if line.startswith(('ATOM  ', 'HETATM'))}
    for contact in report['contacts']:
        distance = math.dist(atoms[contact['res1'],contact['atom1']], atoms[contact['res2'],contact['atom2']])
        contact['distance'] = distance
        contact['violation'] = max(0,contact['min']-distance,distance-contact['max'])
    path = Path(str(destination)+'.json')
    report.setdefault('native_severe_clashes', report.get('severe_clashes'))
    report['severe_clashes'] = None
    report.pop('serialized_sterics', None)
    path.write_text(json.dumps(report, indent=2)+'\n')
    review = inspect(destination)
    audit_failed = False
    report['serialized_geometry'] = {key: value for key, value in review.items()
                                     if key not in ('model', 'optimizer_report')}
    if review['valid']:
        try:
            sterics = inspect_sterics(destination)
        except (OSError, ValueError, KeyError, subprocess.CalledProcessError) as exc:
            audit_failed = True
            review['valid'] = False
            review['failures'].append(f'Final steric audit could not complete: {type(exc).__name__}')
        else:
            report['serialized_sterics'] = sterics
            report['severe_clashes'] = sterics['severe_clashes']
            if sterics['severe_clashes']:
                review['valid'] = False
                review['failures'].append(f'Final serialized model has {sterics["severe_clashes"]} severe clashes')
    if not review['valid'] and not audit_failed and isinstance(model, ReecModel) and target is not None:
        from dohm_recovery import repair_existing, copy_pair
        repair_options = SimpleNamespace(iterations=6000, seed=report.get('seed', 0))
        repaired, work = repair_existing(target, destination, repair_options, reec.ROOT,
                                         contact_report=report['contacts'], orient=False)
        copy_pair(repaired, destination)
        repaired_report = json.loads(path.read_text())
        repaired_report['postprocessing_repair'] = {
            'work': str(work), 'failures': review['failures'],
            'original_serialized_sterics': report.get('serialized_sterics'),
            'original_native_severe_clashes': report.get('native_severe_clashes', report.get('severe_clashes'))}
        path.write_text(json.dumps(repaired_report, indent=2)+'\n')
        return
    if not review['valid']:
        report['valid'] = False
        report['serialized_geometry'] = {key: value for key, value in review.items()
                                         if key not in ('model', 'optimizer_report')}
        failure = reec.ROOT / 'tmp' / f'{target or destination.stem}.postprocessing-failure.json'
        failure.parent.mkdir(parents=True, exist_ok=True)
        failure.write_text(json.dumps(report, indent=2)+'\n')
        raise reec.CandidateRejected(f'Postprocessed model failed geometry review: {review["failures"]}', failure)
    report['valid'] = True
    report['serialized_geometry'] = {key: value for key, value in review.items()
                                     if key not in ('model', 'optimizer_report')}
    path.write_text(json.dumps(report, indent=2)+'\n')
