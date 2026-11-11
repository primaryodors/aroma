#!/usr/bin/env python3
"""Run isolated, paired native REEC measurements without publishing models."""
import argparse
import hashlib
import itertools
import json
import os
import platform
import shutil
import statistics
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from modeling_backend import pir_rows
from dohm import experimental_template
from reec_backend import ROOT, pu
import reec_review

CODES = dict(zip(
    'ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR'.split(),
    'ACDEFGHIKLMNPQRSTVWY'))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atoms(path, chain='A'):
    return {(int(s[22:26]), s[12:16].strip()): [float(s[i:i+8]) for i in (30, 38, 46)]
            for s in path.read_text().splitlines()
            if s.startswith('ATOM  ') and s[16] in (' ', 'A') and s[21] == chain
            and not s[12:16].strip().lstrip('0123456789').startswith('H')}


def reference_atoms(name, template_id):
    """Map an experimental chain onto target numbering through PIR columns."""
    description, row = pir_rows(ROOT / 'hm/experimental.ali')[template_id]
    fields = description.split(':')
    path = ROOT / 'hm/tpl' / f'{template_id}.pdb'
    chain = fields[3].strip()
    points = atoms(path, chain)
    residue_numbers = list(dict.fromkeys(int(s[22:26]) for s in path.read_text().splitlines()
                                        if s.startswith('ATOM  ') and s[12:16].strip() == 'CA'
                                        and s[21] == chain and s[16] in (' ', 'A')
                                        and int(s[22:26]) >= int(fields[2])))
    if len(residue_numbers) != len(row.replace('-', '')):
        raise ValueError('Experimental reference residue count differs from its PIR row')
    sequence = ''.join(CODES[s[17:20]] for s in path.read_text().splitlines()
                       if s.startswith('ATOM  ') and s[12:16].strip() == 'CA'
                       and s[21] == chain and s[16] in (' ', 'A')
                       and int(s[22:26]) >= int(fields[2]))
    if sequence != row.replace('-', ''):
        raise ValueError('Experimental reference sequence differs from its PIR row')
    aligned = ''.join(pu.prots[name]['aligned'].split())
    mapping, target_no, cursor = {}, 0, 0
    for template_aa, target_aa in itertools.zip_longest(row, aligned, fillvalue='-'):
        if target_aa != '-':
            target_no += 1
        if template_aa != '-':
            number = residue_numbers[cursor]
            cursor += 1
            if target_aa == template_aa:
                mapping[number] = target_no
    mapped = {(mapping[number], atom): point for (number, atom), point in points.items() if number in mapping}
    return mapped, {'id': template_id, 'source': str(path.relative_to(ROOT)),
                    'sha256': digest(path), 'matched_residues': len(mapping)}


def pocket_distance(first, second, pocket, atom_names=None):
    """Superpose common CA atoms, then compare common pocket heavy atoms."""
    common = sorted(first.keys() & second.keys())
    ca = [key for key in common if key[1] == 'CA']
    keys = [key for key in common if key[0] in pocket and (atom_names is None or key[1] in atom_names)]
    if len(ca) < 3 or not keys:
        return None
    a, b = (np.array([model[key] for key in ca]) for model in (first, second))
    ac, bc = a.mean(axis=0), b.mean(axis=0)
    u, _, vt = np.linalg.svd((a-ac).T @ (b-bc))
    correction = np.eye(3)
    correction[2, 2] = np.linalg.det(u @ vt)
    rotation = u @ correction @ vt
    delta = (np.array([first[key] for key in keys])-ac) @ rotation + bc
    delta -= np.array([second[key] for key in keys])
    return {'rmsd_angstrom': float(np.sqrt(np.mean(np.sum(delta**2, axis=1)))),
            'common_pocket_atoms': len(keys), 'alignment_ca_atoms': len(ca)}


def reference_accuracy(reference, coordinates, pocket):
    return {'pocket_heavy': pocket_distance(reference, coordinates, pocket),
            'pocket_backbone': pocket_distance(reference, coordinates, pocket, {'N', 'CA', 'C', 'O'}),
            'whole_ca': pocket_distance(reference, coordinates, {key[0] for key in reference}, {'CA'})}


def target_row(name, extend_terminal_tail=False):
    target = ''.join(pu.prots[name]['aligned'].split())
    sequence = pu.prots[name]['sequence']
    core = target.replace('-', '')
    if core == sequence:
        return target
    if extend_terminal_tail and core and sequence.startswith(core):
        return (target.rstrip('-') + sequence[len(core):]).ljust(len(target), '-')
    raise ValueError(f'Target alignment mismatch: {name}')


def prepare_case(name, folder, template_ids=None, extend_terminal_tail=False):
    folder.mkdir()
    experimental = pir_rows(ROOT / 'hm/experimental.ali')
    ids = set(template_ids if template_ids is not None else pu.templates_for_hm(name))
    if name == 'OR51E1' and template_ids is None:
        ids |= {'8uxv', '8hti', '9wg4', '8uxy', '8uy0', '8uyq', '9wpm',
                '9ldv', '9ldw', '9ldx', '9ldz'}
    rows, templates, inputs = {}, [], []
    for template_id, (description, row) in experimental.items():
        if template_id not in ids or any(tag in description.lower() for tag in ('inactive', 'apo ')):
            continue
        source = ROOT / 'hm/tpl' / f'{template_id}.pdb'
        if not source.is_file():
            raise FileNotFoundError(source)
        path, description = experimental_template(template_id, source, description, folder)
        rows[template_id] = description, row
        templates.extend(['--template', f'{template_id}={path}'])
        inputs.append({'source': str(source.relative_to(ROOT)), 'sha256': digest(source)})
    if not rows:
        raise ValueError(f'No experimental templates available for {name}')
    if template_ids is not None and ids != set(rows):
        raise ValueError(f'Explicit template selection is incomplete for {name}')
    target = target_row(name, extend_terminal_tail)
    rows[name] = f'sequence:{name}:1:A::::::', target
    width = max(len(row) for _, row in rows.values())
    alignment = folder / 'protein.ali'
    alignment.write_text(''.join(f'>P1;{key}\n{desc}\n{row.ljust(width, "-")}*\n'
                                 for key, (desc, row) in rows.items()))
    pocket = sorted({int(number) for bw in pu.bsrs
                     if (number := pu.resno_from_bw(name, bw)) and 1 <= int(number) <= len(pu.prots[name]['sequence'])})
    return alignment, templates, inputs, pocket


def prepare_repair(folder, source, template, reference):
    """Snapshot a retained OR1B1 repair without altering its source files."""
    folder.mkdir()
    inputs = []
    paths = {}
    for key, path in {'alignment': source / 'protein.ali',
                      'contacts': source / 'receptor.contacts',
                      'environment': source / 'environment.pdb',
                      'template': template, 'reference': reference}.items():
        target = folder / (key + path.suffix)
        shutil.copyfile(path, target)
        paths[key] = target
        inputs.append({'role': key, 'source': str(path), 'sha256': digest(target)})
    rows = pir_rows(paths['alignment'])
    if 'OR1B1_tpl' not in rows or 'OR1B1' not in rows:
        raise ValueError('Repair alignment must contain OR1B1_tpl and OR1B1')
    pocket = sorted({int(number) for bw in pu.bsrs
                     if (number := pu.resno_from_bw('OR1B1', bw)) and 1 <= int(number) <= 317})
    extra = ['--template', f'OR1B1_tpl={paths["template"]}', '--template', str(paths['reference']),
             '--contacts', str(paths['contacts']), '--environment', str(paths['environment']),
             '--fit-environment', '--environment-hydrogens']
    return paths['alignment'], extra, inputs, pocket


def summarize(runs, models, pocket):
    comparisons = []
    for (label_a, a), (label_b, b) in itertools.combinations(models, 2):
        comparison = pocket_distance(a, b, pocket)
        if comparison:
            comparisons.append(dict(first=label_a, second=label_b, **comparison))
    groups = []
    for steps in sorted({run['anneal_steps'] for run in runs}):
        selected = [run for run in runs if run['anneal_steps'] == steps]
        candidates = [candidate for run in selected for candidate in run['candidates']]
        groups.append({'anneal_steps': steps, 'runs': len(selected),
                       'accepted_runs': sum(run['exit_code'] == 0 for run in selected),
                       'candidates': len(candidates),
                       'accepted_candidates': sum(candidate['valid'] for candidate in candidates),
                       'independently_accepted_candidates': sum(
                           candidate.get('independent_review', {}).get('valid', False)
                           for candidate in candidates),
                       'median_wall_seconds': statistics.median(run['wall_seconds'] for run in selected)})
        for metric in ('severe_clashes', 'max_bond_error', 'max_geometry_distance_error',
                       'max_planarity_volume_error', 'chirality_errors', 'mandatory_violations',
                       'final_energy'):
            values = [candidate[metric] for candidate in candidates if metric in candidate]
            if values:
                groups[-1][metric] = {'min': min(values), 'median': statistics.median(values),
                                      'max': max(values)}
    return {'groups': groups, 'pocket_comparisons': comparisons}


def paired_summary(runs):
    """Compare matched candidates, separating final validity from score changes."""
    by_seed = {}
    for run in runs:
        settings = by_seed.setdefault(run['seed'], {})
        if run['anneal_steps'] in settings:
            raise ValueError('Duplicate seed and annealing setting')
        settings[run['anneal_steps']] = run
    outcomes, best_scores, time_ratios = [], [], []

    def accepted(candidate):
        return candidate['valid'] and candidate.get('independent_review', {}).get('valid', False)

    for seed, settings in by_seed.items():
        if not all(steps in settings for steps in (0, 64)):
            continue
        off, on = settings[0], settings[64]
        if len(off['candidates']) != len(on['candidates']):
            raise ValueError('Paired candidate counts differ')
        time_ratios.append(on['wall_seconds'] / off['wall_seconds'])
        for a, b in zip(off['candidates'], on['candidates']):
            if a['candidate_index'] != b['candidate_index'] or a['start_signature'] != b['start_signature']:
                raise ValueError('Paired starts differ')
            comparable = (accepted(a) and accepted(b)
                          and a.get('released_template_residues', []) == b.get('released_template_residues', []))
            outcomes.append({'seed': seed, 'candidate_index': a['candidate_index'],
                             'off_accepted': accepted(a), 'on_accepted': accepted(b),
                             'off_energy': a['final_energy'], 'on_energy': b['final_energy'],
                             'score_comparable': comparable,
                             'score_improvement_percent': (100 * (a['final_energy']-b['final_energy'])
                                 / max(abs(a['final_energy']), 1.)) if comparable else None})
        a = min((c for c in off['candidates'] if accepted(c)), key=lambda c: c['final_energy'], default=None)
        b = min((c for c in on['candidates'] if accepted(c)), key=lambda c: c['final_energy'], default=None)
        if a is not None and b is not None:
            comparable = a.get('released_template_residues', []) == b.get('released_template_residues', [])
            best_scores.append({'seed': seed, 'off_energy': a['final_energy'], 'on_energy': b['final_energy'],
                                'score_comparable': comparable,
                                'score_improvement_percent': 100*(a['final_energy']-b['final_energy'])
                                    / max(abs(a['final_energy']), 1.) if comparable else None})
    return {'matched_pairs': len(outcomes), 'starts_match': True,
            'off_accepted': sum(p['off_accepted'] for p in outcomes),
            'on_accepted': sum(p['on_accepted'] for p in outcomes),
            'rescued': sum(not p['off_accepted'] and p['on_accepted'] for p in outcomes),
            'lost': sum(p['off_accepted'] and not p['on_accepted'] for p in outcomes),
            'median_runtime_ratio': statistics.median(time_ratios) if time_ratios else None,
            'candidate_outcomes': outcomes, 'best_scores': best_scores}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--receptors', nargs='+', default=['OR1A1', 'OR1B1', 'OR51E1', 'OR51E2'])
    parser.add_argument('--seeds', nargs='+', type=int, default=[0, 1])
    parser.add_argument('--models', type=int, default=3)
    parser.add_argument('--iterations', type=int, default=3000)
    parser.add_argument('--repair-inputs', type=Path, help='Retained OR1B1 engine directory')
    parser.add_argument('--repair-template', type=Path)
    parser.add_argument('--reference-template', type=Path)
    parser.add_argument('--template-ids', nargs='+', help='Explicit experimental template IDs for isolated controls')
    parser.add_argument('--reference-id', help='Experimental template excluded from reconstruction and used for evaluation')
    parser.add_argument('--counterbalance', action='store_true', help='Alternate off/on execution order across seeds and receptors')
    parser.add_argument('--extend-terminal-tail', action='store_true', help='Restore a missing terminal sequence suffix in isolated target rows')
    args = parser.parse_args()
    if args.models < 1 or args.iterations < 1 or any(seed < 0 or seed > 2147483647 for seed in args.seeds):
        parser.error('Invalid model count, iteration budget, or seed')
    if args.repair_inputs:
        if not args.repair_template or not args.reference_template:
            parser.error('Repair measurements require --repair-template and --reference-template')
        args.receptors = ['OR1B1']
    if args.reference_id and (not args.template_ids or args.reference_id in args.template_ids or args.repair_inputs):
        parser.error('Reference evaluation requires explicit template IDs excluding the reference and no repair inputs')
    folder = args.output.resolve()
    folder.mkdir(parents=True, exist_ok=False)
    with (folder / 'build.log').open('w') as log:
        subprocess.run(['make', 'bin/reec'], cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, check=True)
    pu.load_prots()
    report = {'scope': 'Native engine only; no PHEW publication, retries, or docking',
              'accuracy_status': 'Not evaluated; template coordinates are not held-out references',
              'host': platform.node(), 'platform': platform.platform(),
              'python': platform.python_version(), 'numpy': np.__version__,
              'cpu': platform.processor(), 'logical_cpu_count': os.cpu_count(),
              'engine_sha256': digest(ROOT / 'bin/reec'),
              'source_sha256': digest(ROOT / 'src/reec.cpp'),
              'database_sha256': digest(ROOT / 'data/receptor.json'),
              'experimental_alignment_sha256': digest(ROOT / 'hm/experimental.ali'),
              'settings': {'seeds': args.seeds, 'models': args.models, 'iterations': args.iterations,
                           'anneal_steps': [0, 64], 'anneal_temperature': 12, 'clash_scale': 1,
                           'counterbalanced_order': args.counterbalance,
                           'extend_terminal_tail': args.extend_terminal_tail},
              'limitations': ['Experimental-template controls do not include validated subfamily models',
                             'OR1B1 is generation, not replay of the absent maltol repair inputs',
                             'OR51E2 is full-sequence generation, not truncated-source publication',
                             'Workflow contacts and docking environments are not included',
                             'Pocket RMSD uses atom names; symmetric side-chain atoms are not remapped',
                             'Pocket RMSD measures diversity, not experimental accuracy'],
              'cases': {}}

    def save():
        (folder / 'summary.json').write_text(json.dumps(report, indent=2) + '\n')

    save()
    if args.repair_inputs:
        report['limitations'] = [
            'Native repair only; no PHEW publication, retries, or docking',
            'Reference-template identity relative to the historical run is unverified',
            'Pocket RMSD measures diversity, not experimental accuracy',
            'Symmetric side-chain atom names are not remapped']
    for case_index, name in enumerate(args.receptors):
        case_folder = folder / name
        if args.repair_inputs:
            alignment, templates, inputs, pocket = prepare_repair(
                case_folder, args.repair_inputs.resolve(), args.repair_template.resolve(),
                args.reference_template.resolve())
        else:
            alignment, templates, inputs, pocket = prepare_case(
                name, case_folder, args.template_ids, args.extend_terminal_tail)
        case = {'inputs': inputs, 'alignment_sha256': digest(alignment),
                'pocket_residues': pocket, 'runs': []}
        if args.extend_terminal_tail:
            case['terminal_extension_residues'] = len(pu.prots[name]['sequence']) - len(
                ''.join(pu.prots[name]['aligned'].split()).replace('-', ''))
        reference = None
        if args.reference_id:
            reference, case['experimental_reference'] = reference_atoms(name, args.reference_id)
            report['accuracy_status'] = 'Experimental reference excluded from reconstruction inputs; native output only'
        report['cases'][name] = case
        models = []
        for seed_index, seed in enumerate(args.seeds):
            order = (64, 0) if args.counterbalance and (case_index + seed_index) % 2 else (0, 64)
            for steps in order:
                label = f'seed-{seed}-anneal-{steps}'
                work = case_folder / label
                work.mkdir()
                output = work / 'result.pdb'
                command = [str(ROOT / 'bin/reec'), '--alignment', str(alignment), '--target', name,
                           '--output', str(output), '--models', str(args.models),
                           '--iterations', str(args.iterations), '--seed', str(seed),
                           '--anneal-steps', str(steps), '--anneal-temperature', '12',
                           '--chain', 'A', '--clash-scale', '1',
                           '--candidate-prefix', str(work / 'candidate'), *templates]
                if pocket:
                    command.extend(['--anneal-residues', ','.join(map(str, pocket))])
                print(f'{name} {label}: starting', flush=True)
                started = time.perf_counter()
                with (work / 'engine.log').open('w') as log:
                    result = subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
                run = {'label': label, 'seed': seed, 'anneal_steps': steps, 'command': command,
                       'wall_seconds': time.perf_counter()-started, 'exit_code': result.returncode,
                       'candidates': []}
                case['runs'].append(run)
                if result.returncode not in (0, 2):
                    save()
                    raise RuntimeError(f'Unexpected engine exit; inspect {work / "engine.log"}')
                for index in range(args.models):
                    path = work / f'candidate.{index}.pdb'
                    candidate = json.loads(Path(str(path)+'.json').read_text())
                    if candidate['valid']:
                        review = reec_review.inspect(path)
                        candidate['independent_review'] = review
                        if review['valid']:
                            coordinates = atoms(path)
                            models.append((f'{label}/candidate.{index}', coordinates))
                            if reference is not None:
                                candidate['experimental_pocket_accuracy'] = pocket_distance(reference, coordinates, pocket)
                                candidate['experimental_accuracy'] = reference_accuracy(reference, coordinates, pocket)
                    run['candidates'].append(candidate)
                print(f'{name} {label}: exit {result.returncode}, {run["wall_seconds"]:.2f}s', flush=True)
                case['summary'] = summarize(case['runs'], models, pocket)
                save()
    print(f'Report: {folder / "summary.json"}', flush=True)


if __name__ == '__main__':
    main()
