#!/usr/bin/env python3
"""Prepare constrained paired receptor models for isolated docking measurements."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

from reec_baseline import atoms, digest, pu
import reec_backend
from modeling_backend import finalize
import reec_review


def molecular_fingerprint(path):
    records = []
    for line in path.read_text().splitlines():
        if line.startswith(('ATOM  ', 'HETATM')):
            records.append(line[12:54] + line[76:78])
        elif line.startswith('CONECT'):
            records.append(line)
    return hashlib.sha256('\n'.join(records).encode()).hexdigest()


def audit(workspace, path):
    result = subprocess.run([str(workspace / 'bin/reec_steric_audit'), str(path)], cwd=workspace,
                            capture_output=True, text=True, check=True)
    return json.loads(result.stdout)


def annotate(workspace, receptor, candidate, output, clock):
    family = pu.family_from_protid(receptor)
    script = output.with_suffix('.phew')
    script.write_text(f'LOAD "pdbs/{family}/{receptor}.inactive.pdb" A I\n'
                      f'LOAD "{candidate}" A A\nBWCOPY I A\nUNCHAIN I\nSTRAND A\n'
                      f'UPRIGHT\nBWCENTER\nHYDRO\nSAVE "{output}"\n')
    environment = dict(os.environ, AROMA_BASELINE_TIME=str(clock),
                       LD_PRELOAD=str(workspace / 'bin/baseline_clock.so'), OMP_NUM_THREADS='1')
    with output.with_suffix('.annotation.log').open('w') as log:
        subprocess.run([str(workspace / 'bin/phew'), str(script)], cwd=workspace,
                       env=environment, stdout=log, stderr=subprocess.STDOUT, check=True)
    finalize(None, output, candidate, target=receptor)
    sterics = audit(workspace, output)
    Path(str(output)+'.sterics.json').write_text(json.dumps(sterics, indent=2) + '\n')
    if sterics['severe_clashes']:
        raise reec_backend.CandidateRejected('Annotated model has severe clashes')
    if sterics['unmatched_constructed_atoms']:
        raise reec_backend.CandidateRejected('Annotated model omits constructed heavy atoms')
    return sterics


def prepare_pair(workspace, recipe, models, max_attempts):
    receptor = recipe['receptor']
    folder = workspace / 'results/docking-models' / receptor
    folder.mkdir(parents=True, exist_ok=True)
    checkpoint = folder / 'result.json'
    if checkpoint.exists():
        previous = json.loads(checkpoint.read_text())
        if previous.get('status') == 'ready':
            return previous
    record = {'receptor': receptor, 'status': 'running', 'recipe': recipe, 'attempts': []}
    for attempt in range(max_attempts):
        iterations, clash = min(100000, 3000 * 2**attempt), min(8, 2**attempt)
        stage = {'seed': attempt, 'iterations': iterations, 'clash_scale': clash, 'conditions': {}}
        record['attempts'].append(stage)
        selected = {}
        for steps in (0, 64):
            condition = f'anneal-{steps}'
            work = folder / f'attempt-{attempt+1}' / condition
            work.mkdir(parents=True, exist_ok=True)
            output = work / 'result.pdb'
            command = [str(workspace / 'bin/reec'), '--alignment', recipe['alignment'],
                       '--target', receptor, '--output', str(output), '--contacts', recipe['contacts'],
                       '--chain', recipe['chain'], '--iterations', str(iterations), '--models', str(models),
                       '--seed', str(attempt), '--clash-scale', str(clash), '--anneal-steps', str(steps),
                       '--anneal-temperature', '12', '--candidate-prefix', str(work / 'candidate')]
            for template in recipe['templates']:
                command.extend(['--template', template])
            if recipe['anneal_residues']:
                command.extend(['--anneal-residues', recipe['anneal_residues']])
            if recipe['environment']:
                command.extend(['--environment', recipe['environment']])
                if recipe['fit_environment']:
                    command.append('--fit-environment')
                if recipe['environment_hydrogens']:
                    command.append('--environment-hydrogens')
            print(f'{receptor} attempt {attempt+1}, {condition}, {iterations} iterations', flush=True)
            started = time.perf_counter()
            with (work / 'engine.log').open('w') as log:
                result = subprocess.run(command, cwd=workspace, stdout=log, stderr=subprocess.STDOUT)
            details = {'command': command, 'exit_code': result.returncode,
                       'wall_seconds': time.perf_counter()-started, 'annotation_rejections': []}
            stage['conditions'][condition] = details
            if result.returncode not in (0, 2):
                raise RuntimeError(f'{receptor}: engine failed; inspect {work / "engine.log"}')
            candidates = []
            signatures = []
            for index in range(models):
                path = work / f'candidate.{index}.pdb'
                report = json.loads(Path(str(path)+'.json').read_text())
                signatures.append(report['start_signature'])
                if report['valid'] and reec_review.inspect(path)['valid']:
                    candidates.append((report['final_energy'], path))
            details['start_signatures'] = signatures
            details['accepted_native_candidates'] = len(candidates)
            for score, candidate in sorted(candidates):
                annotated = work / f'{receptor}.candidate-{candidate.stem.split(".")[1]}.pdb'
                try:
                    if steps == 64 and 0 in selected and molecular_fingerprint(candidate) == selected[0]['raw_fingerprint']:
                        shutil.copy2(selected[0]['model'], annotated)
                        finalize(None, annotated, candidate, target=receptor)
                        sterics = selected[0]['sterics']
                    else:
                        sterics = annotate(workspace, receptor, candidate, annotated, 1770000000)
                except reec_backend.CandidateRejected as error:
                    details['annotation_rejections'].append({'candidate': str(candidate), 'reason': str(error)})
                    continue
                selected[steps] = {'model': str(annotated), 'candidate': str(candidate), 'score': score,
                                   'sterics': sterics, 'raw_fingerprint': molecular_fingerprint(candidate),
                                   'model_sha256': digest(annotated),
                                   'model_fingerprint': molecular_fingerprint(annotated)}
                break
            checkpoint.write_text(json.dumps(record, indent=2) + '\n')
        if stage['conditions']['anneal-0']['start_signatures'] != stage['conditions']['anneal-64']['start_signatures']:
            raise RuntimeError(f'{receptor}: paired starting coordinates differ')
        if len(selected) == 2:
            record['status'] = 'ready'
            record['selected'] = {str(steps): value for steps, value in selected.items()}
            record['identical_models'] = selected[0]['model_fingerprint'] == selected[64]['model_fingerprint']
            checkpoint.write_text(json.dumps(record, indent=2) + '\n')
            print(f'{receptor}: paired models ready', flush=True)
            return record
    record['status'] = 'no_valid_pair'
    checkpoint.write_text(json.dumps(record, indent=2) + '\n')
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=3)
    parser.add_argument('--models', type=int, default=10)
    parser.add_argument('--max-attempts', type=int, default=5)
    args = parser.parse_args()
    workspace = args.workspace.resolve()
    reec_backend.ROOT, reec_backend.HM = workspace, workspace / 'hm'
    pu.load_prots()
    recipes = [json.loads(path.read_text()) for path in sorted((workspace / 'results/recipes').glob('*/recipe.json'))]
    if len(recipes) != 10:
        parser.error('Expected ten frozen receptor recipes')
    output = workspace / 'results/docking-models/manifest.json'
    output.parent.mkdir(parents=True, exist_ok=True)
    manifest = {'scope': 'Workflow contacts, followed by orientation and hydrogenation; no pocket repointing or IC minimization',
                'annotation_clock': 1770000000, 'models_per_batch': args.models,
                'max_paired_attempts': args.max_attempts, 'workers': args.workers, 'cases': {}}
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        pending = {pool.submit(prepare_pair, workspace, recipe, args.models, args.max_attempts): recipe['receptor'] for recipe in recipes}
        for future in as_completed(pending):
            receptor = pending[future]
            try:
                manifest['cases'][receptor] = future.result()
            except Exception as error:
                manifest['cases'][receptor] = {'status': 'error', 'error': str(error)}
                print(f'{receptor}: {error}', file=sys.stderr, flush=True)
            output.write_text(json.dumps(manifest, indent=2) + '\n')
    return 0 if all(case['status'] == 'ready' for case in manifest['cases'].values()) else 2


if __name__ == '__main__':
    raise SystemExit(main())
