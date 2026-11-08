#!/usr/bin/env python3
"""Prepare standard docking configurations and run frozen baseline panels."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import re
import runpy
import shutil
import subprocess
import sys
import time

from reec_baseline import digest
from reec_workflow_baseline import ROOT, bind_workspace


def parse_dock(text):
    parsed = [float(value) for value in re.findall(r'^Total:\s+(\S+)', text, re.M)]
    energies = [value for value in parsed if math.isfinite(value)]
    nonfinite = len(energies) != len(parsed)
    found = re.findall(r'^(\d+) pose\(s\) found\.', text, re.M)
    disqualified = re.findall(r'^Disqualified because: (.*)', text, re.M)
    count = int(found[-1]) if found else 0
    return {'poses': count, 'first_energy': energies[0] if energies else None,
            'best_energy': min(energies) if energies else None, 'energies': energies,
            'disqualifications': disqualified, 'nonfinite_energy': nonfinite,
            'operational_pass': count >= 4 and bool(energies) and energies[0] < 0 and not disqualified and not nonfinite}


def prepare_case(workspace, case, condition, model, panel_name):
    from data import odorutils
    folder = workspace / 'results/dock-configs' / panel_name / f'{case["receptor"]}-{case["ligand_id"][:8]}' / str(condition)
    folder.mkdir(parents=True, exist_ok=True)
    checkpoint = folder / 'prepared.json'
    if checkpoint.exists():
        previous = json.loads(checkpoint.read_text())
        if previous['model_sha256'] == digest(Path(model)):
            return previous
    module, backend = bind_workspace(workspace, 'dohm')
    backend.pu.load_prots()
    odorutils.load_odors()
    for name in ('run.py', 'example.config'):
        destination = workspace / name
        if not destination.exists():
            destination.symlink_to(ROOT / name)
    (workspace / 'test').mkdir(exist_ok=True)
    source_tool = ROOT / 'test/moiety_test'
    workspace_tool = workspace / 'test/moiety_test'
    if source_tool.resolve() != workspace_tool.resolve():
        shutil.copy2(source_tool, workspace_tool)
    (workspace / 'tmp/nodelete').touch()
    family = backend.pu.family_from_protid(case['receptor'])
    active = workspace / 'pdbs' / family / f'{case["receptor"]}.active.pdb'
    cavity = active.with_suffix('.cvty')
    backups = {path: (path.read_bytes(), path.stat().st_atime_ns, path.stat().st_mtime_ns) if path.exists() else None
               for path in (active, cavity)}
    captured = []
    original_run, original_argv = subprocess.run, sys.argv
    cache_key = hashlib.sha256(json.dumps({
        'model': digest(Path(model)), 'ic': digest(workspace / 'bin/ic'),
        'cavity_search': digest(workspace / 'bin/cavity_search'),
        'cavopts': digest(workspace / 'data/cavopts.json'),
        'parameters': {p.name: digest(p) for p in sorted((workspace / 'data').glob('*.dat'))},
    }, sort_keys=True).encode()).hexdigest()
    model_cache = workspace / 'results/model-preparation' / cache_key
    model_cache.mkdir(parents=True, exist_ok=True)
    cached_cavity = model_cache / 'model.cvty'
    cached_contacts = model_cache / 'contact-pairs.json'
    if not cached_cavity.exists() or not cached_contacts.exists():
        for previous in (workspace / 'results/dock-configs').glob('*/*/*/prepared.json'):
            old = json.loads(previous.read_text())
            if old['model_sha256'] != digest(Path(model)):
                continue
            prior_cavity = Path(old['config']).parent / 'model.cvty'
            if prior_cavity.exists():
                shutil.copy2(prior_cavity, cached_cavity)
            pairs = [line.split()[1:] for line in Path(old['config']).read_text().splitlines()
                     if re.fullmatch(r'CNTCT \d+ \d+', line)]
            cached_contacts.write_text(json.dumps({'origin': 'numeric contacts in prior standard configuration',
                                                   'pairs': pairs}, indent=2) + '\n')
            break
    environment = dict(os.environ, LD_PRELOAD=str(workspace / 'bin/baseline_clock.so'),
                       AROMA_BASELINE_TIME='1770000000', OMP_NUM_THREADS='1')

    def intercept(command, *arguments, **keywords):
        executable = Path(str(command[0])).name if isinstance(command, (list, tuple)) else ''
        if executable == 'aromadock':
            captured.append((workspace / command[1]).resolve())
            return subprocess.CompletedProcess(command, 0)
        if executable == 'ic' and 'nooil' in command and cached_contacts.exists():
            pairs = json.loads(cached_contacts.read_text())['pairs']
            text = ''.join(f'{a}-{b}\n' for a, b in pairs).encode()
            return subprocess.CompletedProcess(command, 0, stdout=text)
        if executable == 'cavity_search' and cached_cavity.exists():
            destination = Path(command[command.index('-o')+1])
            shutil.copy2(cached_cavity, destination)
            return subprocess.CompletedProcess(command, 0)
        keywords.setdefault('env', environment)
        result = original_run(command, *arguments, **keywords)
        if executable == 'cavity_search' and cavity.exists():
            shutil.copy2(cavity, cached_cavity)
        elif executable == 'ic' and 'nooil' in command and result.stdout:
            pairs = []
            for line in result.stdout.decode().splitlines():
                value = re.sub(r'\([0-9.]*\)[.][A-Z0-9]+', '', line.split(':')[0])
                value = re.sub('[^0-9-]', '', value)
                if re.fullmatch(r'\d+-\d+', value):
                    pairs.append(value.split('-'))
            cached_contacts.write_text(json.dumps({'origin': 'standard IC contact extraction', 'pairs': pairs}, indent=2)+'\n')
        return result

    try:
        shutil.copy2(model, active)
        cavity.unlink(missing_ok=True)
        subprocess.run = intercept
        sys.argv = [str(workspace / 'run.py'), case['receptor'], case['ligand_id'], 'noi', 'ns']
        runpy.run_path(str(workspace / 'run.py'), run_name='__main__')
        if len(captured) != 1:
            raise RuntimeError('Expected one active-state docking configuration')
        text = captured[0].read_text()
        if not cavity.exists():
            raise RuntimeError('Docking preparation did not produce a cavity file')
        copied_cavity = folder / 'model.cvty'
        shutil.copy2(cavity, copied_cavity)
        lines, sdfs = [], []
        for line in text.splitlines():
            fields = line.split()
            if not fields:
                continue
            key = fields[0]
            if key == 'PROT':
                line = f'PROT {model}'
            elif key == 'VCVTY':
                line = f'VCVTY {copied_cavity}'
            elif key in ('LIG', 'ISO', 'FORM'):
                source = Path(fields[1])
                if not source.is_absolute():
                    source = workspace / source
                target = folder / source.name
                shutil.copy2(source, target)
                sdfs.append({'path': str(target), 'sha256': digest(target), 'role': key})
                line = f'{key} {target}'
            lines.append(line)
        if case['receptor'] == 'OR2T11' and not any(line.startswith('MCOORD Cu +1 ') for line in lines):
            lines.append('MCOORD Cu +1 C5.42 C5.43 M5.46')
        config = folder / 'base.config'
        config.write_text('\n'.join(lines) + '\n')
        result = {'case': case, 'condition': condition, 'config': str(config), 'model': model,
                  'model_sha256': digest(Path(model)), 'sdfs': sdfs,
                  'cavity_sha256': digest(copied_cavity), 'config_sha256': digest(config)}
        checkpoint.write_text(json.dumps(result, indent=2) + '\n')
        return result
    finally:
        subprocess.run, sys.argv = original_run, original_argv
        for path, previous in backups.items():
            if previous is None:
                path.unlink(missing_ok=True)
            else:
                content, accessed, modified = previous
                path.write_bytes(content)
                os.utime(path, ns=(accessed, modified))
        os.chdir(workspace)


def fingerprint(prepared, settings, engine):
    lines = []
    for line in Path(prepared['config']).read_text().splitlines():
        fields = line.split()
        if not fields or fields[0] in ('OUT', 'OUTPDB', 'PROGRESS'):
            continue
        if fields[0] == 'PROT':
            fields[1] = prepared['model_sha256']
        elif fields[0] in ('LIG', 'ISO', 'FORM', 'VCVTY'):
            fields[1] = digest(Path(fields[1]))
        lines.append(' '.join(fields))
    budget = {key: settings[key] for key in ('poses', 'iterations', 'threads')}
    content = json.dumps({'lines': lines, 'settings': budget, 'engine': engine}, sort_keys=True)
    return hashlib.sha256(content.encode()).hexdigest()


def run_job(workspace, prepared, settings, clock, engine):
    key = fingerprint(prepared, settings, engine) + '-' + str(clock)
    folder = workspace / 'results/docking-cache' / key
    folder.mkdir(parents=True, exist_ok=True)
    with (folder / 'job.lock').open('w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        checkpoint = folder / 'result.json'
        if checkpoint.exists():
            return dict(json.loads(checkpoint.read_text()), reused=True)
        for name in ('data', 'sdf', 'bin'):
            path = folder / name
            if not path.exists():
                path.symlink_to(workspace / name)
        (folder / 'tmp').mkdir(exist_ok=True)
        lines = []
        for line in Path(prepared['config']).read_text().splitlines():
            fields = line.split()
            if fields and fields[0] in ('OUT', 'OUTPDB', 'POSE', 'ITERS', 'ITER', 'THREADS', 'THREAD', 'T', 'RETRY'):
                continue
            lines.append(line)
        output = folder / 'result.dock'
        lines.extend([f'OUT {output}', f'OUTPDB {settings["poses"]} {folder}/pose%o.pdb',
                      f'POSE {settings["poses"]}', f'ITERS {settings["iterations"]}', 'THREADS 1', 'RETRY 0'])
        config = folder / 'job.config'
        config.write_text('\n'.join(lines) + '\n')
        environment = dict(os.environ, LD_PRELOAD=str(workspace / 'bin/baseline_clock.so'),
                           AROMA_BASELINE_TIME=str(clock), OMP_NUM_THREADS='1')
        command = [str(workspace / 'bin/aromadock'), str(config), '--ns', '-t', '1']
        started = time.perf_counter()
        print(f'Docking {prepared["case"]["receptor"]} {prepared["case"]["ligand"]}, anneal {prepared["condition"]}, seed {clock}', flush=True)
        with (folder / 'engine.log').open('w') as log:
            try:
                completed = subprocess.run(command, cwd=folder, env=environment, stdout=log,
                                           stderr=subprocess.STDOUT, timeout=1800)
                code = completed.returncode
            except subprocess.TimeoutExpired:
                code = 'timeout'
        result = {'command': command, 'exit_code': code, 'clock_seed': clock,
                  'wall_seconds': time.perf_counter()-started, 'folder': str(folder),
                  'engine_sha256': engine, 'config_sha256': digest(config),
                  'dock_sha256': digest(output) if output.exists() else None,
                  'metrics': parse_dock(output.read_text()) if output.exists() else None,
                  'pose_files': [str(path) for path in sorted(folder.glob('pose*.pdb'))]}
        checkpoint.write_text(json.dumps(result, indent=2) + '\n')
        return dict(result, reused=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['prepare', 'run'])
    parser.add_argument('--workspace', type=Path, required=True)
    parser.add_argument('--panel', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=3)
    args = parser.parse_args()
    workspace, panel_path = args.workspace.resolve(), args.panel.resolve()
    panel = json.loads(panel_path.read_text())
    name = panel_path.stem
    manifest_path = workspace / 'results' / f'{name}.prepared.json'
    if args.action == 'prepare':
        models = json.loads((workspace / 'results/docking-models/manifest.json').read_text())
        manifest = {'panel': str(panel_path), 'panel_sha256': digest(panel_path), 'cases': [], 'pending': []}
        for case in panel['cases']:
            model_case = models['cases'].get(case['receptor'], {})
            if model_case.get('status') != 'ready':
                manifest['pending'].append(case)
                continue
            for condition in case['conditions']:
                model = model_case['selected'][str(condition)]['model']
                manifest['cases'].append(prepare_case(workspace, case, condition, model, name))
                manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
        manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
        return 0
    manifest = json.loads(manifest_path.read_text())
    output = workspace / 'results' / f'{name}.docking.json'
    engine = digest(workspace / 'bin/aromadock')
    results = {'panel': str(panel_path), 'engine_sha256': engine, 'jobs': []}
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        pending = {}
        for prepared in manifest['cases']:
            for clock in panel['settings']['clock_seeds']:
                future = pool.submit(run_job, workspace, prepared, panel['settings'], clock, engine)
                pending[future] = prepared
        for future in as_completed(pending):
            prepared = pending[future]
            try:
                result = future.result()
            except Exception as error:
                result = {'error': str(error)}
            results['jobs'].append({'case': prepared['case'], 'condition': prepared['condition'], **result})
            output.write_text(json.dumps(results, indent=2) + '\n')
    return 0 if all(job.get('exit_code') == 0 for job in results['jobs']) else 2


if __name__ == '__main__':
    raise SystemExit(main())
