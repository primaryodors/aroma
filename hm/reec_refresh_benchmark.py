#!/usr/bin/env python3
"""Rebuild frozen models and docking jobs before a detached remote batch."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import signal
import shutil
import subprocess
import sys
import threading
import time

from reec_thermal import run_guarded


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def refresh_recipes(root):
    from dohm import experimental_ensemble_alignment
    from reec_backend import pu
    pu.load_prots()
    for path in sorted((root / 'results/recipes').glob('*/recipe.json')):
        recipe = json.loads(path.read_text())
        # The validated-subfamily recipe already retains its original donors.
        if recipe['receptor'] == 'OR51E1':
            continue
        work = path.parent / 'weighted-inputs'
        work.mkdir(exist_ok=True)
        if not (work / 'tpl').exists():
            (work / 'tpl').symlink_to(root / 'hm/tpl')
        shutil.copy2(root / 'hm/experimental.ali', work / 'experimental.ali')
        alignment, templates, provenance = experimental_ensemble_alignment(recipe['receptor'], work)
        recipe['previous_template_inputs'] = recipe['template_inputs']
        recipe['alignment'] = str(alignment)
        recipe['alignment_sha256'] = digest(alignment)
        recipe['templates'] = templates
        recipe['template_inputs'] = [{'source': value.split('=',1)[1],
                                     'sha256': digest(Path(value.split('=',1)[1]))} for value in templates]
        recipe['weighting_version'] = 1
        recipe['template_provenance'] = provenance
        path.write_text(json.dumps(recipe, indent=2)+'\n')


def prepare_jobs(root):
    from reec_dock_panel import fingerprint
    engine = digest(root / 'bin/aromadock')
    jobs, panels = {}, {}
    for name in ('top-panel', 'full-panel'):
        panel = json.loads((root / 'results' / f'{name}.json').read_text())
        prepared = json.loads((root / 'results' / f'{name}.prepared.json').read_text())
        if prepared['pending']:
            raise RuntimeError(f'{name} has unresolved models')
        expected = sum(len(case['conditions']) for case in panel['cases'])
        if len(prepared['cases']) != expected:
            raise RuntimeError(f'{name} has an incomplete configuration manifest')
        panels[name] = panel
        for case in prepared['cases']:
            for clock in panel['settings']['clock_seeds']:
                identifier = fingerprint(case, panel['settings'], engine) + '-' + str(clock)
                alias = {'panel': name, 'case': case['case'], 'condition': case['condition'],
                         'sdf_files': [str(Path(sdf['path']).relative_to(root)) for sdf in case['sdfs']]}
                if identifier in jobs:
                    jobs[identifier]['aliases'].append(alias)
                    continue
                folder = root / 'jobs' / identifier
                folder.mkdir(parents=True, exist_ok=True)
                lines = []
                for line in Path(case['config']).read_text().splitlines():
                    fields = line.split()
                    if fields and fields[0] in ('OUT','OUTPDB','POSE','ITERS','ITER','THREADS','THREAD','T','RETRY'):
                        continue
                    lines.append(line)
                settings = panel['settings']
                lines.extend([f'OUT {folder}/result.dock', f'OUTPDB {settings["poses"]} {folder}/pose%o.pdb',
                              f'POSE {settings["poses"]}', f'ITERS {settings["iterations"]}', 'THREADS 1', 'RETRY 0'])
                (folder / 'job.config').write_text('\n'.join(lines)+'\n')
                jobs[identifier] = {'id': identifier, 'receptor': case['case']['receptor'],
                                    'ligand': case['case']['ligand'], 'clock_seed': clock, 'aliases': [alias]}
    manifest = {'remote_root': str(root), 'jobs': list(jobs.values()), 'panels': panels,
                'template_weighting_version': 1, 'workers': 5, 'thermal_limit_c': 75,
                'reec_sha256': digest(root / 'bin/reec'), 'engine_sha256': engine}
    (root / 'jobs.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return len(jobs)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    os.chdir(root)
    phases = []
    cancel = threading.Event()
    for number in (signal.SIGINT, signal.SIGTERM):
        signal.signal(number, lambda signum, frame: cancel.set())

    def phase(name, command, timeout):
        print(f'Phase {name}', flush=True)
        phases.append({'phase': name, 'status': 'running', 'started_at_unix': time.time()})
        status = root / 'refresh-status.json'
        status.write_text(json.dumps(phases, indent=2)+'\n')
        with (root / f'{name}.log').open('w') as log:
            result = run_guarded(command, cwd=root, environment=os.environ.copy(), log=log,
                                 temperature_log=root / f'{name}.temperature.jsonl',
                                 timeout=timeout, pause_at=71, resume_below=66, cancel_event=cancel)
        phases[-1].update(result, status='complete' if result['exit_code'] == 0 else 'error')
        status.write_text(json.dumps(phases, indent=2)+'\n')
        if result['exit_code'] != 0:
            raise RuntimeError(f'{name} failed: {result}')

    phase('build', ['make','-j2','bin/reec','bin/phew','bin/ic','bin/cavity_search','test/moiety_test','test/reec_test'], 1800)
    objects = [str(root / 'obj' / (name+'.o')) for name in
               ('misc','point','atom','space','intera','molecule','aminoacid','protein','moiety','conj','progress')]
    phase('audit-build', ['g++','src/test/reec_steric_audit.cpp',*objects,'-o','bin/reec_steric_audit',
                          '-O3','-flto','-fopenmp','-std=c++14','-march=native'], 900)
    phase('native-checks', ['./test/reec_test'], 600)
    refresh_recipes(root)
    phase('paired-models', [sys.executable,'-u','hm/reec_docking_baseline.py','--workspace',str(root),
                            '--workers','1','--models','10','--max-attempts','5'], 43200)
    for name in ('top-panel', 'full-panel'):
        phase(name+'-prepare', [sys.executable,'-u','hm/reec_dock_panel.py','prepare','--workspace',str(root),
                                '--panel',str(root / 'results' / f'{name}.json')], 14400)
    count = prepare_jobs(root)
    print(f'Prepared {count} unique docking jobs', flush=True)
    if cancel.is_set():
        raise SystemExit(130)
    # The runner writes per-job checkpoints; replacing this process preserves nohup.
    os.execv(sys.executable, [sys.executable,'-u','hm/reec_remote_runner.py','--root',str(root),
                             '--workers','5','--pause-at','71','--resume-below','66'])


if __name__ == '__main__':
    main()
