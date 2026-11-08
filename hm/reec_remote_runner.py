#!/usr/bin/env python3
"""Run a frozen docking manifest with temperature control and checkpoints."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import threading
import time

from reec_thermal import read_temperature, run_guarded

CANCEL = threading.Event()


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic_json(path, value):
    temporary = path.with_suffix(path.suffix + '.new')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def execute(root, job, engine, pause_at, resume_below, state):
    if CANCEL.is_set():
        return {'exit_code': 'interrupted'}
    if not re.fullmatch(r'[0-9a-f]{64}-[0-9]+', job['id']):
        raise ValueError('Invalid benchmark job identifier')
    folder = root / 'jobs' / job['id']
    folder.mkdir(parents=True, exist_ok=True)
    config = folder / 'job.config'
    result_path = folder / 'result.json'
    configuration_hash = digest(config)
    if result_path.exists():
        result = json.loads(result_path.read_text())
        if result['engine_sha256'] != engine or result['config_sha256'] != configuration_hash:
            raise RuntimeError('Checkpoint does not match the frozen engine and configuration')
        state(job['id'], 'complete')
        return result
    for name in ('data', 'sdf'):
        path = folder / name
        if not path.exists():
            path.symlink_to(root / name)
    (folder / 'tmp').mkdir(exist_ok=True)
    for path in [folder / 'result.dock', *folder.glob('pose*.pdb')]:
        path.unlink(missing_ok=True)
    environment = dict(os.environ, LD_PRELOAD=str(root / 'bin/baseline_clock.so'),
                       AROMA_BASELINE_TIME=str(job['clock_seed']), OMP_NUM_THREADS='1')
    command = [str(root / 'bin/aromadock'), str(config), '--ns', '-t', '1']
    state(job['id'], 'running')
    print(f'Start {job["receptor"]} {job["ligand"]}, seed {job["clock_seed"]}', flush=True)
    with (folder / 'engine.log').open('w') as log:
        thermal = run_guarded(command, cwd=folder, environment=environment, log=log,
                              temperature_log=folder / 'temperature.jsonl', timeout=1800,
                              pause_at=pause_at, resume_below=resume_below, cancel_event=CANCEL)
    output = folder / 'result.dock'
    result = {'job_id': job['id'], 'engine_sha256': engine, 'config_sha256': configuration_hash,
              'clock_seed': job['clock_seed'], 'command': command, **thermal,
              'dock_sha256': digest(output) if output.exists() else None,
              'poses': [str(path) for path in sorted(folder.glob('pose*.pdb'))]}
    if thermal['exit_code'] == 'interrupted':
        atomic_json(folder / 'interrupted.json', result)
        state(job['id'], 'interrupted')
    else:
        atomic_json(result_path, result)
        state(job['id'], 'complete')
    print(f'Finish {job["receptor"]}: exit {thermal["exit_code"]}, {thermal["wall_seconds"]:.1f}s, max {thermal["max_temperature_c"]} C', flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=5)
    parser.add_argument('--pause-at', type=float, default=71.)
    parser.add_argument('--resume-below', type=float, default=66.)
    args = parser.parse_args()
    if not 1 <= args.workers <= 6 or not args.resume_below < args.pause_at < 75:
        parser.error('Use one to six workers and a pause threshold below 75 C')
    root = args.root.resolve()
    manifest = json.loads((root / 'jobs.json').read_text())
    engine = digest(root / 'bin/aromadock')
    signals = {number: signal.getsignal(number) for number in (signal.SIGINT, signal.SIGTERM)}
    for number in signals:
        signal.signal(number, lambda signum, frame: CANCEL.set())
    lock = threading.Lock()
    status = {'host': os.uname().nodename, 'engine_sha256': engine, 'workers': args.workers,
              'pause_at_c': args.pause_at, 'resume_below_c': args.resume_below,
              'initial_temperature_c': read_temperature(), 'states': {job['id']: 'queued' for job in manifest['jobs']}}

    def state(identifier, value):
        with lock:
            status['states'][identifier] = value
            status['updated_at_unix'] = time.time()
            atomic_json(root / 'status.json', status)

    atomic_json(root / 'status.json', status)
    results = {}
    try:
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            pending = {pool.submit(execute, root, job, engine, args.pause_at, args.resume_below, state): job for job in manifest['jobs']}
            for future in as_completed(pending):
                job = pending[future]
                try:
                    results[job['id']] = future.result()
                except Exception as error:
                    results[job['id']] = {'error': str(error)}
                    state(job['id'], 'error')
                    print(f'Job error: {error}', flush=True)
                atomic_json(root / 'results.json', results)
    finally:
        for number, previous in signals.items():
            signal.signal(number, previous)
    if CANCEL.is_set():
        return 130
    return 0 if all(result.get('exit_code') == 0 for result in results.values()) else 2


if __name__ == '__main__':
    raise SystemExit(main())
