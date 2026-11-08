#!/usr/bin/env python3
"""Build the frozen native benchmark on its execution host under thermal control."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import threading

from reec_thermal import run_guarded


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    cancelled = threading.Event()
    for number in (signal.SIGINT, signal.SIGTERM):
        signal.signal(number, lambda signum, frame: cancelled.set())
    for name in ('bin', 'obj', 'tmp', 'out'):
        (root / name).mkdir(exist_ok=True)
    reports = []
    commands = [['make', '-j2', 'bin/aromadock'],
                ['gcc', '-shared', '-fPIC', '-O2', 'test/reec_baseline_clock.c', '-ldl', '-o', 'bin/baseline_clock.so']]
    for index, command in enumerate(commands):
        with (root / f'build-{index}.log').open('w') as log:
            result = run_guarded(command, cwd=root, environment=dict(os.environ), log=log,
                                 temperature_log=root / f'build-{index}.temperature.jsonl',
                                 cancel_event=cancelled)
        reports.append({'command': command, **result})
        if result['exit_code'] != 0:
            (root / 'build.json').write_text(json.dumps({'steps': reports}, indent=2)+'\n')
            return 1
    version = subprocess.run(['g++', '--version'], capture_output=True, text=True, check=True).stdout.splitlines()[0]
    report = {'host': os.uname().nodename, 'compiler': version, 'steps': reports,
              'engine_sha256': hashlib.sha256((root / 'bin/aromadock').read_bytes()).hexdigest()}
    (root / 'build.json').write_text(json.dumps(report, indent=2)+'\n')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
