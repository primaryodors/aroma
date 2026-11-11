"""Temperature-gated subprocess execution for CPU-intensive benchmark jobs."""
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import time


def cpu_temperature(reading):
    values = []
    for chip, features in reading.items():
        if not isinstance(features, dict):
            continue
        cpu_chip = chip.startswith(('k10temp', 'coretemp', 'zenpower'))
        for label, fields in features.items():
            if not isinstance(fields, dict):
                continue
            cpu_label = label in ('CPU', 'CPU Temperature', 'Tctl', 'Tdie') or label.startswith(('Package id ', 'Core '))
            if not (cpu_chip or cpu_label):
                continue
            for key, value in fields.items():
                if key.startswith('temp') and key.endswith('_input') and isinstance(value, (int, float)) and math.isfinite(value):
                    values.append(float(value))
    if not values:
        raise RuntimeError('No CPU temperature is available from sensors')
    return max(values)


def read_temperature():
    result = subprocess.run(['sensors', '-j'], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                            text=True, timeout=1)
    return cpu_temperature(json.loads(result.stdout))


def run_guarded(command, *, cwd, environment, log, temperature_log, timeout=1800,
                pause_at=71., resume_below=66., cancel_event=None):
    process = None
    stopped = False
    started = time.monotonic()
    active = 0.
    maximum = None
    paused_seconds = 0.
    previous = started
    last_report = 0.
    events = Path(temperature_log).open('a')
    try:
        while True:
            if cancel_event is not None and cancel_event.is_set():
                code = 'interrupted'
                break
            now = time.monotonic()
            elapsed = now - previous
            previous = now
            if process is not None and not stopped:
                active += elapsed
            else:
                paused_seconds += elapsed
            try:
                temperature = read_temperature()
                maximum = temperature if maximum is None else max(maximum, temperature)
                error = None
            except (OSError, ValueError, RuntimeError, subprocess.TimeoutExpired) as failure:
                temperature, error = None, str(failure)
            cool = temperature is not None and temperature <= resume_below
            hot = temperature is None or temperature >= pause_at
            if process is None and cool:
                process = subprocess.Popen(command, cwd=cwd, env=environment, stdout=log,
                                           stderr=subprocess.STDOUT, start_new_session=True)
            elif process is not None and process.poll() is None:
                if hot and not stopped:
                    os.killpg(process.pid, signal.SIGSTOP)
                    stopped = True
                elif cool and stopped:
                    os.killpg(process.pid, signal.SIGCONT)
                    stopped = False
            state = 'waiting' if process is None else ('cooling' if stopped else 'running')
            events.write(json.dumps({'elapsed_seconds': now-started, 'temperature_c': temperature,
                                     'state': state, 'sensor_error': error}) + '\n')
            events.flush()
            if state != 'running' and now-last_report >= 30:
                print(f'Thermal control: {state}; CPU {temperature} C', flush=True)
                last_report = now
            if process is not None and process.poll() is not None:
                code = process.returncode
                break
            if active >= timeout:
                code = 'timeout'
                break
            time.sleep(1)
        return {'exit_code': code, 'active_seconds': active, 'cooling_seconds': paused_seconds,
                'wall_seconds': time.monotonic()-started, 'max_temperature_c': maximum}
    finally:
        if process is not None and process.poll() is None:
            if stopped:
                os.killpg(process.pid, signal.SIGCONT)
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
        events.close()
