"""CPU sensor selection and suspend/resume checks without real CPU load."""
import json
from pathlib import Path
import signal
import sys
import tempfile
import unittest
from unittest.mock import Mock, patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'hm'))
import reec_thermal as thermal


class ThermalControl(unittest.TestCase):
    def test_cpu_channels_exclude_gpu_storage_and_sensor_limits(self):
        reading = {'k10temp-pci': {'Tctl': {'temp1_input': 74.5, 'temp1_max': 115}},
                   'asus_wmi': {'CPU Temperature': {'temp1_input': 64},
                                'Tsensor 1 Temperature': {'temp5_input': 216}},
                   'nvme': {'Composite': {'temp1_input': 91}},
                   'amdgpu': {'edge': {'temp1_input': 88}}}
        self.assertEqual(thermal.cpu_temperature(reading), 74.5)
        with self.assertRaises(RuntimeError):
            thermal.cpu_temperature({'amdgpu': reading['amdgpu']})

    def test_valid_cpu_reading_survives_unrelated_sensor_warning(self):
        result = Mock(returncode=1, stdout=json.dumps({'thinkpad': {'CPU': {'temp1_input': 63}}}))
        with patch.object(thermal.subprocess, 'run', return_value=result):
            self.assertEqual(thermal.read_temperature(), 63)

    def test_hot_child_is_suspended_until_resume_threshold(self):
        temperatures = [70., 65., 72., 69., 65., 64.]
        progress = {'step': 0, 'clock': 0}

        def reading():
            value = temperatures[progress['step']]
            progress['step'] += 1
            return value

        def clock():
            progress['clock'] += 1
            return float(progress['clock'])

        child = Mock(pid=4242, returncode=0)
        child.poll.side_effect = lambda: 0 if progress['step'] == 6 else None
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            with (folder / 'engine.log').open('w') as log, \
                    patch.object(thermal, 'read_temperature', side_effect=reading), \
                    patch.object(thermal.time, 'monotonic', side_effect=clock), \
                    patch.object(thermal.time, 'sleep'), \
                    patch.object(thermal.subprocess, 'Popen', return_value=child) as start, \
                    patch.object(thermal.os, 'killpg') as signals:
                result = thermal.run_guarded(['dock'], cwd=folder, environment={}, log=log,
                                             temperature_log=folder/'temperature.jsonl')
            self.assertEqual(result['exit_code'], 0)
            self.assertEqual(result['max_temperature_c'], 72)
            self.assertEqual([call.args for call in signals.call_args_list],
                             [(4242, signal.SIGSTOP), (4242, signal.SIGCONT)])
            start.assert_called_once()


if __name__ == '__main__':
    unittest.main()
