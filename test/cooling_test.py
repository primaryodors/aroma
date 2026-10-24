import types
import unittest
from unittest.mock import patch
import data.globals as globals_data


class CoolingTest(unittest.TestCase):
    def test_cool_processor_does_not_sleep(self):
        output = types.SimpleNamespace(stdout=b"Package id 0: +20.0 C\n")
        with patch.object(globals_data.subprocess, "run", return_value=output), \
             patch.object(globals_data.time, "sleep") as sleep:
            globals_data.wait_cool_cpu()
            sleep.assert_not_called()

    def test_hot_processor_still_waits(self):
        hot = types.SimpleNamespace(stdout=b"Package id 0: +120.0 C (high = +90.0 C)\n")
        cool = types.SimpleNamespace(stdout=b"Package id 0: +20.0 C (high = +90.0 C)\n")
        with patch.object(globals_data.subprocess, "run", side_effect=[hot, cool]), \
             patch.object(globals_data.time, "sleep") as sleep:
            globals_data.wait_cool_cpu()
            sleep.assert_called_once_with(60)

    def test_absent_monitor_does_not_sleep(self):
        with patch.object(globals_data.subprocess, "run", side_effect=FileNotFoundError), \
             patch.object(globals_data.time, "sleep") as sleep:
            globals_data.wait_cool_cpu()
            sleep.assert_not_called()


if __name__ == "__main__":
    unittest.main()
