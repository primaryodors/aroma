import pathlib
import os
import subprocess
import tempfile
import unittest
from unittest.mock import Mock, patch
from data.preparation import PreparationCache, prepared_output_current


class PreparationTest(unittest.TestCase):
    def test_output_refreshes_for_every_dependency(self):
        with tempfile.TemporaryDirectory() as root:
            output = pathlib.Path(root) / "pocket.cvty"
            inputs = [pathlib.Path(root) / name for name in ("model.pdb", "scanner", "options.json", "atoms.dat")]
            for path in inputs:
                path.write_text("input")
                os.utime(path, ns=(1000000000, 1000000000))
            self.assertFalse(prepared_output_current(output, inputs))
            output.write_text("cavity")
            os.utime(output, ns=(2000000000, 2000000000))
            self.assertTrue(prepared_output_current(output, inputs))
            for path in inputs:
                os.utime(path, ns=(3000000000, 3000000000))
                self.assertFalse(prepared_output_current(output, inputs))
                os.utime(path, ns=(1000000000, 1000000000))
            inputs[-1].unlink()
            self.assertFalse(prepared_output_current(output, inputs))

    def test_reuses_scan_and_refreshes_changed_input(self):
        with tempfile.TemporaryDirectory() as root:
            binary = pathlib.Path(root) / "scanner"
            protein = pathlib.Path(root) / "model.pdb"
            binary.write_text("scanner")
            protein.write_text("first")
            command = [str(binary), str(protein)]
            result = subprocess.CompletedProcess(command, 0, b"contacts")
            before = Mock()
            with patch("data.preparation.subprocess.run", return_value=result) as run:
                cache = PreparationCache()
                self.assertIs(cache.run(command, [protein], before), result)
                self.assertIs(cache.run(command, [protein], before), result)
                self.assertEqual(run.call_count, 1)
                protein.write_text("second model")
                cache.run(command, [protein], before)
                self.assertEqual(run.call_count, 2)
                self.assertEqual(before.call_count, 2)

    def test_failed_scans_are_not_cached(self):
        with tempfile.TemporaryDirectory() as root:
            binary = pathlib.Path(root) / "scanner"
            binary.write_text("scanner")
            command = [str(binary)]
            cache = PreparationCache()
            with patch("data.preparation.subprocess.run",
                       side_effect=subprocess.CalledProcessError(1, command)) as run:
                for _ in range(2):
                    with self.assertRaises(subprocess.CalledProcessError):
                        cache.run(command, [])
                self.assertEqual(run.call_count, 2)


if __name__ == "__main__":
    unittest.main()
