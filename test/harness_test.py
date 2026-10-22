import pathlib
import shutil
import subprocess
import tempfile
import unittest


class HarnessTest(unittest.TestCase):
    def test_build_failure(self):
        with tempfile.TemporaryDirectory() as root:
            root = pathlib.Path(root)
            (root / "test").mkdir()
            shutil.copyfile("test/unit_tests.sh", root / "test/unit_tests.sh")
            (root / "makefile").write_text("tests aromadock:\n\t@false\n")
            result = subprocess.run(["bash", str(root / "test/unit_tests.sh")],
                                    capture_output=True, timeout=10)
            self.assertNotEqual(result.returncode, 0)

    def test_missing_executables(self):
        with tempfile.TemporaryDirectory() as root:
            root = pathlib.Path(root)
            (root / "test").mkdir()
            shutil.copyfile("test/unit_tests.sh", root / "test/unit_tests.sh")
            (root / "makefile").write_text("tests aromadock:\n\t@true\n")
            result = subprocess.run(["bash", str(root / "test/unit_tests.sh")],
                                    capture_output=True, timeout=10)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn(b"execution failed", result.stdout)


if __name__ == "__main__":
    unittest.main()
