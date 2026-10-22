import os
import subprocess
import unittest


class ConfigInputTest(unittest.TestCase):
    def test_invalid_options(self):
        binary = os.environ.get("AROMA_BIN", "bin/aromadock")
        cases = [
            ["--prot", "A" * 400], ["--prot"], ["--lig"], ["--smiles"],
            ["--pose", "-1"], ["--pose", "0"], ["--pose", "1000000000"],
            ["--size", "1", "2"], ["--cen", "ABS", "1"],
            ["--path", "-1", "ABS", "1", "2", "3"],
            ["--outpdb", "1"], ["--bridge", "1"], ["--mcoord", "Cu"],
            ["--iters", "bad"],
        ]
        for arguments in cases:
            with self.subTest(arguments=arguments):
                result = subprocess.run([binary, *arguments], capture_output=True, text=True, timeout=10)
                self.assertEqual(result.returncode, 2, result.stderr)
                self.assertIn("Configuration error:", result.stderr)
                self.assertNotIn("AddressSanitizer", result.stderr)


if __name__ == "__main__":
    unittest.main()
