import json
import os
import pathlib
import shutil
import subprocess
import tempfile
import unittest


@unittest.skipUnless(shutil.which("php"), "PHP is unavailable")
class OdorantProcessTest(unittest.TestCase):
    def test_names_are_arguments_and_cwd_is_preserved(self):
        with tempfile.TemporaryDirectory(prefix="aroma-odorant-") as directory:
            root = pathlib.Path(directory)
            capture = root / "arguments.json"
            interpreter = root / "python capture"
            interpreter.write_text(
                "#!/usr/bin/python3\nimport json, os, sys\n"
                "with open(os.environ['AROMA_CAPTURE'], 'w') as output:\n"
                "    json.dump(sys.argv[1:], output)\n"
            )
            interpreter.chmod(0o700)
            environment = dict(os.environ, AROMA_PYTHON=str(interpreter), AROMA_CAPTURE=str(capture))
            name = "quoted ' odorant; $(printf marker) `printf marker` \""
            code = ('require "data/odorutils.php"; $before = getcwd(); '
                    'ensure_sdf_exists($argv[1]); echo json_encode([$before, getcwd()]);')
            result = subprocess.run(["php", "-r", code, "--", name], env=environment,
                                    text=True, capture_output=True, check=True)
            arguments = json.loads(capture.read_text())
            self.assertEqual(arguments[0], "-c")
            self.assertEqual(arguments[2:], [name])
            self.assertNotIn(name, arguments[1])
            self.assertIn("sys.argv[1]", arguments[1])
            before, after = json.loads(result.stdout)
            self.assertEqual(before, after)

    def test_process_failure_is_reported(self):
        environment = dict(os.environ, AROMA_PYTHON="/bin/false")
        code = ('require "data/odorutils.php"; try { ensure_sdf_exists("CCO"); } '
                'catch (RuntimeException $error) { echo $error->getMessage(); exit(7); }')
        result = subprocess.run(["php", "-r", code], env=environment, text=True, capture_output=True)
        self.assertEqual(result.returncode, 7)
        self.assertIn("Odorant generation failed", result.stdout)


if __name__ == "__main__":
    unittest.main()
