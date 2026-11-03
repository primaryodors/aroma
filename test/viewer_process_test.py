import json
import shutil
import subprocess
import unittest


@unittest.skipUnless(shutil.which("php"), "PHP is unavailable")
class ViewerProcessTest(unittest.TestCase):
    def test_disabled_execution_is_reported(self):
        code = '''require "www/chem_process.php";
            try { chem_process_result(["obabel"]); }
            catch (Throwable $error) { echo json_encode([$error->getCode(), $error->getMessage()]); }'''
        result = subprocess.run(["php", "-d", "disable_functions=proc_open", "-r", code],
                                capture_output=True, text=True, check=True, timeout=5)
        status, message = json.loads(result.stdout)
        self.assertEqual(status, 503)
        self.assertIn("disabled", message)

    def test_shell_text_is_passed_literally(self):
        code = '''require "www/chem_process.php";
            $input = '{C}$(printf PROBE); " quote';
            $result = chem_process_result([PHP_BINARY, "-r", 'echo $argv[1];', $input]);
            echo json_encode([$input, $result]);'''
        result = subprocess.run(["php", "-r", code], capture_output=True, text=True,
                                check=True, timeout=5)
        original, output = json.loads(result.stdout)
        self.assertEqual(output["status"], 0)
        self.assertEqual(output["output"], original)

    def test_timeout_is_reported(self):
        code = '''require "www/chem_process.php";
            try { chem_process_result([PHP_BINARY, "-r", "sleep(5);"], null, 0.05); }
            catch (Throwable $error) { echo $error->getCode(); }'''
        result = subprocess.run(["php", "-r", code], capture_output=True, text=True,
                                check=True, timeout=5)
        self.assertEqual(result.stdout, "504")


if __name__ == "__main__":
    unittest.main()
