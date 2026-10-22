import pathlib
import shutil
import socket
import subprocess
import time
import unittest
import urllib.error
import urllib.parse
import urllib.request


@unittest.skipUnless(shutil.which("php"), "PHP is unavailable")
class WebSecurityTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with socket.socket() as sock:
            sock.bind(("127.0.0.1", 0))
            port = sock.getsockname()[1]
        cls.base = f"http://127.0.0.1:{port}"
        cls.server = subprocess.Popen(["php", "-S", f"127.0.0.1:{port}", "-t", "www"],
                                      stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        for _ in range(100):
            try:
                with urllib.request.urlopen(cls.base + "/instcheck.php", timeout=1):
                    return
            except (OSError, urllib.error.URLError):
                time.sleep(0.02)
        cls.server.terminate()
        cls.server.wait(timeout=5)
        raise RuntimeError("Test server did not start")

    @classmethod
    def tearDownClass(cls):
        cls.server.terminate()
        cls.server.wait(timeout=5)

    def test_maintenance_is_terminal_only(self):
        database = pathlib.Path("data/receptor.json")
        before = database.read_bytes()
        with self.assertRaises(urllib.error.HTTPError) as result:
            urllib.request.urlopen(self.base + "/getpdbs.php", timeout=5)
        self.assertEqual(result.exception.code, 403)
        self.assertEqual(database.read_bytes(), before)

    def test_input_is_escaped(self):
        value = '\"><script>alert(1)</script>'
        query = urllib.parse.urlencode({"s": value})
        with urllib.request.urlopen(self.base + "/smiles_view.php?" + query, timeout=5) as response:
            page = response.read().decode()
        self.assertNotIn(value, page)
        self.assertIn("&lt;script&gt;", page)


if __name__ == "__main__":
    unittest.main()
