import pathlib
import json
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

    def convert(self, **values):
        body = urllib.parse.urlencode(values).encode()
        request = urllib.request.Request(self.base + "/smiles_view.php", data=body)
        with urllib.request.urlopen(request, timeout=30) as response:
            self.assertIn("no-store", response.headers.get("Cache-Control", ""))
            return json.load(response)

    def test_each_response_contains_its_own_structure(self):
        stale = pathlib.Path("www/test.sdf")
        before = stale.read_bytes() if stale.exists() else None
        ethanol = self.convert(smiles="CCO")
        ethane = self.convert(smiles="CC")
        self.assertEqual(int(ethanol["sdf"].splitlines()[3][:3]), 9)
        self.assertEqual(int(ethane["sdf"].splitlines()[3][:3]), 8)
        self.assertNotEqual(ethanol["sdf"], ethane["sdf"])
        self.assertIn("RDKit", ethanol["log"])
        self.assertEqual(stale.read_bytes() if stale.exists() else None, before)

    def test_transform_shell_text_is_rejected(self):
        body = urllib.parse.urlencode({"smiles": "CCO", "brot": "C1 C2 90; printf marker"}).encode()
        request = urllib.request.Request(self.base + "/smiles_view.php", data=body)
        with self.assertRaises(urllib.error.HTTPError) as result:
            urllib.request.urlopen(request, timeout=5)
        self.assertEqual(result.exception.code, 400)
        self.assertIn("error", json.load(result.exception))

    def test_malformed_atom_label_is_rejected(self):
        body = urllib.parse.urlencode({"smiles": "C{ABCDEFGHIJKLMNO}"}).encode()
        request = urllib.request.Request(self.base + "/smiles_view.php", data=body)
        with self.assertRaises(urllib.error.HTTPError) as result:
            urllib.request.urlopen(request, timeout=5)
        self.assertEqual(result.exception.code, 400)
        self.assertIn("error", json.load(result.exception))

    def test_conversion_failure_has_no_structure(self):
        body = urllib.parse.urlencode({"smiles": "not a molecule"}).encode()
        request = urllib.request.Request(self.base + "/smiles_view.php", data=body)
        with self.assertRaises(urllib.error.HTTPError) as result:
            urllib.request.urlopen(request, timeout=30)
        payload = json.load(result.exception)
        self.assertIn("error", payload)
        self.assertNotIn("sdf", payload)

    def test_rotations_and_ring_flips_modify_generated_coordinates(self):
        plain = self.convert(smiles="C1CCCCC1CC")
        rotated = self.convert(smiles="C1CCCCC1CC", brot="C6 C7 60")
        flipped = self.convert(smiles="C1CCCCC1CC", rflp="C1")
        combined = self.convert(smiles="C1CCCCC1CC", brot="C6 C7 60", rflp="C1")
        def coordinates(result):
            rows = result["sdf"].splitlines()
            count = int(rows[3][:3])
            return [tuple(float(row[i:i+10]) for i in (0, 10, 20))
                    for row in rows[4:4+count]]
        original = coordinates(plain)
        for result in (rotated, flipped, combined):
            points = coordinates(result)
            self.assertEqual(len(points), len(original))
            self.assertTrue(any(abs(a-b) > 0.05 for p, q in zip(original, points)
                                for a, b in zip(p, q)))
        self.assertNotEqual(coordinates(combined), coordinates(rotated))
        self.assertNotEqual(coordinates(combined), coordinates(flipped))


if __name__ == "__main__":
    unittest.main()
