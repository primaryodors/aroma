import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


class DocklistFallbackTest(unittest.TestCase):
    def test_fallback_is_visible_with_cold_and_warm_cache(self):
        reports = [
            "Pose: 1\nTotal: -58.5\nDisqualified because: rejected\n0 pose(s) found.\n",
            "Pose: 1\nTotal: 0\nNOFAIL: initial placement\n0 pose(s) found.\n",
            "Header\n" + "Details\n" * 80 + "Pose: 1\nTotal: -58.5\n0 pose(s) found.\n",
        ]
        for report in reports:
            with self.subTest(report=report[:30]), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                web = root / "www"
                data = root / "data"
                output = root / "out" / "OR51" / "OR51E2"
                for folder in (web, data, output):
                    folder.mkdir(parents=True)
                shutil.copyfile(ROOT / "www/docklist.php", web / "docklist.php")
                (data / "protutils.php").write_text(
                    "<?php $prots=['OR51E2'=>[]]; $prod='test'; "
                    "function family_from_protid($p){return 'OR51';} "
                    "function json_encode_pretty($v){return json_encode($v);}"
                )
                (data / "odorutils.php").write_text(
                    "<?php function find_odorant($v){return ['oid'=>'test','full_name'=>'Test ligand'];} "
                    "function best_empirical_pair($p,$o,$a){return false;}"
                )
                (data / "statistics.php").write_text("<?php function equilibrium($a,$b){return 0;}")
                (web / "dlmenu.php").write_text(
                    "<?php function output_dlmenu_style(){} function output_dlmenu_script(){}"
                )
                (web / "header.php").write_text("")
                name = "OR51E2~test.active.dock"
                dock = output / name
                dock.write_text(report)
                runner = (
                    "error_reporting(E_ERROR); ob_start(); $_SERVER['SERVER_NAME']='test'; "
                    "$_REQUEST=['r'=>'OR51E2']; if ($argv[1]==='fails') $_REQUEST['f']=1; "
                    "include 'docklist.php';"
                )
                for mode in ("all", "fails", "all"):
                    result = subprocess.run(["php", "-r", runner, mode], cwd=web,
                                            capture_output=True, text=True, timeout=15)
                    self.assertEqual(result.returncode, 0, result.stderr)
                    self.assertIn(
                        'mode=active" target="_dock">(fail)</a>', result.stdout
                    )
                    cache = json.loads((web / "docklist.cache.json").read_text())
                    self.assertEqual(cache[name]["benerg"], 1000)
                    self.assertEqual(cache[name]["nump"], "-")
                    # Force the next render through the cache branch.
                    import os
                    os.utime(dock, (1, 1))


if __name__ == "__main__":
    unittest.main()
