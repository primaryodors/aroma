import json
from pathlib import Path
import subprocess
import shutil
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
CASES = {
    "ethyl phenylacetate": "2935c1cf4d4eb001afcb4827a9c2b2c9",
    "ethylphenyl acetate": "3085d446d6702f5497196af67e16f143",
    "ethyl_phenylacetate": "2935c1cf4d4eb001afcb4827a9c2b2c9",
    "ethylphenyl_acetate": "3085d446d6702f5497196af67e16f143",
    "ETHYLPHENYL ACETATE": "3085d446d6702f5497196af67e16f143",
    "ethylphenylacetate": None,
}


class OdorantIdentityTest(unittest.TestCase):
    def test_php_names_and_filename_tokens(self):
        code = '''require 'data/odorutils.php';
        $result = [];
        foreach (json_decode($argv[1],true) as $query) {
            $odor = find_odorant($query);
            $result[$query] = $odor ? $odor['oid'] : null;
        }
        echo json_encode($result);'''
        result = subprocess.run(["php", "-r", code, json.dumps(list(CASES))],
                                cwd=ROOT, capture_output=True, text=True, check=True)
        self.assertEqual(json.loads(result.stdout), CASES)

    def test_python_names_and_filename_tokens(self):
        from data import odorutils
        odorutils.load_odors()
        for query, expected in CASES.items():
            with self.subTest(query=query):
                odor = odorutils.find_odorant(query)
                self.assertEqual(odor["oid"] if odor else None, expected)
        first = odorutils.find_odorant("ethyl phenylacetate")
        second = odorutils.find_odorant("ethylphenyl acetate")
        self.assertNotEqual(first["smiles"], second["smiles"])

    def test_docklist_keeps_both_esters_separate(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            web, data = root / "www", root / "data"
            output = root / "out/OR1/OR1A1"
            for folder in (web, data, output):
                folder.mkdir(parents=True)
            for name in ("docklist.php",):
                shutil.copyfile(ROOT / "www" / name, web / name)
            shutil.copyfile(ROOT / "data/odorutils.php", data / "odorutils.php")
            database = json.loads((ROOT / "data/odorant.json").read_text())
            selected = {oid: database[oid] for oid in set(CASES.values()) if oid}
            (data / "odorant.json").write_text(json.dumps(selected))
            (data / "reference.json").write_text("{}")
            (data / "protutils.php").write_text(
                "<?php $prots=['OR1A1'=>[]]; $prod='test'; "
                "function family_from_protid($p){return 'OR1';} "
                "function json_encode_pretty($v){return json_encode($v);}"
            )
            (data / "statistics.php").write_text("<?php function equilibrium($a,$b){return 0;}")
            (web / "dlmenu.php").write_text(
                "<?php function output_dlmenu_style(){} function output_dlmenu_script(){}"
            )
            (web / "header.php").write_text("")
            for token in ("ethyl_phenylacetate", "ethylphenyl_acetate"):
                (output / f"OR1A1~{token}.active.dock").write_text(
                    "Header\n" + "Details\n" * 80 + "Pose: 1\nTotal: -5\n1 pose(s) found.\n"
                )
            code = """error_reporting(E_ERROR); ob_start();
            $_SERVER['SERVER_NAME']='test'; $_REQUEST=['r'=>'OR1A1'];
            if ($argv[1]) $_REQUEST['o']=$argv[1]; include 'docklist.php';"""
            for oid in ("", CASES["ethylphenyl acetate"]):
                result = subprocess.run(["php", "-r", code, oid], cwd=web,
                                        capture_output=True, text=True, check=True)
                current_name = database[CASES['ethylphenyl acetate']]['full_name']
                self.assertIn(f'>{current_name}</a>', result.stdout)
                self.assertIn('odor=ethylphenyl_acetate&mode=active', result.stdout)
                if not oid:
                    self.assertIn('>ethyl phenylacetate</a>', result.stdout)
                    self.assertIn('odor=ethyl_phenylacetate&mode=active', result.stdout)
                else:
                    self.assertNotIn('>ethyl phenylacetate</a>', result.stdout)


if __name__ == "__main__":
    unittest.main()
