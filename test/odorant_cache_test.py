import os
import pathlib
import tempfile
import subprocess
import unittest
from unittest import mock

try:
    from rdkit import Chem
except ImportError:
    Chem = None


@unittest.skipIf(Chem is None, "RDKit is unavailable")
class OdorantCacheTest(unittest.TestCase):
    def setUp(self):
        from data import odorutils
        self.utils = odorutils
        self.original_cwd = os.getcwd()
        self.directory = tempfile.TemporaryDirectory(prefix="aroma-sdf-cache-")
        self.path = str(pathlib.Path(self.directory.name) / "ligand.sdf")

    def tearDown(self):
        os.chdir(self.original_cwd)
        self.directory.cleanup()

    def write_structure(self, smiles, path):
        with Chem.SDWriter(path) as writer:
            writer.write(Chem.MolFromSmiles(smiles.split("|", 1)[0]))

    def ensure(self, smiles, generator):
        odorant = {"full_name": "fixture", "smiles": smiles, "sdfname": self.path}
        with mock.patch.object(self.utils, "find_odorant", return_value=odorant), \
                mock.patch.object(self.utils, "check_isomers", return_value=False), \
                mock.patch.object(self.utils, "check_forms", return_value=False), \
                mock.patch.object(self.utils, "smiles_to_sdf", side_effect=generator) as generate:
            result = self.utils.ensure_sdf_exists("fixture")
            self.assertEqual(result, self.path)
            return generate.call_count

    def test_replaces_a_positional_isomer_with_the_same_formula(self):
        self.write_structure("C1CCCCCCOC(=O)/C=C/CCCCC1", self.path)
        corrected = "C1CCCC/C=C/CCCOC(=O)CCCC1"
        self.assertEqual(self.ensure(corrected, self.write_structure), 1)
        actual = Chem.SDMolSupplier(self.path)[0]
        self.assertEqual(Chem.MolToSmiles(actual), Chem.MolToSmiles(Chem.MolFromSmiles(corrected)))

    def test_preserves_an_unchanged_cached_conformation(self):
        self.write_structure("CCO", self.path)
        original = pathlib.Path(self.path).read_bytes()
        self.assertEqual(self.ensure("CCO", self.write_structure), 0)
        self.assertEqual(pathlib.Path(self.path).read_bytes(), original)

    def test_replaces_the_wrong_explicit_alkene_stereochemistry(self):
        self.write_structure("C/C=C\\C", self.path)
        self.assertEqual(self.ensure("C/C=C/C", self.write_structure), 1)

    def test_retains_unspecified_stereocentres(self):
        self.write_structure("C[C@H](O)C(=O)O", self.path)
        self.assertEqual(self.ensure("CC(O)C(=O)O", self.write_structure), 0)

    def test_reapplies_conformer_directives(self):
        self.write_structure("C1CCCCC1", self.path)
        self.assertEqual(self.ensure("C1CCCCC1|rflp:1", self.write_structure), 1)

    def test_rejects_a_failed_replacement(self):
        self.write_structure("CCC", self.path)
        with self.assertRaisesRegex(RuntimeError, "does not match"):
            self.ensure("CCO", lambda *args: None)

    def test_creates_a_missing_file(self):
        self.assertEqual(self.ensure("CCO", self.write_structure), 1)

    def test_native_ringflip_preserves_readable_carvone_and_stereochemistry(self):
        binary = pathlib.Path(__file__).resolve().parents[1] / 'bin/ringflip'
        if not binary.is_file():
            self.skipTest('Build bin/ringflip to exercise native SDF export')
        smiles = 'CC(=C)[C@@H]1CC=C(C(=O)C1)C|rflp:C4'
        self.utils._ensure_sdf_current(smiles, self.path)
        self.assertTrue(self.utils._sdf_matches_smiles(self.path, smiles))
        lines = pathlib.Path(self.path).read_text().splitlines()
        count = int(lines[3][:3])
        bonds = int(lines[3][3:6])
        types = [int(line[6:9]) for line in lines[4+count:4+count+bonds]]
        self.assertIn(2, types)
        self.assertTrue(all(value in (1, 2, 3, 4) for value in types))

    def test_menthol_preparation_applies_both_flips_and_preserves_stereochemistry(self):
        self.utils.load_odors()
        odor = self.utils.find_odorant('(1R,2S,5R)-(-)-menthol')
        self.assertTrue(odor)
        native_run = subprocess.run
        with mock.patch.object(self.utils.subprocess, 'run', wraps=native_run) as run:
            self.utils._ensure_sdf_current(odor['smiles'], self.path)
        flips = [call.args[0] for call in run.call_args_list if call.args[0][0] == 'bin/ringflip']
        self.assertEqual(len(flips), 1)
        self.assertEqual(flips[0][2:], ['C2', 'C3'])
        actual = Chem.SDMolSupplier(self.path)[0]
        expected = Chem.MolFromSmiles(odor['smiles'].split('|')[0])
        self.assertIsNotNone(actual)
        self.assertEqual(Chem.MolToSmiles(actual), Chem.MolToSmiles(expected))
        self.assertEqual(len(Chem.FindMolChiralCenters(actual)), 3)
        lines = pathlib.Path(self.path).read_text().splitlines()
        atoms, bonds = int(lines[3][:3]), int(lines[3][3:6])
        self.assertTrue(all(int(line[6:9]) == 1 for line in lines[4+atoms:4+atoms+bonds]))

    def test_failed_ringflip_preserves_the_cached_file(self):
        self.write_structure('C1CCCCC1', self.path)
        original = pathlib.Path(self.path).read_bytes()
        def command(arguments, **options):
            if arguments[0] == 'make':
                return subprocess.CompletedProcess(arguments, 0)
            pathlib.Path(arguments[1]).write_text('invalid SDF')
            raise subprocess.CalledProcessError(1, arguments)
        with mock.patch.object(self.utils.subprocess, 'run', side_effect=command):
            with self.assertRaises(subprocess.CalledProcessError):
                self.utils._ensure_sdf_current('C1CCCCC1|rflp:C1', self.path)
        self.assertEqual(pathlib.Path(self.path).read_bytes(), original)


if __name__ == "__main__":
    unittest.main()
