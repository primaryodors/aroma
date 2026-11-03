import os
import pathlib
import tempfile
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


if __name__ == "__main__":
    unittest.main()
