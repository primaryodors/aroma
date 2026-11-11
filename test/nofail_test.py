import os
import pathlib
import subprocess
import tempfile
import unittest
from rdkit import Chem
from rdkit.Chem import AllChem


class NoFailTest(unittest.TestCase):
    def check_rejected_export(self, ligand_stem='ligand', placeholders=False):
        rows = pathlib.Path("pdbs/OR51/OR51E2.active.pdb").read_text().splitlines()
        residues = []
        selected = []
        for row in rows:
            if not row.startswith("ATOM  "):
                continue
            residue = row[21:27]
            if residue not in residues:
                if len(residues) == 5:
                    break
                residues.append(residue)
            selected.append(row)
        self.assertEqual(len(residues), 5)
        with tempfile.TemporaryDirectory(prefix="aroma-nofail-") as directory:
            root = pathlib.Path(directory)
            receptor = root / "receptor.pdb"
            receptor.write_text("\n".join(selected) + "\nTER\nEND\n")
            ligand = root / f"{ligand_stem}.sdf"
            molecule = Chem.AddHs(Chem.MolFromSmiles("CCC(=O)O"))
            self.assertEqual(AllChem.EmbedMolecule(molecule, randomSeed=1), 0)
            with Chem.SDWriter(str(ligand)) as writer:
                writer.write(molecule)
            dock = root / (f"receptor~{ligand_stem}.dock" if placeholders else "result.dock")
            pdb = root / (f"receptor~{ligand_stem}.model1.pdb" if placeholders else "result.pdb")
            dock_pattern = root / "%p~%l.dock" if placeholders else dock
            pdb_pattern = root / "%p~%l.model%o.pdb" if placeholders else pdb
            config = root / "run.config"
            config.write_text(
                f"PROT {receptor}\nLIG {ligand}\n"
                "CEN 100 100 100\nSIZE 2 2 2\nPOSE 2\nITERS 1\n"
                f"SEARCH CF\nELIM -1e20\nOUT {dock_pattern}\nOUTPDB 1 {pdb_pattern}\nNOFAIL\nAPPENDPROT\n"
            )
            result = subprocess.run([os.environ.get("AROMA_BIN", "bin/aromadock"), str(config)],
                                    text=True, capture_output=True, timeout=90)
            self.assertEqual(result.returncode, 0, result.stderr)
            report = dock.read_text()
            coordinates = pdb.read_text()
            self.assertIn("Pose: 1", report)
            self.assertIn("\nTotal: 0\n", report)
            self.assertIn("Disqualified because:", report)
            self.assertIn("0 pose(s) found.", report)
            self.assertIn("PDBDAT:", report)
            self.assertNotIn("Missing PDB data", report)
            self.assertIn("ATOM  ", coordinates)
            self.assertIn("HETATM", coordinates)
            original = report.split("\nOriginal PDB:\n", 1)[1]
            def positions(lines):
                return {(line[22:26].strip(), line[12:16].strip()):
                        tuple(float(line[start:start + 8]) for start in (30, 38, 46))
                        for line in lines if line.startswith("ATOM  ")
                        and not line[12:16].strip().lstrip("123").startswith("H")}
            self.assertEqual(positions(selected), positions(original.splitlines()))

    def test_rejected_zero_energy_pose_is_exported(self):
        self.check_rejected_export()

    def test_long_ligand_filename_exports_diagnostic_model(self):
        self.check_rejected_export('ligand_'+'x'*80)

    def test_long_ligand_filename_expands_output_placeholders(self):
        self.check_rejected_export('ligand_'+'x'*80, placeholders=True)


if __name__ == "__main__":
    unittest.main()
