import os
import pathlib
import subprocess
import tempfile
import unittest


class NoFailTest(unittest.TestCase):
    def test_rejected_zero_energy_pose_is_exported(self):
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
            dock = root / "result.dock"
            pdb = root / "result.pdb"
            config = root / "run.config"
            config.write_text(
                f"PROT {receptor}\nLIG sdf/propionic_acid.sdf\n"
                "CEN 100 100 100\nSIZE 2 2 2\nPOSE 2\nITERS 1\n"
                f"SEARCH CF\nELIM -1e20\nOUT {dock}\nOUTPDB 1 {pdb}\nNOFAIL\n"
            )
            result = subprocess.run([os.environ.get("AROMA_BIN", "bin/aromadock"), str(config)],
                                    text=True, capture_output=True, timeout=90)
            self.assertEqual(result.returncode, 0, result.stderr)
            report = dock.read_text()
            coordinates = pdb.read_text()
            self.assertIn("Pose: 1", report)
            self.assertIn("\nTotal: 0\n", report)
            self.assertIn("Disqualified because:", report)
            self.assertIn("PDBDAT:", report)
            self.assertNotIn("Missing PDB data", report)
            self.assertIn("ATOM  ", coordinates)
            self.assertIn("HETATM", coordinates)


if __name__ == "__main__":
    unittest.main()
