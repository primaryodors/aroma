import pathlib
import tempfile
import unittest
from data.disulfides import pdb_disulfides


class DisulfidesTest(unittest.TestCase):
    def test_preserves_all_declared_links_without_a_distance_filter(self):
        atoms = "".join(
            f"ATOM  {serial:5d}  SG  CYS A{residue:4d}    {x:8.3f}{0:8.3f}{0:8.3f}\n"
            for serial, residue, x in ((1, 20, 0), (2, 90, 3.2), (3, 40, 10), (4, 70, 12.05)))
        with tempfile.TemporaryDirectory() as directory:
            path = pathlib.Path(directory) / "model.pdb"
            path.write_text(atoms + "CONECT    1    2\nCONECT    3    4\nCONECT    2    1\n")
            self.assertEqual(pdb_disulfides(path), [(20, 90), (40, 70)])
            self.assertEqual(pdb_disulfides(path, "B"), [])

    def test_reads_ssbond_and_ignores_other_atom_links(self):
        atoms = "".join(
            f"ATOM  {serial:5d}  SG  CYS A{residue:4d}    {x:8.3f}{0:8.3f}{0:8.3f}\n"
            for serial, residue, x in ((1, 10, 0), (2, 80, 2.05)))
        with tempfile.TemporaryDirectory() as directory:
            path = pathlib.Path(directory) / "model.pdb"
            path.write_text(atoms + "SSBOND   1 CYS A   10    CYS A   80\nCONECT    1  999\n")
            self.assertEqual(pdb_disulfides(path), [(10, 80)])


if __name__ == "__main__":
    unittest.main()
