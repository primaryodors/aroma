from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'hm'))
from or52e8_ligand_relax import extract_poses, merge_pose
from or52e8_site_scan import read
from reec_pose_checks import read_pose


LIGAND = 'HETATM 9000  O10 LIG     0       1.000   2.000   3.000  1.00  0.00           O  '


class SavedLigandRelaxationTest(unittest.TestCase):
    def test_extract_embedded_poses_ignores_appended_original(self):
        text = ''.join(f'Pose: {n}\nNode: 0\nTotal: -5\nPDBDAT:\n{LIGAND}\nEND\n'
                       for n in (1, 2, 3))
        text += 'Original PDB:\nATOM original receptor\nEND\n'
        poses = extract_poses(text)
        self.assertEqual([p['pose'] for p in poses], [1, 2, 3])
        self.assertTrue(all(len(p['lines']) == 1 for p in poses))
        with self.assertRaises(ValueError):
            extract_poses(f'Pose: 1\nPDBDAT:\n{LIGAND}\n')
        self.assertEqual(extract_poses(f'Pose: 1\nDisqualified because: rejected\nPDBDAT:\n{LIGAND}\nEND\n'), [])

    def test_translated_blank_chain_pose_restores_intact_backbone_and_keeps_ligand(self):
        lines = [line for line in (ROOT / 'pdbs/OR52/OR52E8.active.pdb').read_text().splitlines()
                 if line.startswith('ATOM  ') and int(line[22:26]) in (161, 208)]
        shifted = []
        for line in lines:
            if int(line[22:26]) != 208:
                continue
            xyz = np.array([float(line[i:i+8]) for i in (30, 38, 46)]) + [1, 2, 3]
            shifted.append(line[:21]+' '+line[22:30]+''.join(f'{v:8.3f}' for v in xyz)+line[54:])
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            canonical, output = root / 'source.pdb', root / 'merged.pdb'
            canonical.write_text('\n'.join(lines)+'\nEND\n')
            result = merge_pose(canonical, {'lines': shifted+[LIGAND]}, output)
            old, merged = read(canonical), read(output)
            self.assertEqual(set(old), set(merged))
            for key in old:
                np.testing.assert_allclose(old[key]['xyz'], merged[key]['xyz'], atol=.00001)
            self.assertGreater(result['maximum_saved_backbone_shift'], 3)
            self.assertEqual(read_pose(output)[1]['O10']['point'], (1, 2, 3))


if __name__ == '__main__':
    unittest.main()
