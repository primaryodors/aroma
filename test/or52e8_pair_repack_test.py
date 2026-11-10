from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'hm'))
from or52e8_pair_repack import SideChain, read, validate


class PairRepackGeometryTest(unittest.TestCase):
    def test_joint_rotations_preserve_all_covalent_geometry(self):
        protein = read(ROOT / 'pdbs/OR52/OR52E8.active.pdb')
        asn, met = SideChain(protein, 208), SideChain(protein, 161)
        altered = {k: dict(v) for k, v in protein.items()}
        for chain, angles in ((asn, (85, 65)), (met, (-60, 120, 30))):
            for key, point in zip(chain.keys, chain.conformer(angles)):
                altered[key]['xyz'] = point
        report = validate(protein, altered, (asn, met))
        self.assertTrue(report['all_other_atoms_unchanged'])
        self.assertLess(report['max_bond_change'], 1e-10)
        self.assertLess(report['max_bond_angle_change'], 1e-8)

    def test_clash_mask_excludes_bonds_but_retains_nonbonded_pairs(self):
        protein = read(ROOT / 'pdbs/OR52/OR52E8.active.pdb')
        chain = SideChain(protein, 161)
        self.assertFalse(chain.self_mask[chain.names.index('CB'), chain.names.index('CG')])
        self.assertFalse(chain.self_mask[chain.names.index('CB'), chain.names.index('SD')])
        self.assertTrue(chain.self_mask[chain.names.index('CB'), chain.names.index('CE')])
        np.testing.assert_allclose(chain.conformer((0, 0, 0)), chain.start, atol=1e-12)


if __name__ == '__main__':
    unittest.main()
