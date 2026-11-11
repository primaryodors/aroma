from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'hm'))
from or52e8_site_scan import fit, rotate


class SiteGeometryTest(unittest.TestCase):
    def test_torsion_rotation_preserves_axis_and_covalent_geometry(self):
        points = np.array([[1., 0, 0], [0, 1, 0], [0, 0, 2], [0, 0, -1]])
        origin, axis = np.zeros(3), np.array([0., 0, 1])
        changed = rotate(points, origin, axis, 90)
        np.testing.assert_allclose(changed[0], [0, 1, 0], atol=1e-12)
        np.testing.assert_allclose(changed[2:], points[2:], atol=1e-12)
        np.testing.assert_allclose(np.linalg.norm(changed[:, None]-changed[None, :], axis=2),
                                   np.linalg.norm(points[:, None]-points[None, :], axis=2), atol=1e-12)
        np.testing.assert_allclose(rotate(changed, origin, axis, -90), points, atol=1e-12)

    def test_local_fit_recovers_rigid_motion_without_mirroring(self):
        moving = np.array([[0., 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
        fixed = moving[:, [1, 2, 0]] + [10, 20, 30]
        transform = fit(moving, fixed)
        np.testing.assert_allclose(transform(moving), fixed, atol=1e-12)
        self.assertGreater(np.linalg.det(transform(moving)[1:]-transform(moving)[0]), 0)


if __name__ == '__main__':
    unittest.main()
