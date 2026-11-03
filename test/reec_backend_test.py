"""Adapter contract checks; imports and help work without MODELLER."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'hm'))
import reec_backend as backend


class ReecBackend(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        backend.pu.load_prots()

    def test_backend_selection(self):
        with patch.dict(os.environ, {'AROMA_HM_BACKEND': 'reec'}):
            self.assertTrue(backend.requested([]))
            self.assertFalse(backend.requested(['--backend', 'modeller']))
        self.assertTrue(backend.requested(['--backend=reec']))
        with self.assertRaises(ValueError):
            backend.requested(['--backend', 'invalid'])

    def test_entrypoints_help_without_modeller(self):
        for name in ('dohm', 'fixfail'):
            result = subprocess.run([sys.executable, str(ROOT / 'hm' / f'{name}.py'),
                                     '--backend', 'reec', '--help'], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn('REEC', result.stdout)

    def test_conserved_contacts(self):
        contacts = backend.receptor_contacts('OR1A1')
        self.assertIn('disulfide 97 179', contacts)
        self.assertIn('disulfide 169 189', contacts)
        self.assertIn('contact hard 258 OH 180 OD1|OD2 2.5 3.2 100', contacts)
        contacts = backend.receptor_contacts('OR51E2')
        self.assertIn('disulfide 4 168', contacts)
        self.assertIn('contact hard 150 NE|NH1|NH2 209 OD1|OD2 2.6 3.6 100', contacts)

    def test_metadata_contains_bw_and_helices(self):
        lines = backend.metadata('OR51E2')
        self.assertIn('REMARK 800 SITE BW 4.50 148', lines)
        self.assertIn('REMARK 650 HELIX TMR4 137 163', lines)

    def test_rigid_receptor_frame_centers_bw_anchors(self):
        import math
        rcp = 'OR51E2'
        points = {n: (n % 7, n % 11, n % 5) for n in range(1, len(backend.pu.prots[rcp]['sequence']) + 1)}
        for helix in range(1, 8):
            region = backend.pu.prots[rcp]['region'][f'TMR{helix}']
            for n in range(region['start'], region['end'] + 1):
                u = (n - region['start']) / (region['end'] - region['start'])
                points[n] = (10*math.cos(helix), (1 if helix % 2 else -1)*(15-30*u),
                             10*math.sin(helix))
        text = ''.join(f'ATOM  {n:5d}  CA  ALA A{n:4d}    {p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}  1.00  0.00           C\n'
                       for n, p in points.items())
        result = backend.receptor_frame(rcp, text)
        ca = {int(s[22:26]): tuple(float(s[i:i+8]) for i in (30, 38, 46)) for s in result.splitlines()}
        anchors = [ca[backend.pu.resno_from_bw(rcp, f'{h}.50')] for h in range(1, 8)]
        for axis in range(3):
            self.assertAlmostEqual(sum(p[axis] for p in anchors)/7, 0, places=3)
        self.assertAlmostEqual(math.dist(points[30], points[210]), math.dist(ca[30], ca[210]), delta=.003)

    def test_publication_rejects_invalid_candidates(self):
        import json
        with tempfile.TemporaryDirectory(prefix='reec-publication-') as temporary:
            work = Path(temporary)
            model = work / 'candidate.pdb'
            model.write_text((ROOT / 'testdata/reec-small.pdb').read_text())
            destination = work / 'published.pdb'
            destination.write_text('previous model\n')
            for validity in (False, True):
                with self.subTest(reported_valid=validity):
                    Path(str(model) + '.json').write_text(json.dumps({'valid': validity, 'contacts': []}))
                    with self.assertRaises(ValueError):
                        backend.publish(model, destination, [])
                    self.assertEqual(destination.read_text(), 'previous model\n')
            self.assertFalse(list(work.glob('*.reec-*')))

    def test_legacy_options_can_follow_backend_selection(self):
        args = backend.parser().parse_intermixed_args(['OR1A1', '--backend', 'reec', 'nodel'])
        self.assertEqual(args.legacy, ['nodel'])

    def test_composite_rows_pad_trailing_gaps(self):
        with tempfile.TemporaryDirectory(prefix='reec-adapter-') as temporary:
            work = Path(temporary)
            (work / 'allgpcr.ali').write_text('>P1;OR1A1\nsequence:OR1A1::::::::\n--AAA----*\n')
            with patch.object(backend, 'HM', work), patch.object(backend, 'command'), \
                    patch.object(backend.pu, 'custom_pdb_template', return_value='tpl\n--AAA*'):
                alignment, _, _ = backend.composite_alignment('OR1A1', work)
            rows = [s for s in alignment.read_text().splitlines() if s.endswith('*')]
            self.assertEqual(rows, ['--AAA----*', '--AAA----*'])


if __name__ == '__main__':
    unittest.main()
