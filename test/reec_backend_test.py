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
