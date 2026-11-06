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

    def test_fixfail_preserves_coordinates_omitted_by_dock_export(self):
        import math
        codes = dict(zip(
            'ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR'.split(),
            'ACDEFGHIKLMNPQRSTVWY'))
        lines = [line for line in (ROOT / 'testdata/reec-small.pdb').read_text().splitlines()
                 if line.startswith('ATOM  ')]
        numbers = sorted({int(line[22:26]) for line in lines})
        numbering = {old: i + 1 for i, old in enumerate(numbers)}
        lines = [line[:22] + f'{numbering[int(line[22:26])]:4d}' + line[26:] for line in lines]
        sequence = ''.join(codes[line[17:20]] for line in lines if line[12:16].strip() == 'CA')
        omitted = {1, 2, len(sequence) // 2}
        with tempfile.TemporaryDirectory(prefix='reec-clipped-repair-') as temporary:
            root = Path(temporary)
            (root / 'tmp').mkdir()
            raw = root / 'raw.pdb'
            raw.write_text('\n'.join(lines + ['TER', 'END']) + '\n')
            alignment = root / 'reference.ali'
            alignment.write_text(f'>P1;OR1A1\nsequence:OR1A1:1:A::::::\n{sequence}*\n')
            normalized = root / 'normalized.pdb'
            subprocess.run([str(ROOT / 'bin/reec'), '--alignment', str(alignment),
                            '--template', str(raw), '--target', 'OR1A1', '--output', str(normalized),
                            '--iterations', '1500', '--models', '1'], cwd=ROOT,
                           capture_output=True, text=True, check=True)
            lines = [line for line in normalized.read_text().splitlines() if line.startswith('ATOM  ')]
            family = root / 'pdbs' / 'OR1'
            family.mkdir(parents=True)
            original = family / 'OR1A1.active.pdb'
            original.write_text('\n'.join(lines + ['TER', 'END']) + '\n')
            clipped = root / 'clipped.pdb'
            clipped.write_text('\n'.join([line for line in lines if int(line[22:26]) not in omitted]
                                         + ['TER', 'END']) + '\n')
            destination = root / 'repaired.pdb'
            baseline = root / 'full-repair.pdb'
            messages = []
            def command(args):
                if str(args[0]) == 'make':
                    return
                actual = [str(item) for item in args]
                actual[0] = str(ROOT / 'bin' / Path(actual[0]).name)
                result = subprocess.run(actual, cwd=ROOT, capture_output=True, text=True, check=True)
                messages.append(result.stdout)
            entry = dict(backend.pu.prots['OR1A1'], sequence=sequence)
            with patch.object(backend, 'ROOT', root), \
                    patch.object(backend, 'receptor', return_value=('OR1A1', 'active')), \
                    patch.dict(backend.pu.prots, {'OR1A1': entry}), \
                    patch.object(backend, 'receptor_contacts', return_value=[]), \
                    patch.object(backend, 'command', side_effect=command):
                backend.run_fixfail(['OR1A1', str(original), '--skip-reshape',
                                     '--iterations', '1500', '--models', '1', '--output', str(baseline)])
                backend.run_fixfail(['OR1A1', str(clipped), '--skip-reshape',
                                     '--iterations', '1500', '--models', '1', '--output', str(destination)])
            reference_ca = {int(line[22:26]): tuple(float(line[i:i+8]) for i in (30, 38, 46))
                            for line in baseline.read_text().splitlines()
                            if line.startswith('ATOM  ') and line[12:16].strip() == 'CA'}
            restored = {int(line[22:26]): tuple(float(line[i:i+8]) for i in (30, 38, 46))
                        for line in destination.read_text().splitlines()
                        if line.startswith('ATOM  ') and line[12:16].strip() == 'CA'}
            for residue in omitted:
                self.assertLess(math.dist(reference_ca[residue], restored[residue]), .5,
                                f'Residue {residue}: ' + '\n'.join(messages))

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
