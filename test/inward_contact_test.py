"""Receptor contact preparation preserves heavy ligand and covalent geometry."""
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'hm'))
from modeling_backend import preserve_environment_labels, inward_contact_failures


def atoms(path):
    return {(line[17:20], int(line[22:26]), line[12:16].strip()):
            tuple(float(line[i:i+8]) for i in (30, 38, 46))
            for line in path.read_text().splitlines()
            if line.startswith(('ATOM  ', 'HETATM'))}


class InwardContactTest(unittest.TestCase):
    def test_automatic_distance_uses_heavy_default_and_typed_hydrogen_radius(self):
        for selector, expected in [('208:OD1', 'Contact restraint O10 208:OD1 2.8'),
                                   ('208:2HD2', 'Contact donor O10 208:2HD2 1.66')]:
            with self.subTest(selector=selector), tempfile.TemporaryDirectory() as folder:
                output = Path(folder) / 'automatic.pdb'
                result = subprocess.run([str(ROOT / 'bin/poslig'), str(ROOT / 'testdata/asparagine_inward_contact.pdb'),
                                         str(output), '--receptor', 'auto', 'O10', selector],
                                        cwd=ROOT, capture_output=True, text=True, timeout=30)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertIn(expected, result.stdout)

    def test_extent_hydrogen_retains_its_distance_with_a_heavy_proxy(self):
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder) / 'arginine.pdb'
            result = subprocess.run([str(ROOT / 'bin/poslig'), str(ROOT / 'testdata/asparagine_inward_contact.pdb'),
                                     str(output), '--receptor', '2.5', 'O9', '268:EXTENT'],
                                    cwd=ROOT, capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            checks = [line.split() for line in result.stdout.splitlines() if line.startswith('Contact donor ')]
            proxies = [line.split() for line in result.stdout.splitlines() if line.startswith('Contact restraint ')]
            self.assertEqual(len(checks), 1)
            self.assertEqual(float(checks[0][4]), 2.5)
            atom = checks[0][3].split(':')[1]
            self.assertTrue(atom.lstrip('0123456789').startswith('H'))
            self.assertFalse(proxies[0][3].split(':')[1].lstrip('0123456789').startswith('H'))
            positions = atoms(output)
            oxygen = next(point for key, point in positions.items() if key[0] == 'LIG' and key[2] == 'O9')
            self.assertAlmostEqual(math.dist(positions['ARG', 268, atom], oxygen), 2.5, delta=.01)

    def test_extent_resolves_to_an_actual_protein_atom(self):
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder) / 'extent.pdb'
            result = subprocess.run([str(ROOT / 'bin/poslig'), str(ROOT / 'testdata/asparagine_inward_contact.pdb'),
                                     str(output), '--receptor', '2.94', 'O10', '208:EXTENT'],
                                    cwd=ROOT, capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertIn('Contact restraint O10 208:', result.stdout)
            self.assertNotIn(':EXTENT', result.stdout)

    def test_contact_review_rejects_poor_donor_direction(self):
        lines = (ROOT / 'testdata/asparagine_inward_contact.pdb').read_text().splitlines()
        n = next(line for line in lines if line.startswith('ATOM  ') and int(line[22:26]) == 208 and line[12:16].strip() == 'ND2')
        h = next(line for line in lines if line.startswith('ATOM  ') and int(line[22:26]) == 208 and line[12:16].strip() == '2HD2')
        o = next(line for line in lines if line.startswith('HETATM') and line[12:16].strip() == 'O10')
        o = o[:21] + 'A 999' + o[26:]
        def positioned(line, xyz):
            return line[:30] + ''.join(f'{v:8.3f}' for v in xyz) + line[54:]
        checks = [('O10', '208:ND2', 2.66), ('O10', '208:2HD2', 1.66)]
        with tempfile.TemporaryDirectory() as folder:
            model = Path(folder) / 'contact.pdb'
            model.write_text('\n'.join([positioned(n,(2.66,0,0)), positioned(h,(1.66,0,0)), positioned(o,(0,0,0)), 'END'])+'\n')
            self.assertEqual(inward_contact_failures(model, checks), [])
            x = (4 + 2.32**2 - 1) / (2*2.32)
            model.write_text('\n'.join([positioned(n,(2.32,0,0)), positioned(h,(x,math.sqrt(4-x*x),0)), positioned(o,(0,0,0)), 'END'])+'\n')
            self.assertTrue(any('direction' in failure for failure in inward_contact_failures(model, checks)))

    def test_environment_labels_preserve_postprocessed_coordinates(self):
        with tempfile.TemporaryDirectory() as folder:
            source, output = (Path(folder) / name for name in ('candidate.pdb', 'output.pdb'))
            line = next(line for line in (ROOT / 'testdata/asparagine_inward_contact.pdb').read_text().splitlines()
                        if line.startswith('HETATM') and line[12:16].strip() == 'O10')
            source.write_text(line[:21] + 'A 999' + line[26:] + '\nEND\n')
            shifted = line[:30] + f'{15.:8.3f}{20.:8.3f}{25.:8.3f}' + line[54:]
            output.write_text(shifted + '\nEND\n')
            preserve_environment_labels(source, output)
            actual = output.read_text().splitlines()[0]
            self.assertEqual(actual[30:54], shifted[30:54])
            self.assertEqual(int(actual[22:26]), 999)
            self.assertTrue(inward_contact_failures(output, [('O10', '208:ND2', 2.66)]))
            output.write_text('END\n')
            with self.assertRaises(ValueError):
                preserve_environment_labels(source, output)

    def test_inward_asparagine_keeps_ligand_and_backbone_fixed(self):
        source = ROOT / 'testdata/asparagine_inward_contact.pdb'
        before = atoms(source)
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder) / 'inward.pdb'
            result = subprocess.run([str(ROOT / 'bin/poslig'), str(source), str(output),
                                     '--receptor', '2.66', 'O10', '208:ND2'],
                                    cwd=ROOT, capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            after = atoms(output)
            self.assertTrue(all(line[21] == 'A' for line in output.read_text().splitlines()
                                if line.startswith('ATOM  ')))
            for key, point in before.items():
                if (key[0] == 'LIG' and not key[2].startswith('H')) or key[2] in ('N', 'CA', 'C', 'O'):
                    self.assertEqual(after[key], point)
            ligand = next(point for key, point in after.items() if key[0] == 'LIG' and key[2] == 'O10')
            donor = after['ASN', 208, 'ND2']
            self.assertAlmostEqual(math.dist(donor, ligand), 2.66, delta=.3)
            self.assertIn('Contact donor O10 208:', result.stdout)
            self.assertIn(' 1.66', result.stdout)
            hydrogen = min((after['ASN', 208, name] for name in ('1HD2', '2HD2')),
                           key=lambda point: math.dist(point, ligand))
            self.assertAlmostEqual(math.dist(ligand, hydrogen), 1.66, delta=.3)
            a = tuple(ligand[i] - hydrogen[i] for i in range(3))
            b = tuple(donor[i] - hydrogen[i] for i in range(3))
            angle = math.degrees(math.acos(max(-1, min(1, sum(x*y for x, y in zip(a, b)) /
                                                        (math.sqrt(sum(x*x for x in a))*math.sqrt(sum(x*x for x in b)))))))
            self.assertGreater(angle, 120)
            for a, b in [('CA', 'CB'), ('CB', 'CG'), ('CG', 'OD1'), ('CG', 'ND2')]:
                pair = [('ASN', 208, name) for name in (a, b)]
                self.assertAlmostEqual(math.dist(*(before[key] for key in pair)),
                                       math.dist(*(after[key] for key in pair)), delta=.003)


if __name__ == '__main__':
    unittest.main()
