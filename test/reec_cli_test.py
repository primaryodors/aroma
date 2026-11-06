"""End-to-end REEC regression checks using a short receptor fragment."""
import json
import subprocess
import tempfile
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'hm'))
import reec_review
CODES = dict(zip(
    'ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR'.split(),
    'ACDEFGHIKLMNPQRSTVWY'))


class ReecCLI(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='reec-')
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)
        self.template = ROOT / 'testdata/reec-small.pdb'
        self.sequence = ''.join(CODES[s[17:20]] for s in self.template.read_text().splitlines()
                                if s.startswith('ATOM  ') and s[12:16].strip() == 'CA')
        self.ali = self.work / 'model.ali'
        self.ali.write_text(f'>P1;model\nsequence:model::::::::\n{self.sequence}*\n')
        self.output = self.work / 'model.pdb'

    def run_model(self, *extra, template=None):
        result = subprocess.run([str(ROOT / 'bin/reec'), '--alignment', str(self.ali),
                               '--template', str(template or self.template), '--output', str(self.output),
                               '--iterations', '1500', *extra], cwd=ROOT, capture_output=True, text=True)
        if result.returncode == 0:
            review = reec_review.inspect(self.output)
            self.assertTrue(review['valid'], review['failures'])
        return result

    def test_relaxation_and_complete_sequence(self):
        result = self.run_model()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        report = json.loads(Path(str(self.output) + '.json').read_text())
        self.assertTrue(report['valid'])
        self.assertLess(report['final_energy'], report['initial_energy'])
        self.assertLess(report['max_bond_error'], .08)
        self.assertEqual(report['chirality_errors'], 0)
        self.assertEqual(report['severe_clashes'], 0)
        residues = {int(s[22:26]): s[17:20] for s in self.output.read_text().splitlines()
                    if s.startswith('ATOM  ')}
        self.assertEqual(''.join(CODES[residues[r]] for r in sorted(residues)), self.sequence)

    def test_impossible_contact_retains_previous_output(self):
        self.output.write_text('previous model\n')
        contacts = self.work / 'contacts.txt'
        contacts.write_text('contact hard 1 N 1 CA 8 9 100\n')
        result = self.run_model('--contacts', str(contacts))
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        self.assertEqual(self.output.read_text(), 'previous model\n')
        self.assertFalse(json.loads(Path(str(self.output) + '.json').read_text())['valid'])

    def test_fixed_ligand_coordinates_and_contact(self):
        first = next(s for s in self.template.read_text().splitlines()
                     if s.startswith('ATOM  ') and s[12:16].strip() == 'CA')
        xyz = [float(first[i:i+8]) for i in (30, 38, 46)]
        xyz[0] += 4
        environment = self.work / 'ligand.pdb'
        line = (f'HETATM    1  C1  LIG A 999    {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}'
                '  1.00  0.00           C')
        second = (f'HETATM    2  C2  LIG A 999    {xyz[0]+1.5:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}'
                  '  1.00  0.00           C')
        environment.write_text(line + '\n' + second + '\nCONECT    1    2\nEND\n')
        contacts = self.work / 'contacts.txt'
        contacts.write_text('contact hard 1 CA 999 C1 3.8 4.2 100\n')
        result = self.run_model('--environment', str(environment), '--contacts', str(contacts))
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        ligand = next(s for s in self.output.read_text().splitlines() if s.startswith('HETATM'))
        self.assertEqual(ligand[30:54], line[30:54])
        text = self.output.read_text()
        ligands = [s for s in text.splitlines() if s.startswith('HETATM')]
        self.assertEqual(len(ligands), 2)
        self.assertEqual(ligands[1][30:54], second[30:54])
        self.assertIn('CONECT' + ligands[0][6:11] + ligands[1][6:11], text)
        self.assertEqual(json.loads(Path(str(self.output) + '.json').read_text())['mandatory_violations'], 0)

    def test_missing_residue_and_sidechains_are_rebuilt(self):
        damaged = self.work / 'damaged.pdb'
        damaged.write_text('\n'.join(s for s in self.template.read_text().splitlines()
                                      if s.startswith('ATOM  ') and int(s[22:26]) != 9
                                      and s[12:16].strip() in ('N', 'CA', 'C', 'O')) + '\nEND\n')
        columns = list(self.sequence)
        columns[5] = '-'
        self.ali.write_text('>P1;damaged\nstructure:damaged:FIRST:A:LAST:A::::\n'
                            + ''.join(columns) + '*\n>P1;model\nsequence:model::::::::\n'
                            + self.sequence + '*\n')
        result = self.run_model('--iterations', '4500', template=f'damaged={damaged}')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        lines = [s for s in self.output.read_text().splitlines() if s.startswith('ATOM  ')]
        self.assertEqual(len({int(s[22:26]) for s in lines}), len(self.sequence))
        self.assertTrue(any(s[12:16].strip() == 'SG' for s in lines))
        self.assertTrue(any(s[12:16].strip() == 'OXT' for s in lines))

    def test_multiple_templates_and_alternate_starts(self):
        rotated = self.work / 'rotated.pdb'
        lines = []
        for line in self.template.read_text().splitlines():
            if line.startswith('ATOM  '):
                x, y, z = (float(line[i:i+8]) for i in (30, 38, 46))
                line = line[:30] + f'{-y+8:8.3f}{x-4:8.3f}{z+3:8.3f}' + line[54:]
            lines.append(line)
        rotated.write_text('\n'.join(lines) + '\n')
        result = self.run_model('--template', str(rotated), '--models', '2')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn('candidate 2/2', result.stdout)

    def test_atom_alternatives(self):
        contacts = self.work / 'alternatives.txt'
        contacts.write_text('contact hard 1 CA 2 N|CA 1 8 100\n')
        result = self.run_model('--contacts', str(contacts))
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        report = json.loads(Path(str(self.output) + '.json').read_text())
        self.assertEqual(report['mandatory_violations'], 0)
        self.assertIn(report['contacts'][0]['atom2'], ('N', 'CA'))

    def test_ligand_superposition_into_template_frame(self):
        original = self.template.read_text()
        first = next(s for s in original.splitlines() if s.startswith('ATOM  ')
                     and s[12:16].strip() == 'CA')
        x, y, z = (float(first[i:i+8]) for i in (30, 38, 46))
        x += 4
        ligand = f'HETATM    1  C1  LIG A 999    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C'
        environment = self.work / 'environment.pdb'
        environment.write_text(original.replace('END', '') + ligand + '\nEND\n')
        rotated = self.work / 'reference.pdb'
        lines = []
        for line in original.splitlines():
            if line.startswith('ATOM  '):
                a, b, c = (float(line[i:i+8]) for i in (30, 38, 46))
                line = line[:30] + f'{-b+8:8.3f}{a-4:8.3f}{c+3:8.3f}' + line[54:]
            lines.append(line)
        rotated.write_text('\n'.join(lines) + '\n')
        result = self.run_model('--environment', str(environment), '--fit-environment', template=rotated)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        result_line = next(s for s in self.output.read_text().splitlines() if s.startswith('HETATM'))
        actual = [float(result_line[i:i+8]) for i in (30, 38, 46)]
        for a, b in zip(actual, (-y+8, x-4, z+3)):
            self.assertAlmostEqual(a, b, places=2)

    def test_arg_and_pro_mutations_use_standard_geometry(self):
        sequence = list(self.sequence)
        sequence[2] = 'R'
        sequence[4] = 'P'
        self.ali.write_text('>P1;model\nsequence:model::::::::\n' + ''.join(sequence) + '*\n')
        result = self.run_model('--iterations', '4500')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        lines = self.output.read_text().splitlines()
        self.assertTrue(any(s.startswith('ATOM  ') and s[17:20] == 'ARG' for s in lines))
        self.assertTrue(any(s.startswith('ATOM  ') and s[17:20] == 'PRO' for s in lines))

    def test_disulfide_is_covalent_and_serialized(self):
        sequence = list(self.sequence)
        sequence[3] = 'C'
        self.ali.write_text('>P1;model\nsequence:model::::::::\n' + ''.join(sequence) + '*\n')
        contacts = self.work / 'disulfide.txt'
        contacts.write_text('disulfide 1 4\n')
        result = self.run_model('--contacts', str(contacts), '--iterations', '9000', '--models', '2')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        report = json.loads(Path(str(self.output) + '.json').read_text())
        self.assertEqual(report['contacts'][0]['kind'], 'disulfide')
        self.assertAlmostEqual(report['contacts'][0]['distance'], 2.05, delta=.03)
        self.assertIn('CONECT', self.output.read_text())
        loaded = subprocess.run([str(ROOT / 'test/reec_test'), str(self.output)], cwd=ROOT, capture_output=True, text=True)
        self.assertEqual(loaded.returncode, 0, loaded.stdout + loaded.stderr)

    def test_explicit_environment_hydrogens(self):
        first = next(s for s in self.template.read_text().splitlines()
                     if s.startswith('ATOM  ') and s[12:16].strip() == 'CA')
        x,y,z = [float(first[i:i+8]) for i in (30,38,46)]
        environment = self.work / 'environment.pdb'
        environment.write_text(
            f'HETATM    1  C1  LIG A 999    {x+4:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n'
            f'HETATM    2  H2  LIG A 999    {x+4.36:8.3f}{y+1.03:8.3f}{z:8.3f}  1.00  0.00           H\nEND\n')
        contacts = self.work / 'contacts.txt'
        contacts.write_text('contact hard 1 CA 999 H2 3 6 100\n')
        result = self.run_model('--environment', str(environment), '--environment-hydrogens', '--contacts', str(contacts))
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        report = json.loads(Path(str(self.output) + '.json').read_text())
        self.assertEqual(report['fixed_atoms'], 2)
        lines = [s for s in self.output.read_text().splitlines() if s.startswith(('ATOM  ', 'HETATM'))]
        self.assertEqual(report['heavy_atoms'], sum(s[76:78].strip() != 'H' for s in lines))
        self.assertTrue(any(s[12:16].strip() == 'H2' for s in lines))

    def test_invalid_atom_fails_explicitly(self):
        contacts = self.work / 'contacts.txt'
        contacts.write_text('contact hard 1 ABSENT 2 CA 2 3 10\n')
        result = self.run_model('--contacts', str(contacts))
        self.assertEqual(result.returncode, 1)
        self.assertIn('Unknown atom', result.stderr)
        self.assertFalse(self.output.exists())

    def test_unterminated_alignment_is_rejected(self):
        self.ali.write_text('>P1;model\nsequence:model::::::::\nAAA\n')
        result = self.run_model()
        self.assertEqual(result.returncode, 1)
        self.assertIn('Incomplete PIR', result.stderr)

    def test_aligned_template_mismatch_is_rejected(self):
        self.ali.write_text('>P1;tpl\nstructure:tpl:FIRST:A:LAST:A::::\n'
                            + 'A' * len(self.sequence) + '*\n>P1;model\nsequence:model::::::::\n'
                            + self.sequence + '*\n')
        result = self.run_model('--template', f'tpl={self.template}')
        self.assertEqual(result.returncode, 1)
        self.assertIn('sequence mismatch', result.stderr)


if __name__ == '__main__':
    unittest.main()
