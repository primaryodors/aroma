"""End-to-end REEC regression checks using a short receptor fragment."""
import json
import re
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
ALIASES = {
    '--alignment': '-a', '--template': '-t', '--output': '-o', '--target': '-T',
    '--contacts': '-c', '--chain': '-C', '--iterations': '-i', '--models': '-m',
    '--report': '-r', '--environment': '-e', '--fit-environment': '-f',
    '--environment-hydrogens': '-H', '--candidate-prefix': '-p', '--seed': '-s',
    '--anneal-steps': '-A', '--anneal-temperature': '-u', '--anneal-residues': '-R',
    '--clash-scale': '-x', '--primary-template': '-P', '--nofail': '-N', '--help': '-h',
    '--backbone-steps': '-b', '--backbone-range': '-B', '--backbone-angle': '-v',
    '--fallback-iterations': '-j', '--fallback-method': '-J',
}


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

    def test_help_documents_unique_single_letter_aliases_for_every_option(self):
        outputs = []
        for flag in ('--help', '-h'):
            result = subprocess.run([str(ROOT / 'bin/reec'), flag], cwd=ROOT,
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            outputs.append(result.stdout)
        self.assertEqual(outputs[0], outputs[1])
        entries = re.findall(r'^  (-[A-Za-z]), (--[a-z-]+)', outputs[0], re.M)
        self.assertEqual({long: short for short, long in entries}, ALIASES)
        self.assertEqual(len({short for short, _ in entries}), len(entries))

    def test_backbone_sampling_is_reproducible_and_bounded(self):
        environment = self.work / 'fixed.pdb'
        environment.write_text('HETATM    1  C1  LIG A 999      40.000  40.000  40.000  1.00  0.00           C\nEND\n')
        command = ('--iterations','3000','--backbone-steps','4',
                   '--backbone-range','2:9','--backbone-angle','2','--seed','7',
                   '--environment',str(environment))
        result = self.run_model(*command)
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        pdb, report = self.output.read_bytes(), Path(str(self.output)+'.json').read_bytes()
        sampling = json.loads(report)['backbone_sampling']
        fixed = next(line for line in self.output.read_text().splitlines() if line.startswith('HETATM'))
        self.assertEqual([float(fixed[i:i+8]) for i in (30,38,46)], [40.,40.,40.])
        self.assertEqual(sampling['trials'], 4)
        self.assertLessEqual(sampling['closure_bond_error'], .08)
        self.assertLessEqual(sampling['max_backbone_displacement'], .75)
        result = self.run_model(*command)
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        self.assertEqual(self.output.read_bytes(), pdb)
        self.assertEqual(Path(str(self.output)+'.json').read_bytes(), report)

    def test_backbone_sampling_rejects_invalid_ranges_and_amplitudes(self):
        for extra in (('--backbone-steps','4'),
                      ('--backbone-steps','4','--backbone-range','2:30'),
                      ('--backbone-steps','4','--backbone-range','2:9','--backbone-angle','6')):
            result = self.run_model(*extra)
            self.assertEqual(result.returncode, 1)
            self.assertFalse(self.output.exists())

    def test_minimization_stop_reasons_are_recorded(self):
        result = self.run_model('--iterations','30')
        self.assertIn(result.returncode, (0,2))
        report = json.loads(Path(str(self.output)+'.json').read_text())
        self.assertEqual(len(report['minimizations']), 3)
        for stage in report['minimizations']:
            self.assertIn(stage['reason'], ('iteration_limit','converged','small_step_stall','line_search_stall'))
            self.assertLessEqual(stage['completed'], stage['budget'])
            self.assertGreaterEqual(stage['max_gradient'], 0)
            self.assertLessEqual(stage['score_after'], stage['score_before'])

    def test_optional_fallback_respects_the_total_iteration_cap(self):
        result = self.run_model('--iterations','1500','--fallback-iterations','300','--fallback-method','cg')
        self.assertIn(result.returncode, (0,2))
        report = json.loads(Path(str(self.output)+'.json').read_text())
        self.assertEqual([s['budget'] for s in report['minimizations']], [500,500,200,300])
        self.assertEqual(report['minimizations'][-1]['method'],
                         'conjugate_gradient' if report['fallback']['eligible'] else 'lbfgs')

    def test_primary_template_alias_and_missing_id(self):
        self.ali.write_text(f'>P1;donor\nstructure:donor:FIRST:A:LAST:A::::\n{self.sequence}*\n'
                            f'>P1;model\nsequence:model::::::::\n{self.sequence}*\n')
        for flag in ('--primary-template', '-P'):
            result = self.run_model(flag, 'donor', template=f'donor={self.template}')
            self.assertEqual(result.returncode, 0, result.stderr)
            report = json.loads(Path(str(self.output)+'.json').read_text())
            self.assertEqual(report['primary_template'], 'donor')
            self.assertEqual(report['template_weighting']['version'], 1)
            self.assertEqual(set(report['template_weighting']['templates'][0]['residue_weights']), {1})
        result = self.run_model('-P', 'absent', template=f'donor={self.template}')
        self.assertEqual(result.returncode, 1)
        self.assertIn('Primary template ID absent', result.stderr)

    def test_all_short_options_match_long_option_models_and_reports(self):
        contacts = self.work / 'contacts.txt'
        contacts.write_text('contact soft 1 CA 2 CA 1 8 10\n')
        original = self.template.read_text()
        first = next(s for s in original.splitlines() if s.startswith('ATOM  ') and s[12:16].strip() == 'CA')
        x, y, z = (float(first[i:i+8]) for i in (30, 38, 46))
        environment = self.work / 'environment.pdb'
        environment.write_text(original.replace('END', '') +
            f'HETATM    1  C1  LIG A 999    {x+100:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n' +
            f'HETATM    2  H2  LIG A 999    {x+100.36:8.3f}{y+1.03:8.3f}{z:8.3f}  1.00  0.00           H\nEND\n')
        report = self.work / 'custom-report.json'
        prefix = self.work / 'export'
        options = [('--alignment', self.ali), ('--template', self.template),
                   ('--template', self.template), ('--output', self.output), ('--target', 'model'),
                   ('--contacts', contacts), ('--contacts', contacts), ('--chain', 'A'),
                   ('--iterations', 1500), ('--models', 2), ('--report', report),
                   ('--environment', environment), ('--fit-environment', None),
                   ('--environment-hydrogens', None), ('--candidate-prefix', prefix),
                   ('--seed', 17), ('--anneal-steps', 4), ('--anneal-temperature', 7),
                   ('--anneal-residues', '1,2,3'), ('--clash-scale', 1.5)]
        artifacts = [self.output, report] + [Path(f'{prefix}.{i}.pdb.json') for i in range(2)]
        snapshots = []
        for short in (False, True):
            for path in artifacts:
                path.unlink(missing_ok=True)
            command = [str(ROOT / 'bin/reec')]
            for flag, value in options:
                command.append(ALIASES[flag] if short else flag)
                if value is not None:
                    command.append(str(value))
            result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            data = json.loads(report.read_text())
            self.assertTrue(data['valid'])
            self.assertEqual(data['fixed_atoms'], 2)
            self.assertEqual(len(data['contacts']), 2)
            self.assertEqual(data['annealing']['trials'], 4)
            Path(str(self.output)+'.json').write_text(report.read_text())
            self.assertTrue(reec_review.inspect(self.output)['valid'])
            snapshots.append([path.read_bytes() for path in artifacts])
        self.assertEqual(snapshots[0], snapshots[1])

    def test_relaxation_and_complete_sequence(self):
        result = self.run_model()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        report = json.loads(Path(str(self.output) + '.json').read_text())
        self.assertTrue(report['valid'])
        self.assertEqual(report['annealing']['trials'], 0)
        self.assertLess(report['final_energy'], report['initial_energy'])
        self.assertLess(report['max_bond_error'], .08)
        self.assertEqual(report['chirality_errors'], 0)
        self.assertEqual(report['severe_clashes'], 0)
        residues = {int(s[22:26]): s[17:20] for s in self.output.read_text().splitlines()
                    if s.startswith('ATOM  ')}
        self.assertEqual(''.join(CODES[residues[r]] for r in sorted(residues)), self.sequence)

    def test_candidate_export_stays_out_of_receptor_storage(self):
        folder = self.work / 'pdbs' / 'OR1'
        folder.mkdir(parents=True)
        result = self.run_model('--candidate-prefix', str(folder / 'candidate'))
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        self.assertEqual(list(folder.iterdir()), [])
        self.assertIn('Candidate export skipped', result.stderr)

    def test_candidate_export(self):
        prefix = self.work / 'candidate'
        result = self.run_model('--models', '2', '--candidate-prefix', str(prefix))
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        scores = []
        for index in range(2):
            path = Path(f'{prefix}.{index}.pdb')
            report = json.loads(Path(str(path)+'.json').read_text())
            if report['valid']:
                self.assertTrue(reec_review.inspect(path)['valid'])
                scores.append(report['final_energy'])
        best = json.loads(Path(str(self.output)+'.json').read_text())
        self.assertTrue(scores)
        self.assertAlmostEqual(min(scores), best['final_energy'])

    def test_ten_candidates_have_distinct_starts(self):
        prefix = self.work / 'distinct'
        result = self.run_model('--models', '10', '--iterations', '0', '--anneal-steps', '0',
                                '--candidate-prefix', str(prefix))
        self.assertIn(result.returncode, (0, 2), result.stdout + result.stderr)
        reports = [json.loads(Path(f'{prefix}.{i}.pdb.json').read_text()) for i in range(10)]
        self.assertEqual(len({report['start_signature'] for report in reports}), 10)
        self.assertEqual([report['candidate_index'] for report in reports], list(range(10)))

    def test_cli_seed_is_reproducible(self):
        energies = []
        for seed in (17, 17, 18):
            result = self.run_model('--seed', str(seed), '--iterations', '0')
            self.assertIn(result.returncode, (0, 2), result.stdout + result.stderr)
            report = json.loads(Path(str(self.output)+'.json').read_text())
            self.assertEqual(report['seed'], seed)
            energies.append(report['final_energy'])
        self.assertEqual(energies[0], energies[1])
        self.assertNotEqual(energies[0], energies[2])

    def test_targeted_annealing_and_disable_switch(self):
        residues = ','.join(str(r) for r in range(1, len(self.sequence)+1))
        for steps in (0, 64):
            result = self.run_model('--anneal-steps', str(steps), '--anneal-residues', residues)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            report = json.loads(Path(str(self.output)+'.json').read_text())['annealing']
            self.assertEqual(report['trials'], steps)
            self.assertLessEqual(report['score_after'], report['score_before'])
            if steps:
                self.assertGreater(report['torsions'], 0)

    def test_annealing_controls_reject_invalid_values(self):
        for flag, value in (('--anneal-steps', '-1'), ('--anneal-temperature', '0'), ('--anneal-residues', '99999')):
            result = self.run_model(flag, value)
            self.assertEqual(result.returncode, 1)

    def test_impossible_contact_retains_previous_output(self):
        self.output.write_text('previous model\n')
        contacts = self.work / 'contacts.txt'
        contacts.write_text('contact hard 1 N 1 CA 8 9 100\n')
        result = self.run_model('--contacts', str(contacts))
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        self.assertEqual(self.output.read_text(), 'previous model\n')
        self.assertFalse(json.loads(Path(str(self.output) + '.json').read_text())['valid'])

    def test_nofail_exports_rejected_model_without_changing_status(self):
        contacts = self.work / 'impossible.contacts'
        contacts.write_text('contact hard 1 N 1 CA 8 9 100\n')
        result = self.run_model('--contacts', str(contacts), '--nofail')
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        self.assertIn('ATOM  ', self.output.read_text())
        self.assertFalse(json.loads(Path(str(self.output)+'.json').read_text())['valid'])

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
        result = self.run_model('--environment', str(environment), '--contacts', str(contacts), '--seed', '19')
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

    def test_dihedral_and_mobile_restraints(self):
        import math
        atoms = {}
        numbering = {}
        for line in self.template.read_text().splitlines():
            if not line.startswith('ATOM  '):
                continue
            number = int(line[22:26])
            numbering.setdefault(number, len(numbering)+1)
            atoms[numbering[number], line[12:16].strip()] = tuple(float(line[i:i+8]) for i in (30,38,46))
        angle = math.radians(reec_review.dihedral(atoms[1,'N'], atoms[1,'CA'], atoms[1,'C'], atoms[2,'N']))
        contacts = self.work / 'restraints.txt'
        contacts.write_text(f'dihedral 1 N 1 CA 1 C 2 N {angle} 0.5\nmobile 1 3\n')
        result = self.run_model('--contacts', str(contacts))
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        report = json.loads(Path(str(self.output)+'.json').read_text())
        self.assertEqual(report['dihedral_restraints'], 1)
        self.assertEqual(report['mobile_ranges'], 1)

    def test_shared_phew_and_ic_postprocessing(self):
        import modeling_backend as engines
        result = self.run_model()
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        processed = self.work / 'postprocessed.pdb'
        script = self.work / 'post.phew'
        script.write_text(f'LOAD {self.output} A A\nHYDRO\nSAVE {processed}\n')
        subprocess.run([str(ROOT/'bin/phew'), str(script)], cwd=ROOT,
                       check=True, capture_output=True, text=True)
        subprocess.run([str(ROOT/'bin/ic'), str(processed), 'minc', 'save'], cwd=ROOT,
                       check=True, capture_output=True, text=True)
        model = object.__new__(engines.ReecModel)
        model.outputs = [{'name': str(self.output)}]
        engines.finalize(model, processed)
        report = json.loads(Path(str(processed)+'.json').read_text())
        self.assertTrue(report['valid'])
        self.assertTrue(report['serialized_geometry']['valid'])

    def test_unterminated_alignment_is_rejected(self):
        self.ali.write_text('>P1;model\nsequence:model::::::::\nAAA\n')
        result = self.run_model()
        self.assertEqual(result.returncode, 1)
        self.assertIn('Incomplete PIR', result.stderr)

    def test_truncated_receptor_does_not_require_terminal_oxygen(self):
        result = self.run_model()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        lines = self.output.read_text().splitlines()
        terminal = len(self.sequence) - 1
        fragment = [s for s in lines if s.startswith('ATOM  ') and int(s[22:26]) <= terminal]
        self.output.write_text('\n'.join(fragment) + '\nEND\n')
        report_path = Path(str(self.output) + '.json')
        report = json.loads(report_path.read_text())
        report['target_sequence_length'] = len(self.sequence)
        report_path.write_text(json.dumps(report))
        review = reec_review.inspect(self.output)
        self.assertTrue(review['valid'], review['failures'])
        report['target_sequence_length'] = terminal
        report_path.write_text(json.dumps(report))
        self.assertTrue(any('OXT' in error for error in reec_review.inspect(self.output)['failures']))
        report['target_sequence_length'] = len(self.sequence)
        report_path.write_text(json.dumps(report))
        fragment = [s for s in fragment if not (int(s[22:26]) == terminal and s[12:16].strip() == 'O')]
        self.output.write_text('\n'.join(fragment) + '\nEND\n')
        self.assertFalse(reec_review.inspect(self.output)['valid'])

    def test_repository_pir_padding_and_unrelated_ligand_row(self):
        self.ali.write_text(
            '>P1;tpl\nstructure:tpl:FIRST:A:LAST:A::::\n'
            + self.sequence + '*\n-----*\n'
            + '>P1;ligand\nstructure:ligand::::::::\n....*\n'
            + '>P1;model\nsequence:model::::::::\n' + self.sequence + '-----*\n')
        result = self.run_model('--target', 'model', template=f'tpl={self.template}')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_selected_ligand_columns_are_rejected(self):
        self.ali.write_text('>P1;model\nsequence:model::::::::\n' + self.sequence + '.*\n')
        result = self.run_model()
        self.assertEqual(result.returncode, 1)
        self.assertIn('omit ligand columns', result.stderr)

    def test_bare_template_id_resolves_repository_template(self):
        result = subprocess.run([
            str(ROOT / 'bin/reec'), '--alignment', str(ROOT / 'hm/allgpcr.ali'),
            '--template', '8f76', '--target', 'OR51E2', '--iterations', '0',
            '--output', str(self.output)], cwd=ROOT, capture_output=True, text=True)
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        self.assertNotIn('Cannot read template', result.stderr)
        self.assertIn('REEC built', result.stdout)

    def test_aligned_template_mismatch_is_rejected(self):
        self.ali.write_text('>P1;tpl\nstructure:tpl:FIRST:A:LAST:A::::\n'
                            + 'A' * len(self.sequence) + '*\n>P1;model\nsequence:model::::::::\n'
                            + self.sequence + '*\n')
        result = self.run_model('--template', f'tpl={self.template}')
        self.assertEqual(result.returncode, 1)
        self.assertIn('sequence mismatch', result.stderr)


if __name__ == '__main__':
    unittest.main()
