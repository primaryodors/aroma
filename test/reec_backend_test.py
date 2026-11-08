"""Adapter contract checks; imports and help work without MODELLER."""
import os
import shutil
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

    def test_or5v1_primary_preparation_keeps_experimental_templates_separate(self):
        from dohm import primary_experimental_alignment
        from modeling_backend import pir_rows
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            (work / 'tpl').symlink_to(ROOT / 'hm/tpl')
            shutil.copy2(ROOT / 'hm/experimental.ali', work / 'experimental.ali')
            self.assertIsNone(primary_experimental_alignment('OR2AT4', work))
            alignment, templates, provenance = primary_experimental_alignment('OR5V1', work)
            self.assertEqual([value.split('=')[0] for value in templates], ['9lkb', '9wpm'])
            rows = pir_rows(alignment)
            self.assertEqual(rows['OR5V1'][1].replace('-', ''), backend.pu.prots['OR5V1']['sequence'])
            codes = dict(zip('ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR'.split(),
                             'ACDEFGHIKLMNPQRSTVWY'))
            for name in ('9lkb', '9wpm'):
                coordinates = ''.join(codes[line[17:20]] for line in (work / f'{name}.pdb').read_text().splitlines()
                                      if line.startswith('ATOM  ') and line[12:16].strip() == 'CA')
                self.assertEqual(rows[name][1].replace('-', ''), coordinates)
            self.assertEqual(provenance, '9lkb 9wpm')

    def test_exit_two_is_recoverable_and_other_exit_codes_are_fatal(self):
        from types import SimpleNamespace
        args = SimpleNamespace(iterations=1, models=1, contacts=[], chain='A', seed=9)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            work = root / 'work'
            work.mkdir()
            (root / 'tmp').mkdir()
            report = Path(str(work / 'result.pdb') + '.json')
            report.write_text('{"valid": false, "contacts": []}')
            for code in (2, 1):
                error = subprocess.CalledProcessError(code, ['reec'])
                with patch.object(backend, 'ROOT', root), patch.object(backend, 'command', side_effect=error) as command:
                    expected = backend.CandidateRejected if code == 2 else subprocess.CalledProcessError
                    with self.assertRaises(expected) as raised:
                        backend.model('OR1A1', args, work, [], work / 'input.ali', [])
                    if code == 2:
                        self.assertEqual(raised.exception.report, report)
                    actual = command.call_args.args[0]
                    self.assertEqual(actual[actual.index('--seed')+1], 9)
                    self.assertEqual(actual[actual.index('--anneal-steps')+1], 0)
                    self.assertEqual((root / 'tmp/OR1A1.reec-failure.json').read_text(), report.read_text())

    def test_annealing_defaults_off_for_generation_and_repair(self):
        for repair in (False, True):
            parser = backend.parser(repair)
            self.assertEqual(parser.parse_args(['OR1B1']).anneal_steps, 0)
            self.assertEqual(parser.parse_args(['OR1B1', '--anneal-steps', '64']).anneal_steps, 64)

    def test_geometry_rejection_report_survives_publication_rollback(self):
        import modeling_backend as engines
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / 'active.pdb'
            output.write_text('previous PDB')
            sidecar = Path(str(output)+'.json')
            sidecar.write_text('previous report')
            review = {'valid': False, 'failures': ['Test geometry failure'], 'residues': 1}
            with patch.object(backend, 'ROOT', root), patch('reec_review.inspect', return_value=review):
                with self.assertRaises(backend.CandidateRejected) as raised:
                    with engines.publication_guard(output):
                        output.write_text('END\n')
                        engines.finalize(None, output)
            self.assertEqual(output.read_text(), 'previous PDB')
            self.assertEqual(sidecar.read_text(), 'previous report')
            report = raised.exception.report.read_text()
            self.assertIn('Test geometry failure', report)
            self.assertIn('"valid": false', report)

    def test_backend_selection(self):
        with patch.dict(os.environ, {}, clear=True):
            self.assertTrue(backend.requested([], default='reec'))
            self.assertFalse(backend.requested([]))
            for args in (['--backend', 'modeller'], ['--backend=modeller']):
                self.assertFalse(backend.requested(args, default='reec'))
        with patch.dict(os.environ, {'AROMA_HM_BACKEND': 'reec'}):
            self.assertTrue(backend.requested([]))
            self.assertFalse(backend.requested(['--backend', 'modeller']))
        with patch.dict(os.environ, {'AROMA_HM_BACKEND': 'modeller'}):
            self.assertFalse(backend.requested([], default='reec'))
            self.assertTrue(backend.requested(['--backend=reec'], default='reec'))
        self.assertTrue(backend.requested(['--backend=reec']))
        with self.assertRaises(ValueError):
            backend.requested(['--backend', 'invalid'])

    def test_modeller_backend_options_leave_legacy_arguments_intact(self):
        for args in (['--backend', 'modeller', 'OR5K2', 'nodel'],
                     ['OR5K2', '--backend=modeller', 'nodel']):
            self.assertEqual(backend.legacy_backend_args(args), ['OR5K2', 'nodel'])

    def test_dohm_help_defaults_to_reec(self):
        environment = dict(os.environ)
        environment.pop('AROMA_HM_BACKEND', None)
        result = subprocess.run([sys.executable, str(ROOT / 'hm/dohm.py'), '--help'],
                                env=environment, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('Receptor homology modeling', result.stdout)
        self.assertIn('REEC or MODELLER', result.stdout)

    def test_validated_or5k_alignments(self):
        codes = dict(zip(
            'ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR'.split(),
            'ACDEFGHIKLMNPQRSTVWY'))
        source = ROOT / 'pdbs/_validated/OR5K1.active.pdb'
        template_sequence = ''.join(codes[line[17:20]] for line in source.read_text().splitlines()
                                    if line.startswith('ATOM  ') and line[12:16].strip() == 'CA')
        with tempfile.TemporaryDirectory(prefix='reec-validated-') as temporary:
            work = Path(temporary)
            for receptor in ('OR5K1', 'OR5K2', 'OR5K3', 'OR5K4'):
                with self.subTest(receptor=receptor):
                    alignment, templates, provenance = backend.validated_alignment(receptor, work)
                    rows = [s[:-1] for s in alignment.read_text().splitlines() if s.endswith('*')]
                    self.assertEqual(len(rows[0]), len(rows[1]))
                    self.assertEqual(rows[0].replace('-', ''), template_sequence)
                    self.assertEqual(rows[1].replace('-', ''), backend.pu.prots[receptor]['sequence'])
                    self.assertEqual(templates, [f'{receptor}_tpl={work / (receptor + "_tpl.pdb")}'])
                    self.assertEqual((work / f'{receptor}_tpl.pdb').read_bytes(), source.read_bytes())
                    self.assertEqual(provenance, 'pdbs/_validated/OR5K1.active.pdb')
            self.assertIsNone(backend.validated_alignment('OR5A1', work))
            self.assertIsNone(backend.validated_alignment('OR5K20P', work))

    def test_validated_or5k_template_is_mandatory(self):
        with tempfile.TemporaryDirectory(prefix='reec-validated-missing-') as temporary:
            with patch.object(backend, 'ROOT', Path(temporary)):
                with self.assertRaisesRegex(ValueError, 'Validated OR5K1 template is missing'):
                    backend.validated_alignment('OR5K2', Path(temporary))

    def test_entrypoints_help_without_modeller(self):
        for name in ('dohm', 'fixfail'):
            result = subprocess.run([sys.executable, str(ROOT / 'hm' / f'{name}.py'),
                                     '--backend', 'reec', '--help'], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn('REEC', result.stdout)

    def test_workflow_compatibility_calls_original_scripts(self):
        import dohm, fixfail
        for function, module in ((backend.run_dohm, dohm), (backend.run_fixfail, fixfail)):
            with patch.object(module, 'main', return_value=17) as shared:
                self.assertEqual(function(['OR1A1']), 17)
                shared.assert_called_once_with(['OR1A1'])

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
            import modeling_backend as engines
            with patch.object(backend, 'ROOT', root), \
                    patch.dict(backend.pu.prots, {'OR1A1': entry}), \
                    patch.object(backend, 'command', side_effect=command), backend.directory(root):
                for source, target in ((original, baseline), (clipped, destination)):
                    shutil.copyfile(source, root / 'repair.pdb')
                    args = engines.options(['OR1A1', '--iterations', '1500', '--models', '1'])
                    env = engines.Environment(args)
                    env.io.atom_files_directory = [str(root)]
                    model = engines.construct(engines.ReecModel, env, alnfile=alignment,
                                              knowns='repair', sequence='OR1A1', references=[original])
                    model.max_var_iterations, model.ending_model = 1500, 0
                    model.make()
                    shutil.copyfile(root / model.outputs[0]['name'], target)
            reference_ca = {int(line[22:26]): tuple(float(line[i:i+8]) for i in (30, 38, 46))
                            for line in baseline.read_text().splitlines()
                            if line.startswith('ATOM  ') and line[12:16].strip() == 'CA'}
            restored = {int(line[22:26]): tuple(float(line[i:i+8]) for i in (30, 38, 46))
                        for line in destination.read_text().splitlines()
                        if line.startswith('ATOM  ') and line[12:16].strip() == 'CA'}
            for residue in omitted:
                self.assertLess(math.dist(reference_ca[residue], restored[residue]), .5,
                                f'Residue {residue}: ' + '\n'.join(messages))

    def test_shared_alignment_rows_pad_trailing_gaps(self):
        import modeling_backend as engines
        with tempfile.TemporaryDirectory(prefix='reec-adapter-') as temporary:
            work = Path(temporary)
            source = work / 'source.ali'
            source.write_text('>P1;tpl\nstructure:tpl:FIRST:A:LAST:A::::\n--AAA*\n'
                              '>P1;OR1A1\nsequence:OR1A1::::::::\n--AAA----*\n')
            destination = engines.protein_alignment(source, ['tpl'], 'OR1A1', work / 'protein.ali')
            rows = [line for line in destination.read_text().splitlines() if line.endswith('*')]
            self.assertEqual(rows, ['--AAA----*', '--AAA----*'])

    def test_ligand_repair_preserves_environment_alignment_frame(self):
        import modeling_backend as engines
        codes = dict(zip(
            'ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR'.split(),
            'ACDEFGHIKLMNPQRSTVWY'))
        protein = [line for line in (ROOT / 'testdata/reec-small.pdb').read_text().splitlines()
                   if line.startswith('ATOM  ')]
        sequence = ''.join(codes[line[17:20]] for line in protein if line[12:16].strip() == 'CA')
        first = next(line for line in protein if line[12:16].strip() == 'CA')
        x, y, z = (float(first[i:i+8]) for i in (30, 38, 46))
        ligand = [f'HETATM 9001  C1  LIG A 999    {x+4:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C',
                  f'HETATM 9002  C2  LIG A 999    {x+5.5:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C']
        with tempfile.TemporaryDirectory(prefix='reec-ligand-repair-') as temporary:
            work = Path(temporary)
            (work / 'tmp').mkdir()
            template = work / 'repair.pdb'
            template.write_text('\n'.join(protein + ligand + ['CONECT 9001 9002', 'END']) + '\n')
            alignment = work / 'repair.ali'
            alignment.write_text(f'>P1;repair\nstructure:repair:FIRST:A:LAST:A::::\n{sequence}.*\n'
                                 f'>P1;fragment\nsequence:fragment::::::::\n{sequence}.*\n')
            captured = []

            def command(args):
                if str(args[0]) == 'make':
                    return
                actual = [str(item) for item in args]
                environment = Path(actual[actual.index('--environment') + 1])
                captured.append(environment.read_text())
                actual[0] = str(ROOT / 'bin/reec')
                subprocess.run(actual, cwd=ROOT, capture_output=True, text=True, check=True)

            args = engines.options(['OR1A2', '--iterations', '1500', '--models', '1'])
            env = engines.Environment(args)
            env.io.atom_files_directory = [str(work)]
            env.io.hetatm = True
            with patch.object(backend, 'ROOT', work), \
                    patch.object(backend, 'command', side_effect=command), backend.directory(work):
                model = engines.ReecModel(env, alignment, 'repair', 'fragment')
                model.max_var_iterations, model.ending_model = 1500, 0
                model.make()
                output = (work / model.outputs[0]['name']).read_text()
            self.assertTrue(any(line.startswith('ATOM  ') for line in captured[0].splitlines()))
            atoms = [line for line in output.splitlines() if line.startswith('HETATM')]
            self.assertEqual(len(atoms), 2)
            for actual, original in zip(atoms, ligand):
                self.assertEqual(actual[30:54], original[30:54])
            self.assertIn('CONECT' + atoms[0][6:11] + atoms[1][6:11], output)


if __name__ == '__main__':
    unittest.main()
