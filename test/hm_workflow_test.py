"""Both numerical engines execute the original scripts' shared workflow tails."""
from contextlib import ExitStack, redirect_stdout, redirect_stderr
import io
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'hm'))
import dohm
import fixfail
import modeling_backend as engines
import reec_backend as backend
import data.globals
import data.odorutils


class SharedWorkflow(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        backend.pu.load_prots()
        data.odorutils.load_odors()

    def exercise(self, module, backend_name, actions=(), fail_post=False, trusted=False,
                 reject_model=0, reject_post=0, max_attempts=None, fatal_model=False,
                 stale_dock=False, reject_reset=False):
        with tempfile.TemporaryDirectory(prefix='aroma-workflow-') as temporary:
            root = Path(temporary)
            hm = root / 'hm'
            for folder in (hm, root / 'tmp', root / 'data', root / 'pdbs/OR1', root / 'out/OR1/OR1A1'):
                folder.mkdir(parents=True, exist_ok=True)
            original = ROOT / 'pdbs/OR1/OR1A1.inactive.pdb'
            for state in ('active', 'inactive'):
                shutil.copyfile(original, root / f'pdbs/OR1/OR1A1.{state}.pdb')
            (root / 'data/cavopts.json').write_text('[]')
            source = root / 'template.pdb'
            shutil.copyfile(original, source)
            alignment = root / 'template.ali'
            sequence = backend.pu.prots['OR1A1']['sequence']
            alignment.write_text(f'>P1;tpl\nstructure:tpl:FIRST:A:LAST:A::::\n{sequence}*\n'
                                 f'>P1;OR1A1\nsequence:OR1A1:1:A::::::\n{sequence}*\n')
            events, restraint_sets, scripts = [], [], []
            builds = []
            seeds = []
            destination = root / 'result.pdb'
            destination.write_text('previous model\n')
            raw_template = original.read_text()
            if trusted:
                validated_source = root / 'pdbs/_validated/OR1A1.active.pdb'
                validated_source.parent.mkdir()
                validated_source.write_text(raw_template)

            class Model(engines.ReecModel):
                def make(self):
                    assert self.max_var_iterations == (7 * 2**events.count('model') if module is dohm else 7)
                    assert self.ending_model == 1
                    self.restraints = engines.Restraints()
                    self.disulfides = []
                    self.special_patches(None)
                    self.special_restraints(None)
                    restraint_sets.append(self.contact_lines())
                    events.append('model')
                    seeds.append(self.env.options.seed)
                    if fatal_model:
                        raise ValueError('Invalid modeling input')
                    if events.count('model') <= reject_model:
                        events.append('rejected')
                        report = root / 'tmp' / f'rejection-{events.count("model")}.json'
                        report.write_text('{"valid": false}')
                        raise backend.CandidateRejected('Test candidate rejected', report)
                    output = Path('raw-model.pdb')
                    output.write_text(raw_template)
                    alternate = Path('selected-model.pdb')
                    alternate.write_text(raw_template)
                    self.outputs[:] = [{'name': output.name, 'molpdf': 2, 'failure': None},
                                      {'name': alternate.name, 'molpdf': 1, 'failure': None}]

            def api(args):
                self.assertFalse(trusted, 'Validated models bypass backend construction')
                selected = args.backend
                args.backend = 'reec'
                result = engines.api(args)
                args.backend = selected
                result.AutoModel = result.DOPEHRLoopModel = Model
                return result

            def run(command, **kwargs):
                command = [str(value) for value in command]
                name = Path(command[0]).name
                if name == 'make' and command[1:] == ['bin/reec']:
                    self.assertFalse(events)
                    self.assertFalse(scripts)
                    builds.append(command)
                if name == 'phew':
                    events.append('phew')
                    text = Path(command[1]).read_text()
                    scripts.append(text)
                    if fail_post:
                        destination.write_text('partial output\n')
                        raise subprocess.CalledProcessError(1, command)
                    import re
                    output = re.search(r'LET \$outf = "([^"]+)"', text)
                    if output:
                        Path(output[1]).write_text(raw_template+'REMARK 999 SHARED PHEW\n')
                elif name == 'ic':
                    events.append('ic')
                    path = Path(command[1])
                    path.write_text(path.read_text()+'REMARK 998 SHARED IC\n')
                elif name == 'poslig':
                    shutil.copyfile(command[1], command[2])
                elif name == 'cavity_fit':
                    shutil.copyfile(root / 'pdbs/OR1/OR1A1.active.pdb', command[command.index('-o')+1])
                    events.append('fit')
                elif name == 'cavity_search':
                    Path(command[command.index('-o')+1]).write_text('test cavity\n')
                elif name == 'bash':
                    events.append('dock')
                    count = events.count('model')
                    dock = root / 'out/OR1/OR1A1/OR1A1~linalool.active.dock'
                    if not stale_dock:
                        dock.write_text(f'Total: -10\n{4 if count >= 2 else 0} pose(s) found.\n')
                elif name.startswith('python'):
                    self.assertIn('--backend', command)
                    self.assertEqual(command[command.index('--backend')+1], backend_name)
                    events.append('reset')
                    if reject_reset:
                        raise subprocess.CalledProcessError(2, command)
                elif name != 'make':
                    raise AssertionError(command)
                return subprocess.CompletedProcess(command, 0)

            def finalize(model, path, candidate=None, *, target=None):
                self.assertEqual(target, 'OR1A1')
                if trusted:
                    self.assertIsNone(model)
                    self.assertEqual(candidate, validated_source)
                else:
                    if module is dohm and events.count('finalize'):
                        self.assertIn(candidate, ('selected-model.pdb', 'raw-model.pdb'))
                    else:
                        self.assertEqual(candidate, 'selected-model.pdb')
                events.append('finalize')
                text = Path(path).read_text()
                self.assertIn('SHARED PHEW', text)
                self.assertIn('SHARED IC', text)
                if events.count('finalize') <= reject_post:
                    events.append('rejected')
                    raise backend.CandidateRejected('Test postprocessing rejected')

            argv = ['OR1A1', '--backend', backend_name, '--models', '2', '--iterations', '7',
                    '--skip-reshape', '--template', f'tpl={source}', '--alignment', str(alignment)]
            if module is dohm:
                argv += ['nodel', '--output', str(destination)]
            elif actions:
                docked = root / 'out/OR1/OR1A1/OR1A1~linalool.active.model1.pdb'
                docked.write_text(raw_template)
                argv.insert(1, str(docked))
                argv += list(actions)
                if stale_dock:
                    (docked.parent / 'OR1A1~linalool.active.dock').write_text('Total: -10\n4 pose(s) found.\n')
            else:
                argv += ['--output', str(destination)]
            if max_attempts is not None:
                argv += ['--max-attempts', str(max_attempts)]
            with ExitStack() as stack:
                stack.enter_context(backend.directory(root))
                stack.enter_context(patch.object(sys, 'argv', ['workflow-test']))
                stack.enter_context(patch.object(module, '__file__', str(hm / module.__name__)+'.py'))
                stack.enter_context(patch.object(module, 'script_dir', str(hm)))
                stack.enter_context(patch.object(module, 'root_dir', str(root)))
                stack.enter_context(patch.object(module, 'backend_api', side_effect=api))
                stack.enter_context(patch.object(module, 'finalize', side_effect=finalize))
                stack.enter_context(patch.object(module.subprocess, 'run', side_effect=run))
                stack.enter_context(patch.object(backend.pu, 'load_prots'))
                stack.enter_context(patch.object(data.odorutils, 'load_odors'))
                stack.enter_context(patch.object(data.globals, 'wait_cool_cpu'))
                stack.enter_context(patch.object(data.odorutils, 'ensure_sdf_exists', return_value=str(root / 'ligand.sdf')))
                stack.enter_context(redirect_stdout(io.StringIO()))
                stack.enter_context(redirect_stderr(io.StringIO()))
                if fail_post:
                    with self.assertRaises(subprocess.CalledProcessError):
                        module.main(argv)
                    self.assertEqual(destination.read_text(), 'previous model\n')
                elif fatal_model:
                    with self.assertRaisesRegex(ValueError, 'Invalid modeling input'):
                        module.main(argv)
                    self.assertEqual(events.count('model'), 1)
                    self.assertNotIn('dock', events)
                    return events, restraint_sets
                elif module is dohm and max_attempts is not None and reject_model >= max_attempts:
                    with self.assertRaises(backend.CandidateRejected):
                        module.main(argv)
                    self.assertEqual(events.count('model'), max_attempts)
                    self.assertEqual(destination.read_text(), 'previous model\n')
                    return events, restraint_sets
                else:
                    status = module.main(argv)
                    self.assertEqual(len(builds), int(backend_name == 'reec'))
                    if module is dohm and not trusted:
                        logs = list((root / 'tmp').glob('OR1A1-dohm-*/attempts.json'))
                        self.assertEqual(len(logs), 1)
                        records = json.loads(logs[0].read_text())
                        self.assertEqual(records[-1]['status'], 'accepted')
                        self.assertEqual(sorted({entry['seed'] for entry in records}), seeds)
                        self.assertEqual(sorted({entry['iterations'] for entry in records}),
                                         [7*2**i for i in range(len(seeds))])
                    if module is fixfail:
                        logs = list((root / 'tmp').glob('OR1A1-fixfail-*/attempts.json'))
                        self.assertEqual(len(logs), 1)
                        records = json.loads(logs[0].read_text())
                        self.assertEqual([entry['seed'] for entry in records], seeds)
                        for entry in records:
                            if 'report' in entry:
                                self.assertTrue(Path(entry['report']).is_file())
                        if max_attempts is not None:
                            self.assertLessEqual(len(records), max_attempts)
                        if status == 2:
                            self.assertEqual(status, 2)
                            self.assertIn(records[-1]['status'], ('rejected', 'dock_failed'))
                            if records[-1]['status'] == 'rejected':
                                self.assertNotIn('dock', events)
                                if actions:
                                    self.assertEqual((root / 'pdbs/OR1/OR1A1.active.pdb').read_text(), raw_template)
                                else:
                                    self.assertEqual(destination.read_text(), 'previous model\n')
                            return events, restraint_sets
                        self.assertEqual(status, 0)
                        self.assertEqual(records[-1]['status'], 'success' if actions else 'accepted')
            if trusted:
                self.assertEqual(validated_source.read_text(), raw_template)
                self.assertNotIn('model', events)
                self.assertIn(str(validated_source), scripts[0])
                self.assertIn('UPRIGHT', scripts[0])
                self.assertIn('BWCENTER', scripts[0])
            else:
                self.assertTrue(restraint_sets)
                self.assertTrue(any(line.startswith('disulfide ') for line in restraint_sets[0]))
                self.assertTrue(any(' CA ' in line for line in restraint_sets[0]))
            self.assertTrue(scripts)
            self.assertIn('BWCOPY I A', scripts[0])
            if not fail_post:
                self.assertEqual(events.count('model') - min(reject_model, events.count('model'))
                                 + int(trusted), events.count('finalize')
                                 - (min(reject_post, events.count('finalize')) if module is dohm else 0))
                self.assertGreaterEqual(events.count('phew'), events.count('model') - reject_model)
            return events, restraint_sets

    def test_dohm_shared_tail_and_restraints(self):
        results = [self.exercise(dohm, name) for name in ('reec', 'modeller')]
        self.assertEqual(results[0], results[1])

    def test_dohm_retries_generation_with_larger_budget_and_fresh_seed(self):
        events, _ = self.exercise(dohm, 'reec', reject_model=2)
        self.assertEqual(events.count('model'), 3)
        self.assertEqual(events.count('finalize'), 1)

    def test_dohm_tries_other_candidates_after_processing_rejection(self):
        events, _ = self.exercise(dohm, 'reec', reject_post=1)
        self.assertEqual(events.count('model'), 1)
        self.assertEqual(events.count('finalize'), 2)

    def test_dohm_explicit_attempt_limit_preserves_previous_output(self):
        self.exercise(dohm, 'reec', reject_model=10, max_attempts=2)

    def test_dohm_input_error_stops_immediately(self):
        self.exercise(dohm, 'reec', fatal_model=True)

    def test_loop_retries_rejection_without_docking_stale_model(self):
        for phase in ('model', 'post'):
            kwargs = {'reject_model': 1} if phase == 'model' else {'reject_post': 1}
            events, _ = self.exercise(fixfail, 'reec', ('loop',), **kwargs)
            self.assertEqual(events.count('model'), 2)
            self.assertEqual(events.count('dock'), 1)
            self.assertLess(events.index('rejected'), events.index('dock'))

    def test_rejection_limit_and_single_attempt_exit_cleanly(self):
        self.exercise(fixfail, 'reec', ('loop',), reject_model=10, max_attempts=2)
        self.exercise(fixfail, 'reec', ('loop',), reject_post=10, max_attempts=2)
        self.exercise(fixfail, 'reec', reject_model=1)

    def test_input_error_is_fatal_in_loop_mode(self):
        self.exercise(fixfail, 'reec', ('loop',), fatal_model=True)

    def test_stale_dock_success_is_not_accepted(self):
        events, _ = self.exercise(fixfail, 'reec', ('loop',), stale_dock=True, max_attempts=2)
        self.assertEqual(events.count('model'), 2)
        self.assertEqual(events.count('dock'), 2)

    def test_rejected_reset_continues_loop_without_bootstrap_dock(self):
        events, _ = self.exercise(fixfail, 'reec', ('reset', 'loop'), reject_reset=True)
        self.assertLess(events.index('model'), events.index('dock'))

    def test_validated_model_uses_shared_processing_without_modeling(self):
        for name in ('reec', 'modeller'):
            events, _ = self.exercise(dohm, name, trusted=True)
            self.assertIn('phew', events)
            self.assertIn('ic', events)
            self.assertIn('finalize', events)

    def test_validated_processing_failure_restores_previous_output(self):
        self.exercise(dohm, 'reec', fail_post=True, trusted=True)

    def test_fixfail_shared_tail_and_restraints(self):
        results = [self.exercise(fixfail, name) for name in ('reec', 'modeller')]
        self.assertEqual(results[0], results[1])

    def test_fixfail_reset_fit_predock_and_retries(self):
        for name in ('reec', 'modeller'):
            events, _ = self.exercise(fixfail, name, ('reset', 'fit', 'predock', 'loop'))
            self.assertEqual(events.count('model'), 2)
            self.assertIn('reset', events)
            self.assertIn('fit', events)
            self.assertEqual(events.count('dock'), 3)

    def test_shared_postprocessing_failure_preserves_previous_output(self):
        for module in (dohm, fixfail):
            for name in ('reec', 'modeller'):
                self.exercise(module, name, fail_post=True)

    def test_reec_api_does_not_import_modeller(self):
        args = engines.options(['OR1A1', '--backend', 'reec'])
        with patch.object(engines.importlib, 'import_module') as imported:
            engines.api(args)
        imported.assert_not_called()

    def test_complete_reec_workflows_with_modeller_unavailable(self):
        with patch.dict(sys.modules, {'modeller': None, 'modeller.automodel': None}):
            for module in (dohm, fixfail):
                self.exercise(module, 'reec')
            self.exercise(dohm, 'reec', trusted=True)

    def test_build_failure_stops_before_workflow_initialization(self):
        for module in (dohm, fixfail):
            with patch.object(backend, 'command', side_effect=subprocess.CalledProcessError(2, ['make'])), \
                    patch.object(module, 'backend_api') as api:
                with self.assertRaises(subprocess.CalledProcessError):
                    module.main(['OR1A1', '--backend', 'reec'])
                api.assert_not_called()

    def test_prepare_builds_once_and_only_for_reec(self):
        for name in ('reec', 'modeller'):
            args = engines.options(['OR1A1', '--backend', name])
            with patch.object(backend, 'command') as command:
                engines.prepare(args)
                engines.prepare(args)
                self.assertEqual(command.call_count, int(name == 'reec'))
                if name == 'reec':
                    command.assert_called_once_with(['make', 'bin/reec'])


if __name__ == '__main__':
    unittest.main()
