"""Protection and exit-status checks for isolated workflow measurements."""
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'hm'))
import reec_workflow_baseline as workflow


class WorkflowIsolation(unittest.TestCase):
    def test_outputs_cannot_escape_through_paths_or_symlinks(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            workspace = root / 'workspace'
            workspace.mkdir()
            (workspace / 'escape').symlink_to(root)
            self.assertEqual(workflow.isolated_output(workspace, workspace / 'model.pdb'),
                             workspace / 'model.pdb')
            for output in (root / 'model.pdb', workspace / 'escape/model.pdb'):
                with self.assertRaises(ValueError):
                    workflow.isolated_output(workspace, output)

    def test_repair_rejection_status_is_preserved(self):
        class Rejected(ValueError):
            pass
        with tempfile.TemporaryDirectory() as temporary:
            workspace = Path(temporary)
            (workspace / 'results').mkdir()
            (workspace / 'snapshot.json').write_text('{}\n')
            module = SimpleNamespace(main=Mock(return_value=2))
            backend = SimpleNamespace(CandidateRejected=Rejected)
            arguments = ['baseline', 'run', '--workspace', str(workspace), '--entry', 'fixfail', '--', 'OR1B1']
            with patch.object(sys, 'argv', arguments), patch.object(workflow, 'bind_workspace', return_value=(module, backend)):
                self.assertEqual(workflow.main(), 2)
            module.main.assert_called_once_with(['OR1B1'])

    def test_pdb_and_binary_snapshots_are_independent_files(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source, workspace = root / 'source', root / 'workspace'
            for name in ('data', 'src', 'hm/tpl', 'bin', 'obj', 'sdf', 'pdbs/OR1'):
                (source / name).mkdir(parents=True, exist_ok=True)
            for name in ('dohm', 'fixfail', 'modeling_backend', 'reec_backend'):
                (source / 'hm' / f'{name}.py').write_text('pass\n')
            for name in ('reec', 'aromadock'):
                (source / 'bin' / name).write_bytes(b'original binary')
            (source / 'src/reec.cpp').write_text('// fixture\n')
            (source / 'makefile').write_text('# fixture\n')
            original = source / 'pdbs/OR1/OR1A1.active.pdb'
            original.write_text('original model\n')
            with patch.object(workflow, 'ROOT', source):
                workflow.prepare_workspace(workspace)
            self.assertFalse((workspace / 'data').is_symlink())
            previous_directory = Path.cwd()
            try:
                import os
                os.chdir(workspace / 'data')
                os.chdir('..')
                self.assertEqual(Path.cwd(), workspace)
            finally:
                os.chdir(previous_directory)
            clone = workspace / 'pdbs/OR1/OR1A1.active.pdb'
            self.assertFalse(clone.is_symlink())
            clone.write_text('changed benchmark model\n')
            (workspace / 'bin/reec').write_bytes(b'changed benchmark binary')
            self.assertEqual(original.read_text(), 'original model\n')
            self.assertEqual((source / 'bin/reec').read_bytes(), b'original binary')
            manifest = json.loads((workspace / 'snapshot.json').read_text())
            self.assertEqual(manifest['canonical_pdb_hashes']['pdbs/OR1/OR1A1.active.pdb'],
                             workflow.digest(original))


if __name__ == '__main__':
    unittest.main()
