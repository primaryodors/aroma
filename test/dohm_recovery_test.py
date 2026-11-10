import json
import hashlib
import os
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'hm'))
import dohm
import dohm_recovery as recovery
from reec_backend import CandidateRejected, parser


class BatchRecoveryTest(unittest.TestCase):
    def test_default_batch_budget_is_bounded_and_repair_budget_is_preserved(self):
        self.assertEqual(parser().parse_args(['OR52E8']).max_attempts, 3)
        self.assertEqual(parser(True).parse_args(['OR52E8']).max_attempts, 100)
        self.assertEqual(parser().parse_args(['OR52E8', '--max-attempts', '0']).max_attempts, 0)

    def test_rejected_native_candidate_is_written_and_labeled(self):
        dohm.pu.load_prots()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source, output = root / 'candidate.pdb', root / 'active.pdb'
            source.write_bytes((ROOT / 'testdata/reec-small.pdb').read_bytes())
            Path(str(source)+'.json').write_text(json.dumps({'valid': False, 'contacts': []}))
            review = {'valid': True, 'failures': [], 'residues': 12}
            with patch.object(recovery, 'review_source', return_value=review):
                recovery.fallback_model('OR52E8', output, SimpleNamespace(last_rejected_model=source),
                                        CandidateRejected('Test rejection'), root)
            self.assertIn('DOHM FALLBACK', output.read_text())
            report = json.loads(Path(str(output)+'.json').read_text())
            self.assertFalse(report['valid'])
            self.assertFalse(report['refresh_accepted'])
            self.assertEqual(report['status'], 'fallback')

    def test_batch_entry_point_recovers_after_candidate_rejection(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / 'active.pdb'
            output.write_bytes((ROOT / 'testdata/reec-small.pdb').read_bytes())
            with patch.object(dohm, '_main', side_effect=CandidateRejected('Test rejection')), \
                    patch.object(dohm, 'root_dir', str(root)):
                dohm.main(['OR52E8', '--output', str(output)])
            self.assertIn('ATOM  ', output.read_text())
            self.assertEqual(json.loads(Path(str(output)+'.json').read_text())['status'], 'fallback')

    def test_invalid_validated_source_is_repaired_before_publication(self):
        old_argv, old_cwd = sys.argv[:], Path.cwd()
        try:
            with tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                (root / 'hm').mkdir()
                trusted = root / 'pdbs/_validated/OR52E8.active.pdb'
                trusted.parent.mkdir(parents=True)
                trusted.write_bytes((ROOT / 'testdata/reec-small.pdb').read_bytes())
                repaired, output = root / 'repaired.pdb', root / 'active.pdb'
                repaired.write_text('REMARK repaired fixture\n'+trusted.read_text())
                Path(str(repaired)+'.json').write_text(json.dumps({'valid': True, 'contacts': []}))
                with patch.object(dohm, 'root_dir', str(root)), patch.object(dohm, 'script_dir', str(root / 'hm')), \
                        patch.object(dohm, 'prepare'), patch.object(dohm, 'process_model') as process, \
                        patch.object(recovery, 'review_source', return_value={'valid': False, 'failures': ['Test geometry']}), \
                        patch.object(recovery, 'repair_existing', return_value=(repaired, root)) as repair:
                    dohm.main(['OR52E8', '--output', str(output)])
                process.assert_not_called()
                repair.assert_called_once()
                self.assertEqual(output.read_bytes(), repaired.read_bytes())
                self.assertEqual(trusted.read_bytes(), repaired.read_bytes())
                self.assertTrue(json.loads(Path(str(trusted)+'.json').read_text())['valid'])
        finally:
            sys.argv = old_argv
            os.chdir(old_cwd)

    def test_ready_validated_source_keeps_its_conformation_without_reprocessing(self):
        old_argv, old_cwd = sys.argv[:], Path.cwd()
        try:
            with tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                (root / 'hm').mkdir()
                trusted = root / 'pdbs/_validated/OR52E8.active.pdb'
                trusted.parent.mkdir(parents=True)
                trusted.write_bytes((ROOT / 'testdata/reec-small.pdb').read_bytes())
                status = {'valid': True, 'validated_ready': True,
                          'validated_sha256': hashlib.sha256(trusted.read_bytes()).hexdigest()}
                Path(str(trusted)+'.json').write_text(json.dumps(status))
                output = root / 'active.pdb'
                with patch.object(dohm, 'root_dir', str(root)), patch.object(dohm, 'script_dir', str(root / 'hm')), \
                        patch.object(dohm, 'prepare'), patch.object(dohm, 'process_model') as process, \
                        patch.object(recovery, 'review_source', return_value={'valid': True, 'failures': []}), \
                        patch.object(recovery, 'repair_existing') as repair:
                    dohm.main(['OR52E8', '--output', str(output)])
                process.assert_not_called()
                repair.assert_not_called()
                self.assertEqual(output.read_bytes(), trusted.read_bytes())
        finally:
            sys.argv = old_argv
            os.chdir(old_cwd)


if __name__ == '__main__':
    unittest.main()
