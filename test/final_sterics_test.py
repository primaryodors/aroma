"""Publication checks for final-coordinate sterics and bounded cleanup."""
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'hm'))
import modeling_backend as backend
import reec_backend
from reec_review import inspect_sterics


class FinalStericsTest(unittest.TestCase):
    def test_rejection_preserves_prior_files_and_final_diagnostic(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output, candidate = root / 'active.pdb', root / 'candidate.pdb'
            output.write_text('previous model')
            Path(str(output)+'.json').write_text('previous report')
            candidate.write_bytes((ROOT / 'testdata/reec-small.pdb').read_bytes())
            Path(str(candidate)+'.json').write_text(json.dumps({'contacts': [], 'severe_clashes': 0}))
            geometry = {'valid': True, 'failures': [], 'residues': 12}
            sterics = {'severe_clashes': 1, 'max_overlap_angstrom': 1.07202}
            with patch.object(reec_backend, 'ROOT', root), \
                    patch('reec_review.inspect', return_value=geometry), \
                    patch('reec_review.inspect_sterics', return_value=sterics):
                with self.assertRaises(reec_backend.CandidateRejected) as caught:
                    with backend.publication_guard(output):
                        shutil.copyfile(candidate, output)
                        backend.finalize(None, output, candidate)
            self.assertEqual(output.read_text(), 'previous model')
            self.assertEqual(Path(str(output)+'.json').read_text(), 'previous report')
            report = json.loads(caught.exception.report.read_text())
            self.assertFalse(report['valid'])
            self.assertEqual(report['native_severe_clashes'], 0)
            self.assertEqual(report['severe_clashes'], 1)
            self.assertEqual(report['serialized_sterics'], sterics)

    def test_success_records_final_measurement(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / 'active.pdb'
            output.write_bytes((ROOT / 'testdata/reec-small.pdb').read_bytes())
            Path(str(output)+'.json').write_text(json.dumps({'contacts': [], 'severe_clashes': 3}))
            with patch('reec_review.inspect', return_value={'valid': True, 'failures': []}), \
                    patch('reec_review.inspect_sterics', return_value={'severe_clashes': 0, 'max_overlap_angstrom': .8}):
                backend.finalize(None, output, output)
            report = json.loads(Path(str(output)+'.json').read_text())
            self.assertTrue(report['valid'])
            self.assertEqual(report['native_severe_clashes'], 3)
            self.assertEqual(report['severe_clashes'], 0)

    def test_audit_error_cannot_publish(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / 'active.pdb'
            output.write_bytes((ROOT / 'testdata/reec-small.pdb').read_bytes())
            Path(str(output)+'.json').write_text(json.dumps({'contacts': []}))
            with patch.object(reec_backend, 'ROOT', root), \
                    patch('reec_review.inspect', return_value={'valid': True, 'failures': []}), \
                    patch('reec_review.inspect_sterics', side_effect=ValueError('Audit unavailable')):
                with self.assertRaises(reec_backend.CandidateRejected):
                    backend.finalize(None, output, output)

    def test_failed_cleanup_restores_prior_publication(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output, candidate = root / 'active.pdb', root / 'candidate.pdb'
            output.write_text('previous model')
            Path(str(output)+'.json').write_text('previous report')
            candidate.write_bytes((ROOT / 'testdata/reec-small.pdb').read_bytes())
            contacts = [{'kind': 'contact', 'hard': True, 'res1': 4, 'atom1': 'CA',
                         'res2': 5, 'atom2': 'CA', 'min': 0, 'max': 100}]
            Path(str(candidate)+'.json').write_text(json.dumps({'contacts': contacts, 'seed': 7}))
            reec_backend.pu.load_prots()
            model = object.__new__(backend.ReecModel)
            with patch('reec_review.inspect', return_value={'valid': True, 'failures': []}), \
                    patch('reec_review.inspect_sterics', return_value={'severe_clashes': 1}), \
                    patch('dohm_recovery.repair_existing', side_effect=reec_backend.CandidateRejected('Cleanup exhausted')) as repair:
                with self.assertRaises(reec_backend.CandidateRejected):
                    with backend.publication_guard(output):
                        shutil.copyfile(candidate, output)
                        backend.finalize(model, output, candidate, target='OR1B1')
            self.assertEqual(output.read_text(), 'previous model')
            self.assertEqual(Path(str(output)+'.json').read_text(), 'previous report')
            self.assertEqual(repair.call_args.kwargs['contact_report'][0]['hard'], True)
            self.assertFalse(repair.call_args.kwargs['orient'])
            self.assertEqual(repair.call_count, 1)

    def test_real_audit_detects_deliberate_nonbonded_overlap(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / 'clash.pdb'
            lines = (ROOT / 'testdata/reec-small.pdb').read_text().splitlines()
            numbers = sorted({int(line[22:26]) for line in lines if line.startswith('ATOM  ')})
            mapping = {number: index+1 for index, number in enumerate(numbers)}
            lines = [line[:22]+f'{mapping[int(line[22:26])]:4d}'+line[26:]
                     if line.startswith('ATOM  ') else line for line in lines]
            ca = [i for i, line in enumerate(lines) if line.startswith('ATOM  ') and line[12:16].strip() == 'CA']
            self.assertGreaterEqual(len(ca), 5)
            a, b = ca[0], ca[-1]
            lines[b] = lines[b][:30] + lines[a][30:54] + lines[b][54:]
            output.write_text('\n'.join(lines)+'\n')
            self.assertGreater(inspect_sterics(output)['severe_clashes'], 0)


if __name__ == '__main__':
    unittest.main()
