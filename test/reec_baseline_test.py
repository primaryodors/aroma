"""Checks for baseline superposition and acceptance accounting."""
import sys
import unittest
import tempfile
from unittest.mock import patch
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'hm'))
from reec_baseline import pocket_distance, summarize, prepare_repair, paired_summary, target_row


class BaselineMetrics(unittest.TestCase):
    def setUp(self):
        self.first = {(1, 'CA'): [0., 0., 0.], (2, 'CA'): [2., 0., 0.],
                      (3, 'CA'): [0., 3., 0.], (1, 'CB'): [0., 0., 1.]}

    def test_rigid_motion_is_removed(self):
        moved = {key: [-y+5, x-2, z+9] for key, (x, y, z) in self.first.items()}
        result = pocket_distance(self.first, moved, {1})
        self.assertAlmostEqual(result['rmsd_angstrom'], 0., places=10)
        self.assertEqual(result['common_pocket_atoms'], 2)

    def test_sidechain_motion_survives_superposition(self):
        moved = dict(self.first)
        moved[1, 'CB'] = [0., 0., 3.]
        result = pocket_distance(self.first, moved, {1})
        self.assertAlmostEqual(result['rmsd_angstrom'], 2**.5)
        backbone = pocket_distance(self.first, moved, {1}, {'CA'})
        self.assertAlmostEqual(backbone['rmsd_angstrom'], 0., places=10)

    def test_insufficient_correspondence_is_explicit(self):
        self.assertIsNone(pocket_distance(self.first, {}, {1}))
        self.assertIsNone(pocket_distance(self.first, self.first, {99}))

    def test_rejected_runs_remain_in_acceptance_denominator(self):
        runs = [{'anneal_steps': 0, 'exit_code': 2, 'wall_seconds': 2.,
                 'candidates': [{'valid': False}, {'valid': False}]},
                {'anneal_steps': 0, 'exit_code': 0, 'wall_seconds': 4.,
                 'candidates': [{'valid': True}, {'valid': False}]}]
        result = summarize(runs, [], {1})['groups'][0]
        self.assertEqual(result['runs'], 2)
        self.assertEqual(result['accepted_runs'], 1)
        self.assertEqual(result['candidates'], 4)
        self.assertEqual(result['accepted_candidates'], 1)
        self.assertEqual(result['median_wall_seconds'], 3.)

    def test_repair_snapshots_are_isolated_from_original_inputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            source = base / 'retained'
            source.mkdir()
            alignment = '>P1;OR1B1_tpl\nstructure:x:FIRST:A:LAST:A::::\nA*\n'
            alignment += '>P1;OR1B1\nsequence:x::::::::\nA*\n'
            (source / 'protein.ali').write_text(alignment)
            (source / 'receptor.contacts').write_text('disulfide 98 184\n')
            (source / 'environment.pdb').write_text('END\n')
            template, reference = base / 'template.pdb', base / 'reference.pdb'
            template.write_text('template snapshot\n')
            reference.write_text('reference snapshot\n')
            with patch('reec_baseline.pu.resno_from_bw', return_value=None):
                output, flags, inputs, pocket = prepare_repair(base / 'case', source, template, reference)
            self.assertEqual(output.read_text(), alignment)
            output.write_text('changed local copy\n')
            self.assertEqual((source / 'protein.ali').read_text(), alignment)
            reference.write_text('changed external reference\n')
            self.assertEqual((base / 'case/reference.pdb').read_text(), 'reference snapshot\n')
            self.assertIn('--environment-hydrogens', flags)
            self.assertEqual(len(inputs), 5)

    def test_paired_outcomes_distinguish_rescue_loss_and_changed_weights(self):
        def candidate(index, valid, energy, released=()):
            return {'candidate_index': index, 'start_signature': str(index), 'valid': valid,
                    'independent_review': {'valid': valid}, 'final_energy': energy,
                    'released_template_residues': list(released)}
        runs = [{'seed': 0, 'anneal_steps': 0, 'wall_seconds': 2.,
                 'candidates': [candidate(0, False, 10), candidate(1, True, 100), candidate(2, True, 200)]},
                {'seed': 0, 'anneal_steps': 64, 'wall_seconds': 3.,
                 'candidates': [candidate(0, True, 100), candidate(1, False, 10), candidate(2, True, 50, (7,))]}]
        result = paired_summary(runs)
        self.assertEqual(result['off_accepted'], result['on_accepted'])
        self.assertEqual(result['rescued'], 1)
        self.assertEqual(result['lost'], 1)
        self.assertEqual(result['median_runtime_ratio'], 1.5)
        self.assertFalse(result['candidate_outcomes'][2]['score_comparable'])
        self.assertIsNone(result['best_scores'][0]['score_improvement_percent'])
        runs[1]['candidates'][2]['start_signature'] = 'different'
        with self.assertRaises(ValueError):
            paired_summary(runs)

    def test_explicit_terminal_extension_preserves_core_and_rejects_internal_changes(self):
        records = {'model': {'aligned': 'A-C-D----', 'sequence': 'ACDEFG'}}
        with patch('reec_baseline.pu.prots', records):
            with self.assertRaises(ValueError):
                target_row('model')
            extended = target_row('model', True)
            self.assertTrue(extended.startswith('A-C-D'))
            self.assertEqual(extended.replace('-', ''), 'ACDEFG')
            records['model']['sequence'] = 'ACFDEFG'
            with self.assertRaises(ValueError):
                target_row('model', True)


if __name__ == '__main__':
    unittest.main()
