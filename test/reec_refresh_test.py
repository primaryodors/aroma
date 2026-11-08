"""Frozen docking-job deduplication preserves conditions and replication."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'hm'))
from reec_refresh_benchmark import prepare_jobs, digest


class RefreshJobs(unittest.TestCase):
    def test_identical_conditions_share_work_without_losing_seed_aliases(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / 'bin').mkdir()
            (root / 'results').mkdir()
            (root / 'bin/aromadock').write_bytes(b'frozen engine')
            (root / 'bin/reec').write_bytes(b'frozen model engine')
            model = root / 'model.pdb'
            model.write_text('model')
            ligand = root / 'ligand.sdf'
            ligand.write_text('ligand')
            config = root / 'base.config'
            config.write_text(f'PROT {model}\nLIG {ligand}\nPOSE 10\nITERS 50\n')
            case = {'receptor': 'OR2AT4', 'ligand': 'sandalore', 'conditions': [0, 64]}
            prepared = {'case': case, 'config': str(config), 'model_sha256': digest(model),
                        'sdfs': [{'path': str(ligand)}]}
            for name, seeds, conditions in [('top-panel', [101, 202, 303], [0, 64]),
                                             ('full-panel', [101], [0])]:
                panel = {'cases': [dict(case, conditions=conditions)],
                         'settings': {'clock_seeds': seeds, 'poses': 10, 'iterations': 50, 'threads': 1}}
                manifest = {'pending': [], 'cases': [dict(prepared, condition=c) for c in conditions]}
                (root / 'results' / f'{name}.json').write_text(json.dumps(panel))
                (root / 'results' / f'{name}.prepared.json').write_text(json.dumps(manifest))
            self.assertEqual(prepare_jobs(root), 3)
            jobs = json.loads((root / 'jobs.json').read_text())['jobs']
            self.assertEqual({job['clock_seed'] for job in jobs}, {101, 202, 303})
            self.assertEqual(sum(len(job['aliases']) for job in jobs), 7)
            for job in jobs:
                text = (root / 'jobs' / job['id'] / 'job.config').read_text()
                self.assertEqual(text.count('\nPOSE '), 1)
                self.assertIn('THREADS 1\nRETRY 0\n', text)
                self.assertIn(f'OUT {root}/jobs/{job["id"]}/result.dock', text)
            path = root / 'results/top-panel.prepared.json'
            incomplete = json.loads(path.read_text())
            incomplete['pending'] = [case]
            path.write_text(json.dumps(incomplete))
            with self.assertRaisesRegex(RuntimeError, 'unresolved models'):
                prepare_jobs(root)


if __name__ == '__main__':
    unittest.main()
