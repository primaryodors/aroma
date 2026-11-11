import json
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from data import odorutils as ou, protutils as pu


class Taar5ActivityTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ou.load_odors()
        pu.load_prots()

    def test_taar5_top_agonist_is_trimethylamine(self):
        ligand = ou.find_odorant(pu.prots['TAAR5']['best_agonist'])
        self.assertEqual(ligand['full_name'], 'trimethylamine')
        self.assertIn(ligand['oid'], ou.empirical_pairs('TAAR5', agonists_only=True))

    def test_thyronamine_is_sourced_as_an_inverse_agonist(self):
        ligand = ou.find_odorant('3-iodo-thyronamine')
        activity = ligand['activity']
        self.assertNotIn('TAAR5', activity['https://doi.org/10.1111/cbdd.12903'])
        record = activity['https://doi.org/10.1371/journal.pone.0117774']['TAAR5']
        self.assertEqual(record['type'], 'ia')
        self.assertNotIn('adjusted_curve_top', record)
        self.assertNotIn(ligand['oid'], ou.empirical_pairs('TAAR5', agonists_only=True))

    def test_inverse_agonist_without_numeric_curve_is_excluded(self):
        records = {'inverse': {'activity': {'source': {'TAAR5': {'type': 'ia'}}}},
                   'positive': {'activity': {'source': {'TAAR5': {'type': 'a'}}}}}
        with patch.dict(ou.odors, records, clear=True):
            self.assertEqual(set(ou.empirical_pairs('TAAR5', agonists_only=True)), {'positive'})

    def test_activity_source_has_a_registered_citation(self):
        source = 'https://doi.org/10.1371/journal.pone.0117774'
        references = json.loads((ROOT / 'data/reference.json').read_text())
        ligand = ou.find_odorant('3-iodo-thyronamine')
        self.assertIn(source, ligand['activity'])
        self.assertEqual(references[source]['short_name'], 'Dinter et al (2015)')
        self.assertIn('10.1371/journal.pone.0117774', references[source]['citation'])


if __name__ == '__main__':
    unittest.main()
