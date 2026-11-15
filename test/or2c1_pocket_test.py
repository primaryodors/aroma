import json
from pathlib import Path
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]


class OR2C1PocketTest(unittest.TestCase):
    def test_exact_override_is_unbiased_and_keeps_aromatic_candidates(self):
        from data import dyncenter
        odor={'full_name':'octanethiol','sdfname':'sdf/octanethiol.sdf'}
        with patch.object(dyncenter.data.odorutils,'find_odorant',return_value=odor):
            pocket=dyncenter.get_pocket('OR2C1','octanethiol')
        self.assertNotIn('!',pocket['pocket'])
        self.assertNotIn('mcoord',pocket)
        for site in ('3.32','4.60','5.46','6.55','7.41'):
            self.assertIn(site,pocket['pocket'].split())
        catalog=json.loads((ROOT/'data/binding_pocket.json').read_text())
        self.assertIn('Cu +1 C5.42 C5.43 M5.46',catalog['OR2[MTV][0-9]{1,2}']['odorophores']['S']['mcoord'])

    def test_panel_labels_preserve_measured_controls_and_unknown_nonane(self):
        catalog=json.loads((ROOT/'data/odorant.json').read_text())
        by_name={row['full_name']:row for row in catalog.values()}
        for name in ('octanethiol','nonanethiol','hexanethiol','heptanethiol'):
            self.assertTrue(any('OR2C1' in study and study['OR2C1'].get('type') in ('a','sa','ma','wa','vsa')
                                for study in by_name[name]['activity'].values()),name)
        for name in ('1-octanol','1-nonanol','nonanal'):
            self.assertEqual(by_name[name]['activity']['https://doi.org/10.1126/scisignal.2000016']['OR2C1']['type'],'na')
        self.assertFalse(any('OR2C1' in study for study in by_name.get('nonane',{}).get('activity',{}).values()))


if __name__=='__main__':
    unittest.main()
