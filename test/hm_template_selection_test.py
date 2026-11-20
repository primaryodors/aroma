"""Template selection follows the receptor phylogram rather than family names."""
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
sys.path.insert(0,str(ROOT/'hm'))
from data import protutils as pu


class TemplateSelection(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pu.load_prots()

    def test_specific_clades_and_cross_family_members(self):
        cases = {
            'OR1A1':{'8uxy','9w45'}, 'OR7A2':{'8uxy','9w45'},
            'OR2A1':{'8uy0'}, 'OR13J1':{'8uy0'},
            'OR4A4':{'8uyq'}, 'OR8B2':{'9wpm'},
            'OR2AP1':{'9ldv','9ldw','9ldx','9ldz'},
            'OR9A2':{'9ldv','9ldw','9ldx','9ldz'},
            'OR9A4':{'9ldv','9ldw','9ldx','9ldz'},
            'OR5V1':{'9lkb','9wpm'},
        }
        for receptor,expected in cases.items():
            with self.subTest(receptor=receptor):
                self.assertEqual(set(pu.templates_for_hm(receptor)),expected)

    def test_all_schiff_members_keep_the_three_specific_templates(self):
        for receptor in ('OR6A2','OR6B1','OR6B2','OR6B3','OR6P1','OR6Y1'):
            with self.subTest(receptor=receptor):
                self.assertEqual(pu.templates_for_hm(receptor),['9le0','9le1','9le2'])

    def test_distant_class_ii_members_use_the_broad_pool(self):
        expected = {'8uxy','8uy0','8uyq','9wpm','9ldv','9ldw','9ldx','9ldz'}
        for receptor in ('OR2AT4','OR3A1','OR12D1','OR13A1','OR13G1',
                         'OR10AD1','OR10A2','OR11A1','OR14A2','OR6K2','OR8S1'):
            with self.subTest(receptor=receptor):
                self.assertEqual(set(pu.templates_for_hm(receptor)),expected)

    def test_binary_branch_order_does_not_change_selection(self):
        expected = {name:pu.templates_for_hm(name) for name in pu.prots}
        swapped = {name:dict(record,btree=record.get('btree','').translate(str.maketrans('01','10')))
                   for name,record in pu.prots.items()}
        with patch.object(pu,'prots',swapped):
            for name,result in expected.items():
                self.assertEqual(pu.templates_for_hm(name),result,name)

    def test_ensemble_does_not_broaden_or6a2_selection(self):
        from dohm import experimental_ensemble_alignment
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            (work/'tpl').symlink_to(ROOT/'hm/tpl')
            shutil.copyfile(ROOT/'hm/experimental.ali',work/'experimental.ali')
            _,templates,_ = experimental_ensemble_alignment('OR6A2',work)
            self.assertEqual([value.split('=',1)[0] for value in templates],['9le0','9le1','9le2'])


if __name__=='__main__':
    unittest.main()
