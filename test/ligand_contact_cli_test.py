"""Ligand contact triples preserve independent distances and resolve BW selectors."""
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'hm'))
import fixfail
import modeling_backend
from data import protutils


class LigandContacts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        protutils.load_prots()

    def test_repeated_flags_have_independent_distances(self):
        args = modeling_backend.options(['OR52E8', 'example.pdb', '-lc', '208:OD1', 'O10', '2.6',
                                         '-lc', '6.59:EXTENT', 'O9', '2.5'], repair=True)
        contacts = fixfail.ligand_contact_triples('OR52E8', args.ligand_contact)
        self.assertEqual(contacts, [('O10', '208:OD1', 2.6), ('O9', '268:EXTENT', 2.5)])

    def test_optional_distances_preserve_actions_and_mixed_contacts(self):
        args = modeling_backend.options(['OR52E8', 'example.pdb', '-lc', '208', 'O10',
                                         '-lc', '6.59', 'O9', '2.5', 'dock'], repair=True)
        self.assertEqual(args.legacy, ['example.pdb', 'dock'])
        self.assertEqual(fixfail.ligand_contact_triples('OR52E8', args.ligand_contact),
                         [('O10', '208:EXTENT', None), ('O9', '268:EXTENT', 2.5)])
        for option in ('-lc', '--ligand-contact', '--ligand-contact=208:OD1'):
            values = [option, 'O10'] if '=' in option else [option, '208:OD1', 'O10']
            args = modeling_backend.options(['OR52E8', 'example.pdb', *values, '-m', '1'], repair=True)
            self.assertEqual(args.ligand_contact, [['208:OD1', 'O10', 'auto']])
            self.assertEqual(args.models, 1)

    def test_unflagged_contacts_are_rejected_and_extent_defaults(self):
        values = ['208:OD1', 'O10', '2.6', '6.59:EXTENT', 'O9', '2.5']
        for prefix in (['OR52E8'], ['OR52E8', 'example.pdb']):
            with self.assertRaisesRegex(ValueError, 'require -lc'):
                modeling_backend.options([*prefix, *values], repair=True)
        for selector in ('208', '5.42'):
            self.assertTrue(fixfail.protein_contact_selector(selector))
            self.assertEqual(fixfail.ligand_contact_triples('OR52E8', [(selector, 'O10', '2.6')]),
                             [('O10', '208:EXTENT', 2.6)])
        self.assertFalse(fixfail.protein_contact_selector('ligand.pdb'))

    def test_invalid_ligand_atom_and_distances(self):
        for ligand, distance in [('09', '2.5'), ('O9', 'nan'), ('O9', '-2'), ('O9', 'infinity')]:
            with self.subTest(ligand=ligand, distance=distance), self.assertRaises(ValueError):
                fixfail.ligand_contact_triples('OR52E8', [('208', ligand, distance)])

    def test_contacts_require_an_explicit_ligand_input(self):
        for prefix in (['OR52E8'], ['OR52E8', 'dock']):
            with self.assertRaisesRegex(ValueError, 'explicit ligand'):
                modeling_backend.options([*prefix, '-lc', '208', 'O10', '2.6'], repair=True)
        for source in ('example.pdb', '3-hydroxy-3-methyl_heptanoic_acid'):
            args = modeling_backend.options(['OR52E8', source, '-lc', '208', 'O10', '2.6'], repair=True)
            self.assertEqual(args.legacy, [source])


if __name__ == '__main__':
    unittest.main()
