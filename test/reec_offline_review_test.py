"""Offline review excludes failed scores and identifies the requested atoms."""
from pathlib import Path
import sys
import tempfile
import unittest
from rdkit import Chem

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'hm'))
from reec_offline_review import metrics, paired, ligand_geometry
from reec_pose_checks import hydroxy_atom_labels, asparagine_check, copper_check


class OfflineReview(unittest.TestCase):
    def test_accepted_pose_bond_compression_is_measurable(self):
        molecule = Chem.MolFromSmiles('CC')
        conformer = Chem.Conformer(2)
        conformer.SetAtomPosition(0, (0,0,0))
        conformer.SetAtomPosition(1, (1.53,0,0))
        molecule.AddConformer(conformer)
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            sdf, pose = folder/'ligand.sdf', folder/'pose.pdb'
            writer = Chem.SDWriter(str(sdf))
            writer.write(molecule)
            writer.close()
            molecule.GetConformer().SetAtomPosition(1,(1.19,0,0))
            pose.write_text(Chem.MolToPDBBlock(molecule).replace('UNL','LIG'))
            result = ligand_geometry(pose,[sdf])
            self.assertEqual(result['heavy_bonds'],1)
            self.assertAlmostEqual(result['max_bond_change_angstrom'],.34,places=5)

    def test_zero_pose_diagnostic_score_is_not_a_success(self):
        row = metrics('Total: -58.5\nDisqualified because: rejected\n0 pose(s) found.\n', 0, 1)
        self.assertEqual(row['diagnostic_first_energy'], -58.5)
        self.assertIsNone(row['first_energy'])
        self.assertFalse(row['operational_pass'])
        self.assertFalse(row['pose_count_agrees'])
        self.assertIsNone(metrics('Total: -5\n4 pose(s) found.\n', 'timeout', 4)['first_energy'])
        self.assertIsNone(metrics('Total: nan\n4 pose(s) found.\n', 0, 4)['first_energy'])

    def test_few_poses_allow_exploratory_scores_but_fail_the_frozen_rule(self):
        row = metrics('Total: -5\n3 pose(s) found.\n', 0, 3)
        self.assertEqual(row['first_energy'], -5)
        self.assertFalse(row['operational_pass'])
        self.assertTrue(metrics('Total: -5\n4 pose(s) found.\n', 0, 4)['operational_pass'])

    def test_pairing_retains_shared_input_and_censors_failed_scores(self):
        row = {'receptor':'OR6P1', 'ligand':'phenylacetaldehyde', 'seed':1, 'job_id':'same',
               'exit_code':'timeout', 'poses':0, 'operational_pass':False, 'first_energy':None}
        pair = paired([dict(row,condition=0),dict(row,condition=64)])[0]
        self.assertTrue(pair['same_native_job'])
        self.assertIsNone(pair['delta_on_minus_off'])
        self.assertFalse(pair['both_operational_pass'])

    def test_hydroxy_oxygen_is_distinct_from_carboxyl_oxygens(self):
        molecule = Chem.AddHs(Chem.MolFromSmiles('O=C(O)CC(O)(C)CCC'))
        self.assertEqual(hydroxy_atom_labels(molecule), ['O6'])
        protein = {(208,'OD1'):{'point':(0,0,0),'element':'O'}}
        ligand = {'O1':{'point':(0,0,2.8),'element':'O'},
                  'O6':{'point':(0,0,5),'element':'O'}}
        self.assertFalse(asparagine_check(protein, ligand, ['O6'])['contact_pass'])
        ligand['O6']['point']=(0,0,2.8)
        ligand['H7']={'point':(0,0,1.8),'element':'H'}
        result=asparagine_check(protein,ligand,['O6'])
        self.assertTrue(result['contact_pass'])
        self.assertTrue(result['hydrogen_bond_pass'])

    def test_copper_ligand_contact_does_not_substitute_for_two_cysteines(self):
        protein={(199,'SG'):{'point':(-2.2,0,0)}, (200,'SG'):{'point':(2.3,0,0)}}
        copper=[{'point':(0,0,0)}]
        ligand={'S1':{'point':(0,2.5,0),'element':'S'}}
        result=copper_check(protein,ligand,copper)
        self.assertTrue(result['coordination_pass'])
        self.assertTrue(result['ligand_copper_contact'])
        protein[200,'SG']['point']=(3.4,0,0)
        result=copper_check(protein,ligand,copper)
        self.assertFalse(result['coordination_pass'])
        self.assertTrue(result['ligand_copper_contact'])


if __name__ == '__main__':
    unittest.main()
