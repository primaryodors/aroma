import contextlib
import io
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
sys.path.insert(0,str(ROOT/'hm'))
from data import dyncenter, odorutils, protutils
from data.preparation import generate_prepared_output
import fixfail
from reec_backend import receptor


class WorkflowAuditTest(unittest.TestCase):
    def setUp(self):
        self.previous = Path.cwd()

    def tearDown(self):
        os.chdir(self.previous)

    def test_repair_uses_canonical_filename_and_requested_state(self):
        odor = {'full_name':'ethyl phenylacetate'}
        for state, flag in [('active','noi'),('inactive','noa')]:
            self.assertEqual(fixfail.repair_dock_path('/work','OR1','OR1D2',odor,state),
                             Path(f'/work/out/OR1/OR1D2/OR1D2~ethyl_phenylacetate.{state}.dock'))
            self.assertEqual(fixfail.repair_dock_command('OR1D2',odor,state),
                             ['/bin/bash','./dock.sh','OR1D2','ethyl phenylacetate',flag])

    def test_state_suffixes_are_case_insensitive(self):
        for name in ['taar5i','TAAR5I']:
            self.assertEqual(receptor(name),('TAAR5','inactive'))
        self.assertEqual(receptor('or1a1A'),('OR1A1','active'))

    def test_full_receptor_alignment_restores_terminal_tails(self):
        records={'fixture':{'sequence':'ACDEF ','aligned':'A--CD---'}}
        with patch.object(protutils,'prots',records):
            self.assertEqual(protutils.complete_alignment('fixture'),'A--CD---EF')
            records['fixture']['aligned']='A--CE'
            with self.assertRaisesRegex(ValueError,'alignment does not match'):
                protutils.complete_alignment('fixture')
        protutils.load_prots()
        for name,record in protutils.prots.items():
            if name.startswith(('OR','TAAR')):
                self.assertEqual(protutils.complete_alignment(name).replace('-',''),record['sequence'],name)

    def test_cavity_failure_or_empty_output_preserves_cache(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)/'site.cvty'
            output.write_text('previous cavity\n')
            for failed in [True,False]:
                def producer(command, **options):
                    self.assertTrue(options['check'])
                    staged = Path(command[command.index('-o')+1])
                    if failed:
                        staged.write_text('partial cavity\n')
                        raise subprocess.CalledProcessError(1,command)
                    return subprocess.CompletedProcess(command,0)
                with patch('data.preparation.subprocess.run',side_effect=producer):
                    with self.assertRaises(subprocess.CalledProcessError if failed else RuntimeError):
                        generate_prepared_output(['scanner','-o',str(output)],output)
                self.assertEqual(output.read_text(),'previous cavity\n')

    def test_completed_cavity_replaces_cache(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)/'site.cvty'
            def producer(command, **options):
                Path(command[command.index('-o')+1]).write_text('complete cavity\n')
                return subprocess.CompletedProcess(command,0)
            with patch('data.preparation.subprocess.run',side_effect=producer):
                generate_prepared_output(['scanner','-o',str(output)],output)
            self.assertEqual(output.read_text(),'complete cavity\n')

    def test_native_cavity_fit_reports_truncated_ligand(self):
        program=ROOT/'bin/cavity_fit'
        if not program.is_file():self.skipTest('Build bin/cavity_fit for this check')
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'broken.sdf';path.write_text('truncated ligand\n')
            result=subprocess.run([str(program),str(path)],capture_output=True,text=True,timeout=10)
            self.assertNotEqual(result.returncode,0)
            self.assertIn('Invalid or incomplete ligand SDF',result.stdout+result.stderr)

    def test_ring_digits_do_not_turn_zero_moiety_matches_into_matches(self):
        data = {'OR1*':{'cen':'original','odorophores':{'C1CCCCC1':{'cen':'wrong'}}}}
        for count, expected in [(0,'original'),(1,'wrong')]:
            result = subprocess.CompletedProcess([],0,f'Pattern C1CCCCC1 occurs {count} times in molecule.\n'.encode())
            with patch.object(dyncenter,'bsrdat',data), \
                 patch.object(odorutils,'find_odorant',return_value={'sdfname':'fixture.sdf'}), \
                 patch.object(dyncenter.subprocess,'run',return_value=result) as run:
                self.assertEqual(dyncenter.get_pocket('OR1A1','fixture')['cen'],expected)
                self.assertTrue(run.call_args.kwargs['check'])

    def test_unknown_ligand_has_a_clear_lookup_error(self):
        with patch.object(odorutils,'find_odorant',return_value=False):
            with self.assertRaisesRegex(ValueError,'Unknown ligand'):
                dyncenter.get_pocket('OR1A1','absent')

    def test_invalid_structure_cannot_validate_an_arbitrary_cached_ligand(self):
        from rdkit import Chem,rdBase
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'ligand.sdf'
            with Chem.SDWriter(str(path)) as writer:
                writer.write(Chem.MolFromSmiles('CCO'))
            with rdBase.BlockLogs():
                self.assertFalse(odorutils._sdf_matches_smiles(str(path),'not a structure'))

    def test_named_atom_labels_and_branch_stereochemistry_survive_normalization(self):
        from rdkit import Chem
        molecule=odorutils.structure_molecule('CCC[C@H]([C@H](CC(=O)[O-])O)C{C13}')
        self.assertEqual(Chem.GetAtomAlias(molecule.GetAtomWithIdx(10)),'C13')
        for structure,expected in [('c1ccccc1/C=C\\(C)C=O','STEREOE'),
                                   ('c1ccccc1/C=C/(C)C=O','STEREOZ')]:
            molecule=odorutils.structure_molecule(structure)
            stereo=[str(b.GetStereo()) for b in molecule.GetBonds() if b.GetStereo()!=Chem.BondStereo.STEREONONE]
            self.assertEqual(stereo,[expected])

    def test_embedding_retries_preserve_stereochemistry(self):
        from rdkit import Chem
        original=odorutils.AllChem.EmbedMolecule
        calls=[]
        def embed(molecule,**options):
            calls.append(options)
            return -1 if len(calls)==1 else original(molecule,**options)
        with tempfile.TemporaryDirectory() as directory:
            path=str(Path(directory)/'ligand.sdf')
            structure='C[C@H](O)C(=O)O'
            with patch.object(odorutils.AllChem,'EmbedMolecule',side_effect=embed):
                odorutils.smiles_to_sdf(structure,path)
            self.assertEqual(len(calls),2)
            self.assertTrue(calls[1]['useRandomCoords'])
            self.assertTrue(odorutils._sdf_matches_smiles(path,structure))

    def test_coupled_preparation_default_arguments_and_isolated_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            os.chdir(directory)
            receptor_dir=Path('pdbs/OR1');receptor_dir.mkdir(parents=True)
            (receptor_dir/'OR1A1.active.pdb').write_text('REMARK 800 SITE BW 3.50 122\n')
            scripts=[]
            def producer(command, **options):
                self.assertTrue(options['check'])
                script=Path(command[1]);scripts.append(script)
                text=script.read_text()
                saved=re.search(r'^SAVE "([^"]+)"',text,re.M).group(1)
                Path(saved).write_text('complete model\n')
            with patch.object(protutils.subprocess,'run',side_effect=producer):
                for _ in range(2):protutils.prepare_coupled('input.pdb','output.pdb','OR1A1')
            self.assertNotEqual(scripts[0],scripts[1])
            self.assertTrue(all(not path.exists() for path in scripts))
            self.assertEqual(Path('output.pdb').read_text(),'complete model\n')

    def test_coupled_preparation_failure_preserves_previous_model(self):
        with tempfile.TemporaryDirectory() as directory:
            os.chdir(directory)
            receptor_dir=Path('pdbs/OR1');receptor_dir.mkdir(parents=True)
            (receptor_dir/'OR1A1.active.pdb').write_text('REMARK 800 SITE BW 3.50 122\n')
            output=Path('output.pdb');output.write_text('previous model\n')
            def producer(command, **options):
                text=Path(command[1]).read_text()
                saved=re.search(r'^SAVE "([^"]+)"',text,re.M).group(1)
                Path(saved).write_text('partial model\n')
                raise subprocess.CalledProcessError(1,command)
            with patch.object(protutils.subprocess,'run',side_effect=producer):
                with self.assertRaises(subprocess.CalledProcessError):
                    protutils.prepare_coupled('input.pdb',str(output),'OR1A1')
            self.assertEqual(output.read_text(),'previous model\n')

    def test_catalog_structures_and_variants_parse_and_nonadienal_has_two_E_bonds(self):
        from rdkit import Chem
        odors=json.loads((ROOT/'data/odorant.json').read_text())
        for odor in odors.values():
            with self.subTest(odor=odor['full_name']):
                mol=odorutils.structure_molecule(odor['smiles'])
                self.assertIsNotNone(mol)
                for field in ('isomers','forms'):
                    for name,value in odor.get(field,{}).items():
                        self.assertIsNotNone(odorutils.structure_molecule(value),(odor['full_name'],name))
                if odor['full_name']=='(E,E)-2,4-nonadienal':
                    self.assertEqual([str(b.GetStereo()) for b in mol.GetBonds() if b.GetStereo()!=Chem.BondStereo.STEREONONE],
                                     ['STEREOE','STEREOE'])


if __name__=='__main__':
    unittest.main()
