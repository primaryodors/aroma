"""Metal fitting is shared by config and PHEW, while dohm placement is opt-in."""
import json
import os
import re
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'hm'))
from data import protutils as pu
from data.metal_sites import configured_metal_sites, inspect_metal_sites
from modeling_backend import options, forwarded_options


class MetalCoordination(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pu.load_prots()

    def test_dohm_defaults_to_ligand_dependent_placement(self):
        args=options(['OR2T11'])
        self.assertFalse(args.coordinate_metals)
        self.assertNotIn('--coordinate-metals',forwarded_options(args))
        enabled=options(['OR2T11','--coordinate-metals'])
        self.assertTrue(enabled.coordinate_metals)
        self.assertIn('--coordinate-metals',forwarded_options(enabled))

    def test_cavity_frontend_rejects_an_absent_donor_without_looping(self):
        result=subprocess.run([str(ROOT/'bin/cavity_search'),'-p',str(ROOT/'testdata/metal-site-OR2T11.pdb'),
                               '--metal','Cu','+1','199','200','999'],cwd=ROOT,capture_output=True,text=True,timeout=10)
        self.assertNotEqual(result.returncode,0)
        self.assertIn('absent or has no donor',result.stdout+result.stderr)

    def test_policy_uses_the_configured_donor_motif(self):
        sites=configured_metal_sites('OR2T11')
        self.assertEqual(len(sites),1)
        self.assertEqual(sites[0]['residues'],[199,200,203])
        self.assertEqual(sites[0]['charge'],1)
        self.assertEqual(configured_metal_sites('OR2AG1'),[])
        with tempfile.TemporaryDirectory() as temporary:
            policy=Path(temporary)/'policy.json'
            policy.write_text(json.dumps({'OR2T11':{'mcoord':['Cu +1 M5.42 C5.43 M5.46']}}))
            self.assertEqual(configured_metal_sites('OR2T11',policy),[])
            policy.write_text(json.dumps({'OR2T11':{'mcoord':'Cu +1 C5.42 C5.43 M5.46'}}))
            self.assertEqual(configured_metal_sites('OR2T11',policy)[0]['residues'],[199,200,203])

    def test_docking_config_calls_shared_fitter_and_writes_one_ion(self):
        from rdkit import Chem
        from rdkit.Chem import AllChem
        with tempfile.TemporaryDirectory() as temporary:
            folder=Path(temporary)
            (folder/'data').symlink_to(ROOT/'data')
            (folder/'sdf').symlink_to(ROOT/'sdf')
            (folder/'tmp').mkdir()
            molecule=Chem.AddHs(Chem.MolFromSmiles('CS'))
            AllChem.EmbedMolecule(molecule,randomSeed=1)
            writer=Chem.SDWriter(str(folder/'ligand.sdf'))
            writer.write(molecule);writer.close()
            config=folder/'dock.config'
            config.write_text(f'PROT {ROOT}/testdata/metal-site-OR2T11.pdb\nLIG {folder}/ligand.sdf\n'
                              'CEN ABS 4.3 11.4 10.7\nSIZE 5\nMCOORD Cu +1 C5.42 C5.43 M5.46\n'
                              'SEARCH TS\nPOSE 1\nITERS 1\nTHREADS 1\nRETRY 0\nNOFAIL\nMOVIE\n'
                              f'OUT {folder}/result.dock\nOUTPDB 1 {folder}/pose%o.pdb\n')
            result=subprocess.run([str(ROOT/'bin/aromadock'),str(config),'--ns','-t','1'],
                                  cwd=folder,capture_output=True,text=True,timeout=60)
            self.assertEqual(result.returncode,0,result.stdout[-2000:]+result.stderr)
            pose=folder/'pose1.pdb'
            self.assertTrue(pose.is_file())
            inspect_metal_sites(pose,configured_metal_sites('OR2T11'))
            self.assertEqual(sum(line.startswith(('ATOM  ','HETATM')) and line[76:78].strip()=='Cu'
                                 for line in pose.read_text().splitlines()),1)
            movie=(folder/'tmp/_iters.dock').read_text()
            frames=re.split(r'^Pose: [0-9]+\n',movie,flags=re.M)[1:]
            self.assertGreaterEqual(len(frames),2)
            for frame in frames:
                records=frame.split('\nEND\n')[0].splitlines()
                self.assertEqual(sum(line.startswith(('ATOM  ','HETATM')) and line[76:78].strip()=='Cu'
                                     for line in records),1)
                self.assertTrue(any(line.startswith('HETATM') and line[17:20]=='LIG' for line in records))

    def test_load_plus_keeps_ligand_separate_from_protein_ion(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder=Path(temporary)
            output=folder/'metal.pdb'
            subprocess.run([str(ROOT/'test/metal_fit_test'),str(ROOT/'testdata/metal-site-OR2T11.pdb'),str(output)],
                           cwd=ROOT,check=True,capture_output=True,timeout=30)
            ligand='HETATM 8000  C1  LIG A 999      30.000  30.000  30.000  1.00  0.00           C  \n'
            output.write_text(output.read_text().replace('TER\n',ligand+'TER\n'))
            first_model=output.read_text().replace('END\n','')
            second_model=first_model.replace('30.000','60.000')
            output.write_text('MODEL 1\n'+first_model+'ENDMDL\nMODEL 2\n'+second_model+'ENDMDL\nEND\n')
            for suffix,count in ((' +',1),('',0)):
                saved=folder/f'reload-{count}.pdb'
                script=folder/'reload.phew'
                script.write_text(f'LOAD "{output}" A A{suffix}\nSTRAND A\nSAVE "{saved}"\n')
                result=subprocess.run([str(ROOT/'bin/phew'),str(script)],cwd=ROOT,capture_output=True,text=True,timeout=30)
                self.assertEqual(result.returncode,0,result.stdout+result.stderr)
                lines=saved.read_text().splitlines()
                self.assertEqual(sum(line.startswith(('ATOM  ','HETATM')) and line[76:78].strip()=='Cu' for line in lines),1)
                self.assertEqual(sum(line.startswith('HETATM') and line[17:20]=='LIG' for line in lines),count)

    def test_best_binding_retains_a_scored_metal_contact(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder=Path(temporary)
            (folder/'data').symlink_to(ROOT/'data')
            (folder/'sdf').symlink_to(ROOT/'sdf')
            (folder/'tmp').mkdir()
            config=folder/'dock.config'
            clock=folder/'clock.so'
            subprocess.run(['gcc','-shared','-fPIC','-O2',str(ROOT/'test/reec_baseline_clock.c'),
                            '-ldl','-o',str(clock)],check=True,capture_output=True)
            environment=dict(os.environ,LD_PRELOAD=str(clock),AROMA_BASELINE_TIME='1770000000',OMP_NUM_THREADS='1')
            config.write_text(f'PROT {ROOT}/testdata/metal-site-OR2T11.pdb\n'
                              f'LIG {ROOT}/sdf/diallyl_disulfide.sdf\n'
                              'CEN ABS 4.3 11.4 10.7\nSIZE 5\nMCOORD Cu +1 C5.42 C5.43 M5.46\n'
                              'SEARCH BB\nPOSE 2\nITERS 3\nTHREADS 1\nRETRY 0\nNOFAIL\n'
                              f'OUT {folder}/result.dock\nOUTPDB 1 {folder}/pose%o.pdb\n')
            result=subprocess.run([str(ROOT/'bin/aromadock'),str(config),'--ns','-t','1'],
                                  cwd=folder,capture_output=True,text=True,timeout=90,env=environment)
            self.assertEqual(result.returncode,0,result.stdout[-2000:]+result.stderr)
            report=(folder/'result.dock').read_text()
            energies=re.findall(r'^Total metal coordination: ([-+0-9.eE]+)',report,re.M)
            self.assertTrue(energies,report[:2000])
            self.assertTrue(all(float(energy)<0 for energy in energies))
            self.assertNotRegex(result.stdout,r'Staying (?:C|H)[0-9]+ near (?:199|200|203):Cu1')
            import math
            for pose in folder.glob('pose*.pdb'):
                inspect_metal_sites(pose,configured_metal_sites('OR2T11'))
                records=[line for line in pose.read_text().splitlines() if line.startswith(('ATOM  ','HETATM'))]
                ions=[line for line in records if line[76:78].strip()=='Cu']
                self.assertEqual(len(ions),1)
                coordinates=lambda line: tuple(float(line[i:i+8]) for i in (30,38,46))
                distance=min(math.dist(coordinates(line),coordinates(ions[0])) for line in records
                             if line[17:20]=='LIG' and line[76:78].strip()=='S')
                self.assertLess(distance,2.9)

    def test_native_fit_preserves_geometry_and_roundtrips(self):
        with tempfile.TemporaryDirectory() as temporary:
            output=Path(temporary)/'roundtrip.pdb'
            result=subprocess.run([str(ROOT/'test/metal_fit_test'),str(ROOT/'testdata/metal-site-OR2T11.pdb'),str(output)],
                                  cwd=ROOT,capture_output=True,text=True,timeout=30)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            sites=inspect_metal_sites(output,configured_metal_sites('OR2T11'))
            self.assertEqual(len(sites),1)
            self.assertEqual(sum(line.startswith(('ATOM  ','HETATM')) and line[76:78].strip()=='Cu'
                                 for line in output.read_text().splitlines()),1)

    def test_native_scoring_counts_each_ion_once(self):
        result=subprocess.run([str(ROOT/'test/metal_scoring_test'),str(ROOT/'testdata/metal-site-OR2T11.pdb')],
                              cwd=ROOT,capture_output=True,text=True,timeout=30)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)

    def test_saves_preserve_two_ions_and_collapse_legacy_duplicate_records(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder=Path(temporary)
            source=folder/'source.pdb'
            subprocess.run([str(ROOT/'test/metal_fit_test'),str(ROOT/'testdata/metal-site-OR2T11.pdb'),str(source)],
                           cwd=ROOT,check=True,capture_output=True,timeout=30)
            copper=next(line for line in source.read_text().splitlines() if line.startswith('HETATM'))
            duplicate='ATOM  '+copper[6:]
            second=list(copper)
            second[6:11]=f'{9002:5d}';second[12:16]=' Cu2';second[21]='A';second[22:26]=f'{299:4d}'
            second[30:38]=f'{float(copper[30:38])+20:8.3f}'
            source.write_text(source.read_text().replace('TER\n',duplicate+'\n'+''.join(second)+'\nTER\n'))
            output=folder/'output.pdb'
            script=folder/'save.phew'
            script.write_text(f'LOAD "{source}" A A\nSTRAND A\nSAVE "{output}"\n')
            result=subprocess.run([str(ROOT/'bin/phew'),str(script)],cwd=ROOT,capture_output=True,text=True,timeout=30)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            ions=[line for line in output.read_text().splitlines() if line.startswith('HETATM')]
            self.assertEqual(len(ions),2)
            self.assertEqual({line[12:16].strip() for line in ions},{'Cu1','Cu2'})
            script.write_text(f'LOAD "{source}" A A\nLOAD "{ROOT}/testdata/metal-site-OR2T11.pdb" A A\nSTRAND A\nSAVE "{output}"\n')
            result=subprocess.run([str(ROOT/'bin/phew'),str(script)],cwd=ROOT,capture_output=True,text=True,timeout=30)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            self.assertFalse(any(line.startswith('HETATM') for line in output.read_text().splitlines()))

    def test_phew_accepts_config_style_residues_and_retains_one_ion(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder=Path(temporary)
            output=folder/'coordinated.pdb'
            script=folder/'metal.phew'
            script.write_text(f'LOAD "{ROOT}/testdata/metal-site-OR2T11.pdb" A A\nSTRAND A\n'
                              'MCOORD Cu 1 C5.42 C5.43 M5.46\n'
                              'MCOORD Cu 1 C5.42 C5.43 M5.46\n'
                              f'SAVE "{output}"\n')
            result=subprocess.run([str(ROOT/'bin/phew'),str(script)],cwd=ROOT,capture_output=True,text=True,timeout=30)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            self.assertEqual(len(inspect_metal_sites(output,configured_metal_sites('OR2T11'))),1)
            self.assertEqual(sum(line.startswith(('ATOM  ','HETATM')) and line[76:78].strip()=='Cu' for line in output.read_text().splitlines()),1)
            script.write_text(f'LOAD "{output}" A A\nSTRAND A\nMCOORD Cu 1 199 199 203\n')
            result=subprocess.run([str(ROOT/'bin/phew'),str(script)],cwd=ROOT,capture_output=True,text=True,timeout=30)
            self.assertNotEqual(result.returncode,0)
            self.assertIn('absent or repeated',result.stdout+result.stderr)


if __name__=='__main__':
    unittest.main()
