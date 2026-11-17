"""Schiff-adduct repair regressions from isolated OR6A2 and OR6P1 fragments."""
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'hm'))
import schiff
from modeling_backend import preserve_environment_labels

CODES = dict(zip('ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR'.split(),
                 'ACDEFGHIKLMNPQRSTVWY'))


class SchiffRepairTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.work = Path(self.temporary.name)
        self.source = self.work/'input.pdb'
        self.source.write_text((ROOT/'testdata/schiff_octanal.pdb').read_text())
        self.sdf = ROOT/'sdf/octanal.sdf'

    def prepare(self):
        adducts = schiff.discover(self.source,self.sdf)
        schiff.annotate(self.source,adducts,self.sdf)
        return adducts

    def test_legacy_inputs_identify_the_iminium_state(self):
        for odor in ('octanal','phenylacetaldehyde'):
            with self.subTest(odor=odor):
                self.source.write_text((ROOT/f'testdata/schiff_{odor}.pdb').read_text())
                adducts = schiff.discover(self.source,ROOT/f'sdf/{odor}.sdf')
                self.assertEqual(adducts,[schiff.Link(2,'NZ',999,'C8','C7','O9','HZ3')])

    def test_missing_reference_and_ambiguous_assignment_are_rejected(self):
        with self.assertRaisesRegex(ValueError,'SDF reference'):
            schiff.discover(self.source)
        text = self.source.read_text()
        line = next(l for l in text.splitlines() if l.startswith('ATOM  ') and l[12:16].strip()=='NZ')
        self.source.write_text(text.replace('END\n','')+line[:6]+'90000'+line[11:22]+'  50'+line[26:]+'\nEND\n')
        with self.assertRaisesRegex(ValueError,'Ambiguous Schiff assignment'):
            schiff.discover(self.source,self.sdf)

    def test_intact_carbonyl_is_not_an_adduct(self):
        text = self.source.read_text()
        line = next(l for l in text.splitlines() if l.startswith('HETATM') and l[12:16].strip()=='C8')
        oxygen = line[:6]+'90000'+'   O9'+line[16:76]+' O'
        self.source.write_text(text.replace('END\n','')+oxygen+'\nEND\n')
        self.assertEqual(schiff.discover(self.source,self.sdf),[])

    def test_ordinary_protein_is_unchanged(self):
        self.source.write_text((ROOT/'testdata/reec-small.pdb').read_text())
        original = self.source.read_bytes()
        self.assertEqual(schiff.discover(self.source,self.sdf),[])
        schiff.annotate(self.source,[])
        self.assertEqual(original,self.source.read_bytes())

    def test_neutral_imine_retains_no_nitrogen_hydrogens(self):
        rows = []
        for line in self.source.read_text().splitlines():
            if line.startswith('ATOM  ') and line[12:16].strip()=='HZ3':
                continue
            if line.startswith('ATOM  ') and line[12:16].strip()=='NZ':
                line = line.ljust(80)[:78]+'  '
            rows.append(line)
        self.source.write_text('\n'.join(rows)+'\n')
        adducts = self.prepare()
        self.assertEqual(adducts[0].hydrogen,'-')
        self.assertEqual(schiff.failures(self.source,require_hydrogens=True),[])

    def test_connectivity_survives_serial_changes_and_hydrogen_cleanup(self):
        self.prepare()
        candidate = self.source.read_text()
        target = self.work/'postprocessed.pdb'
        rows = []
        for line in candidate.splitlines():
            if line.startswith(('ATOM  ','HETATM')):
                line = line[:6]+f'{int(line[6:11])+100:5d}'+line[11:]
                rows.append(line)
        nz = next(l for l in rows if l[12:16].strip()=='NZ')
        rows.append(nz[:6]+'95000'+nz[11:12]+' HZ1'+nz[16:76]+' H')
        target.write_text('\n'.join(rows)+'\nEND\n')
        preserve_environment_labels(self.source,target)
        self.assertEqual(schiff.failures(target),[])
        self.assertEqual(schiff.links(target.read_text()),schiff.links(candidate))
        inventory = schiff.atoms(target.read_text())
        self.assertIn((2,'HZ3'),inventory)
        self.assertNotIn((2,'HZ1'),inventory)
        self.assertEqual(inventory[2,'NZ'][0][78:80],'1+')
        self.assertAlmostEqual(math.dist(inventory[2,'NZ'][1],inventory[2,'HZ3'][1]),1.01,places=2)
        self.assertTrue(any(l.startswith('CONECT') for l in target.read_text().splitlines()))

    def test_explicit_connectivity_can_recover_a_stretched_link(self):
        self.prepare()
        rows = [l for l in self.source.read_text().splitlines() if not l.startswith(schiff.PREFIX)]
        for i,line in enumerate(rows):
            if line.startswith('ATOM  ') and line[12:16].strip() in ('NZ','HZ3'):
                rows[i] = line[:30]+f'{float(line[30:38])+1:8.3f}'+line[38:]
        self.source.write_text('\n'.join(rows)+'\n')
        self.assertEqual(len(schiff.discover(self.source,self.sdf)),1)

    def test_validation_catches_broken_bond_and_restored_oxygen(self):
        self.prepare()
        rows = self.source.read_text().splitlines()
        for i,line in enumerate(rows):
            if line.startswith('HETATM') and line[12:16].strip()=='C8':
                rows[i] = line[:30]+f'{float(line[30:38])+3:8.3f}'+line[38:]
                rows.append(line[:6]+'90000'+'   O9'+line[16:76]+' O')
                break
        self.source.write_text('\n'.join(rows)+'\n')
        failures = schiff.failures(self.source)
        self.assertTrue(any('distance' in s for s in failures))
        self.assertIn('Departed Schiff oxygen reappeared',failures)

    def test_fresh_dock_input_omits_adduct_but_keeps_other_environment(self):
        adducts = self.prepare()
        text = self.source.read_text()
        metal = 'HETATM90000  Cu1 MTL A 999      40.000  40.000  40.000  1.00  0.00          Cu'
        self.source.write_text(text.replace('END\n','')+metal+'\nEND\n')
        schiff.receptor_only(self.source,adducts)
        text = self.source.read_text()
        self.assertFalse(schiff.links(text))
        self.assertNotIn(' LIG ',text)
        self.assertIn(metal,text)
        self.assertIn(' NZ ',text)

    def test_native_repair_and_serialized_audit_on_both_fragments(self):
        for odor in ('octanal','phenylacetaldehyde'):
            with self.subTest(odor=odor):
                self.source.write_text((ROOT/f'testdata/schiff_{odor}.pdb').read_text())
                sdf = ROOT/f'sdf/{odor}.sdf'
                schiff.annotate(self.source,schiff.discover(self.source,sdf),sdf)
                sequence = ''.join(CODES[l[17:20]] for l in self.source.read_text().splitlines()
                                   if l.startswith('ATOM  ') and l[12:16].strip()=='CA')
                alignment = self.work/'model.ali'
                alignment.write_text(f'>P1;model\nsequence:model::::::::\n{sequence}*\n')
                output = self.work/'output.pdb'
                result = subprocess.run([str(ROOT/'bin/reec'),'-a',str(alignment),'-t',str(self.source),
                                         '-e',str(self.source),'-H','-o',str(output),'-i','3000'],
                                        cwd=ROOT,capture_output=True,text=True)
                self.assertEqual(result.returncode,0,result.stderr+result.stdout[-1000:])
                self.assertEqual(schiff.failures(output),[])
                report = json.loads(Path(str(output)+'.json').read_text())
                self.assertTrue(report['valid'])
                self.assertEqual(report['severe_clashes'],0)
                self.assertEqual(report['mandatory_violations'],0)
                self.assertTrue(any(c['kind']=='schiff' for c in report['contacts']))
                before,after = schiff.atoms(self.source.read_text()),schiff.atoms(output.read_text())
                for key,(_,point) in before.items():
                    if key[0]==999:
                        self.assertEqual(point,after[key][1])
                adducts = schiff.links(output.read_text())
                schiff.annotate(output,adducts)
                audit = subprocess.run([str(ROOT/'bin/reec_steric_audit'),str(output)],
                                       cwd=ROOT,capture_output=True,text=True,check=True)
                audit = json.loads(audit.stdout)
                self.assertEqual(audit['severe_clashes'],0)
                self.assertTrue(audit['covalent_geometry']['valid'])
                self.assertEqual(audit['fixed_heavy_atoms'],len([name for r,name in before if r==999 and not name.startswith('H')]))
                hydrated = self.work/'hydrated.pdb'
                script = self.work/'hydrate.phew'
                script.write_text(f'LOAD "{output}" A A +\nSTRAND A\nHYDRO\nSAVE "{hydrated}"\n')
                subprocess.run([str(ROOT/'bin/phew'),str(script)],cwd=ROOT,
                               capture_output=True,text=True,check=True)
                preserve_environment_labels(output,hydrated)
                self.assertEqual(schiff.failures(hydrated,require_hydrogens=True),[])
                inventory = schiff.atoms(hydrated.read_text())
                self.assertEqual({n for r,n in inventory if r==2 and n.startswith('HZ')},{'HZ3'})
                schiff.receptor_only(hydrated,adducts)
                receptor = self.work/'receptor.pdb'
                script.write_text(f'LOAD "{hydrated}" A A +\nSTRAND A\nHYDRO\nSAVE "{receptor}"\n')
                subprocess.run([str(ROOT/'bin/phew'),str(script)],cwd=ROOT,
                               capture_output=True,text=True,check=True)
                inventory = schiff.atoms(receptor.read_text())
                self.assertGreaterEqual(len([n for r,n in inventory if r==2 and n.startswith('HZ')]),2)
                self.assertFalse(any(r==999 for r,n in inventory))
                self.assertFalse(schiff.links(receptor.read_text()))
                reaction = subprocess.run([str(ROOT/'test/schiff_repair_test'),str(receptor),str(sdf),'2'],
                                          cwd=ROOT,capture_output=True,text=True)
                self.assertEqual(reaction.returncode,0,reaction.stdout+reaction.stderr)


if __name__=='__main__':
    unittest.main()
