"""Run an isolated OR2C1 ligand panel; docking scores are not activity labels."""
import argparse
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys

from rdkit import Chem
from rdkit.Chem import AllChem

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))

PANEL=[
    ('octanethiol','CCCCCCCCS','agonist'),
    ('nonanethiol','CCCCCCCCCS','agonist'),
    ('hexanethiol','CCCCCCS','agonist'),
    ('heptanethiol','CCCCCCCS','agonist'),
    ('1-octanol','CCCCCCCCO','nonagonist'),
    ('1-nonanol','CCCCCCCCCO','nonagonist'),
    ('nonanal','CCCCCCCCC=O','nonagonist'),
    ('nonane','CCCCCCCCC','unknown'),
]


def pose_geometry(path,molecule):
    if not path.is_file():
        return None
    points={line[12:16].strip():tuple(float(line[i:i+8]) for i in (30,38,46))
            for line in path.read_text().splitlines() if line.startswith('HETATM') and line[17:20]=='LIG'}
    conformer=molecule.GetConformer()
    errors=[]
    for bond in molecule.GetBonds():
        a,b=bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()
        names=[molecule.GetAtomWithIdx(i).GetSymbol()+str(i+1) for i in (a,b)]
        if not all(name in points for name in names):
            raise ValueError('Output pose has an incomplete ligand atom inventory')
        reference=(conformer.GetAtomPosition(a)-conformer.GetAtomPosition(b)).Length()
        errors.append(abs(math.dist(points[names[0]],points[names[1]])-reference))
    return {'max_bond_change_angstrom':max(errors,default=0)}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'tmp/or2c1-thiol-panel')
    parser.add_argument('--iterations',type=int,default=25)
    parser.add_argument('--resume',action='store_true',help='Keep successful runs and retry execution failures')
    parser.add_argument('--states',nargs='+',choices=['active','inactive'],default=['active','inactive'])
    args=parser.parse_args()
    work=args.output.resolve()
    work.mkdir(parents=True,exist_ok=True)
    pocket=json.loads((ROOT/'data/binding_pocket.json').read_text())['OR2C1']['pocket']
    environment=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
    results_path=work/'results.json'
    prior=json.loads(results_path.read_text()) if args.resume and results_path.is_file() else None
    if prior and prior['iterations']!=args.iterations:
        raise ValueError('Resume requires the same iteration count')
    rows=[row for row in prior['rows'] if row['exit_code']==0] if prior else []
    completed={(row['ligand'],row['state']) for row in rows}
    for name,smiles,activity in PANEL:
        molecule=Chem.AddHs(Chem.MolFromSmiles(smiles))
        if AllChem.EmbedMolecule(molecule,randomSeed=11973)<0:
            raise ValueError(f'Embedding failed: {name}')
        AllChem.UFFOptimizeMolecule(molecule,maxIters=1000)
        sdf=work/(name+'.sdf')
        with Chem.SDWriter(str(sdf)) as writer:
            writer.write(molecule)
        for state in args.states:
            if (name,state) in completed:
                continue
            prefix=work/(name+'.'+state)
            dock=Path(str(prefix)+'.dock')
            pose=Path(str(prefix)+'.model1.pdb')
            config=Path(str(prefix)+'.config')
            for artifact in (dock,pose):
                artifact.unlink(missing_ok=True)
            config.write_text(
                f'PROT {ROOT}/pdbs/OR2/OR2C1.{state}.pdb\nLIG {sdf}\n'
                f'CEN RES {pocket}\nSIZE 8\nSEARCH BB\nPOSE 1\nITERS {args.iterations}\n'
                f'THREADS 1\nRETRY 0\nNOFAIL\nOUT {dock}\nOUTPDB 1 {prefix}.model%o.pdb\n')
            command=[str(ROOT/'bin/aromadock'),str(config),'--ns','-t','1']
            with Path(str(prefix)+'.log').open('w') as log:
                try:
                    result=subprocess.run(command,cwd=ROOT,env=environment,stdout=log,stderr=subprocess.STDOUT,timeout=120)
                    exit_code=result.returncode
                except subprocess.TimeoutExpired:
                    exit_code='timeout'
            text=dock.read_text() if dock.is_file() else ''
            counts=re.findall(r'^([0-9]+) pose\(s\) found[.]',text,re.M)
            totals=re.findall(r'^Total: ([-+0-9.eE]+)',text,re.M)
            geometry=pose_geometry(pose,molecule)
            if geometry is not None and geometry['max_bond_change_angstrom']>.02:
                raise ValueError(f'Covalent geometry changed: {name} {state}: {geometry}')
            rows.append({'ligand':name,'state':state,'reference_activity':activity,
                         'exit_code':exit_code,'accepted_poses':int(counts[-1]) if counts else 0,
                         'energy_units':'kJ/mol',
                         'reported_energy':float(totals[0]) if totals else None,
                         'primary_target':next((line[5:] for line in text.splitlines() if line.startswith('Pri: ')),None),
                         'geometry':geometry,'dock':str(dock)})
            results_path.write_text(json.dumps({
                'purpose':'Mechanistic docking smoke test; scores do not establish agonism',
                'iterations':args.iterations,'rows':rows},indent=2)+'\n')
            print(name,state,'reference='+activity,'accepted='+str(rows[-1]['accepted_poses']),flush=True)
    return int(any(row['exit_code']!=0 for row in rows))


if __name__=='__main__':
    raise SystemExit(main())
