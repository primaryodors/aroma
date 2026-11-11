"""Resolve configured metal sites independently of ligand selection for tests."""
import json
import math
from pathlib import Path
import re

from data import protutils as pu

DONORS = {'C':('SG',), 'M':('SD',), 'H':('ND1','NE2'), 'K':('NZ',)}


def configured_metal_sites(receptor, policy_file=None):
    policy_file = policy_file or Path(__file__).with_name('binding_pocket.json')
    policies = json.loads(Path(policy_file).read_text())
    group = next((value for pattern,value in policies.items()
                  if pattern == receptor or (pattern.endswith('*') and receptor.startswith(pattern[:-1]))
                  or (not pattern.endswith('*') and re.search(pattern,receptor))), {})
    definitions = []
    for rule in [group, *group.get('odorophores',{}).values()]:
        values = rule.get('mcoord',[])
        definitions.extend([values] if isinstance(values,str) else values)
    result = []
    sequence = pu.prots[receptor]['sequence']
    for definition in dict.fromkeys(definitions):
        if definition.startswith('#'):
            continue
        element,charge,*tokens = definition.split()
        residues,donors = [],[]
        for token in tokens:
            match = re.fullmatch(r'([A-Z]*)([0-9]+(?:\.[0-9]+)?)',token)
            if not match:
                raise ValueError('Invalid configured metal residue')
            letters,position = match.groups()
            number = pu.resno_from_bw(receptor,position) if '.' in position else int(position)
            if not number or not 1 <= number <= len(sequence) or (letters and sequence[number-1] not in letters):
                break
            amino = sequence[number-1]
            if amino not in DONORS:
                raise ValueError('No donor atom definition for the configured metal residue')
            residues.append(number)
            donors.append({'residue':number,'amino':amino,'atoms':DONORS[amino]})
        else:
            if len(residues)<2 or len(set(residues))!=len(residues):
                raise ValueError('Configured metal site requires distinct coordinating residues')
            result.append({'element':element,'charge':int(charge),'residues':residues,'donors':donors})
    return result


def inspect_metal_sites(path, sites):
    atoms,ions = {},[]
    for line in Path(path).read_text().splitlines():
        if not line.startswith(('ATOM  ','HETATM')):
            continue
        name,number = line[12:16].strip(),int(line[22:26])
        point = tuple(float(line[i:i+8]) for i in (30,38,46))
        if not all(math.isfinite(x) for x in point):
            raise ValueError('Nonfinite coordinate after metal placement')
        atoms[number,name]=point
        if line[76:78].strip().lower() in {site['element'].lower() for site in sites}:
            ions.append((line[76:78].strip(),point))
    result,used = [],set()
    for site in sites:
        candidates=[]
        for index,(element,point) in enumerate(ions):
            if index in used or element.lower()!=site['element'].lower():
                continue
            distances=[]
            for donor in site['donors']:
                options=[(math.dist(point,atoms[donor['residue'],name]),name)
                         for name in donor['atoms'] if (donor['residue'],name) in atoms]
                if not options:
                    break
                distance,name=min(options)
                distances.append({'residue':donor['residue'],'amino':donor['amino'],
                                  'atom':name,'distance_angstrom':distance})
            else:
                candidates.append((sum(d['distance_angstrom'] for d in distances),index,point,distances))
        if not candidates:
            raise ValueError('Configured metal or donor atoms are absent from the saved PDB')
        _,index,point,distances=min(candidates)
        used.add(index)
        if site['element']=='Cu' and all(d['atom'] in ('SG','SD') for d in distances):
            if any(not 1.8<=d['distance_angstrom']<=2.8 for d in distances):
                raise ValueError('Saved Cu site failed its sulfur-distance check')
        record={'element':site['element'],'charge':site['charge'],'donors':distances}
        cysteines=[atoms[d['residue'],d['atom']] for d in distances if d['amino']=='C']
        if len(cysteines)==2:
            a,b=[tuple(p[i]-point[i] for i in range(3)) for p in cysteines]
            cosine=sum(x*y for x,y in zip(a,b))/math.sqrt(sum(x*x for x in a)*sum(x*x for x in b))
            theta=math.degrees(math.acos(max(-1,min(1,cosine))))
            record['cysteine_metal_cysteine_degrees']=theta
            if site['element']=='Cu' and not 80<=theta<=180:
                raise ValueError('Saved Cu site failed its cysteine-angle check')
        result.append(record)
    return result
