"""Schiff-adduct identity and validation across receptor repair stages."""
from dataclasses import dataclass, replace
import math
from pathlib import Path

PREFIX = 'REMARK 950 SCHIFF '


@dataclass(frozen=True)
class Link:
    residue: int
    nitrogen: str
    ligand: int
    carbon: str
    neighbor: str
    oxygen: str
    hydrogen: str = '-'
    configuration: str = 'trans'

    def remark(self):
        return PREFIX + f'{self.residue} {self.nitrogen} {self.ligand} {self.carbon} {self.neighbor} {self.oxygen} {self.hydrogen} {self.configuration}'


def atoms(text):
    result = {}
    for line in text.splitlines():
        if line.startswith(('ATOM  ', 'HETATM')):
            key = (int(line[22:26]), line[12:16].strip())
            if key in result:
                if key[1].lstrip('0123456789').startswith('H') and not key[1].startswith('HZ'):
                    continue
                raise ValueError(f'Ambiguous Schiff atom identity: {key}')
            result[key] = (line, tuple(float(line[i:i+8]) for i in (30, 38, 46)))
    return result


def links(text):
    result = []
    for line in text.splitlines():
        if line.startswith(PREFIX):
            fields = line[len(PREFIX):].split()
            r, n, l, c, b, o = fields[:6]
            if len(fields)>8 or (len(fields)>7 and fields[7] not in ('cis','trans')):
                raise ValueError('Invalid Schiff configuration metadata')
            link = Link(int(r), n, int(l), c, b, o, fields[6] if len(fields)>6 else '-',fields[7] if len(fields)>7 else 'trans')
            if link not in result:
                result.append(link)
    return result


def topology(sdf):
    rows = Path(sdf).read_text().splitlines()
    if 'V2000' not in rows[3]:
        raise ValueError('Schiff identification requires a V2000 ligand reference')
    count, total = int(rows[3][:3]), int(rows[3][3:6])
    names = [rows[4+i][31:34].strip()+str(i+1) for i in range(count)]
    bonds = []
    for row in rows[4+count:4+count+total]:
        bonds.append((names[int(row[:3])-1], names[int(row[3:6])-1], float(row[6:].split()[0])))
    return names, bonds


def discover(path, sdf=None):
    text = Path(path).read_text()
    tagged = links(text)
    if tagged:
        return tagged
    inventory = atoms(text)
    by_serial = {int(line[6:11]):key for key,(line,_) in inventory.items()}
    explicit = set()
    for line in text.splitlines():
        if not line.startswith('CONECT'):
            continue
        serials = [int(line[i:i+5]) for i in range(6,len(line),5) if line[i:i+5].strip()]
        for other in serials[1:]:
            a,b = by_serial.get(serials[0]),by_serial.get(other)
            if a is None or b is None:
                continue
            if inventory[b][0].startswith('ATOM  '):
                a,b = b,a
            if a[1]=='NZ' and inventory[a][0][17:20]=='LYS' and inventory[b][0].startswith('HETATM') and b[1].startswith('C'):
                explicit.add((a,b))
    ligand = {key: value for key, value in inventory.items() if value[0].startswith('HETATM') and value[0][17:20]=='LIG'}
    if not ligand:
        return []
    if sdf is None:
        if explicit or any(1.05<=math.dist(p,q)<=1.65 for (r,n),(line,p) in inventory.items()
                           if n=='NZ' and line[17:20]=='LYS' for (lr,c),(_,q) in ligand.items() if c.startswith('C')):
            raise ValueError('A possible Schiff linkage requires a ligand SDF reference or saved Schiff metadata')
        return []
    names, bonds = topology(sdf)
    result = []
    for a, b, order in bonds:
        if order != 2:
            continue
        carbon, oxygen = (a, b) if a.startswith('C') and b.startswith('O') else (b, a)
        if not carbon.startswith('C') or not oxygen.startswith('O'):
            continue
        for (lr, name), (_, point) in ligand.items():
            if name != carbon or (lr, oxygen) in ligand:
                continue
            neighbors = [y if x==carbon else x for x,y,_ in bonds if carbon in (x,y) and oxygen not in (x,y)]
            heavy = [n for n in neighbors if not n.startswith('H')]
            if len(heavy)!=1 or not any(n.startswith('H') for n in neighbors):
                raise ValueError('Missing carbonyl oxygen has no supported aldehyde Schiff topology')
            candidates = [r for (r,n),(line,p) in inventory.items() if n=='NZ' and line[17:20]=='LYS'
                          and (((r,n),(lr,carbon)) in explicit or 1.05 <= math.dist(p, point) <= 1.65)]
            if len(candidates)!=1:
                raise ValueError(f'Ambiguous Schiff assignment for {carbon}: {len(candidates)} nearby lysine nitrogens')
            r = candidates[0]
            protons = [n for (number,n),(_,p) in inventory.items() if number==r and n.lstrip('0123456789').startswith('H')
                       and math.dist(p,inventory[r,'NZ'][1])<1.25]
            if len(protons)>1:
                raise ValueError('Schiff nitrogen has an unreconciled amine hydrogen inventory')
            charged = '+' in inventory[r,'NZ'][0][78:80]
            from reec_review import dihedral
            torsion = abs(dihedral(inventory[r,'CE'][1],inventory[r,'NZ'][1],point,inventory[lr,heavy[0]][1]))
            result.append(Link(r, 'NZ', lr, carbon, heavy[0], oxygen, protons[0] if protons else ('HZ3' if charged else '-'),
                               'trans' if torsion>90 else 'cis'))
    if len({(x.residue,x.nitrogen) for x in result})!=len(result):
        raise ValueError('Multiple Schiff carbons assigned to one nitrogen')
    if explicit - {((x.residue,x.nitrogen),(x.ligand,x.carbon)) for x in result}:
        raise ValueError('Protein-ligand connectivity is not a supported aldehyde Schiff linkage')
    return result


def annotate(path, adducts, sdf=None, ligand_number=None):
    """Restore identities and connectivity after serial changes; remove amine hydrogens."""
    if not adducts:
        return
    path = Path(path)
    text = path.read_text()
    inventory = atoms(text)
    updated = [replace(link, ligand=ligand_number) if ligand_number is not None else link for link in adducts]
    required = {(x.residue,x.nitrogen) for x in updated} | {(x.ligand,x.carbon) for x in updated} | {(x.ligand,x.neighbor) for x in updated}
    if not required.issubset(inventory):
        raise ValueError(f'Schiff atoms were lost during repair: {sorted(required-set(inventory))}')
    removed = set()
    for key,(line,point) in inventory.items():
        name = key[1].lstrip('0123456789')
        if name.startswith('H') and any(key[0]==x.residue and (name.startswith('HZ') or math.dist(point,inventory[x.residue,x.nitrogen][1])<1.25) for x in updated):
            removed.add(int(line[6:11]))
    rows = []
    for line in text.splitlines():
        if line.startswith(PREFIX):
            continue
        if line.startswith(('ATOM  ','HETATM')) and int(line[6:11]) in removed:
            continue
        if line.startswith('CONECT'):
            serials = [int(line[i:i+5]) for i in range(6,len(line),5) if line[i:i+5].strip()]
            if serials[0] in removed:
                continue
            serials = [n for n in serials if n not in removed]
            if len(serials)<2:
                continue
            line = 'CONECT'+''.join(f'{n:5d}' for n in serials)
        if line.strip()=='END':
            continue
        rows.append(line)
    serial = {key:int(value[0][6:11]) for key,value in inventory.items()}
    for link in updated:
        rows.append(link.remark())
        rows.append(f'CONECT{serial[link.residue,link.nitrogen]:5d}{serial[link.ligand,link.carbon]:5d}{serial[link.ligand,link.carbon]:5d}')
        if link.hydrogen!='-':
            n = inventory[link.residue,link.nitrogen][1]
            ends = [inventory[link.residue,'CE'][1],inventory[link.ligand,link.carbon][1]]
            direction = [-sum((point[k]-n[k])/math.dist(point,n) for point in ends) for k in range(3)]
            length = math.sqrt(sum(v*v for v in direction))
            if length<.1:
                raise ValueError('Schiff nitrogen has degenerate bond geometry')
            point = [n[k]+1.01*direction[k]/length for k in range(3)]
            last = max(serial.values())+1
            serial[link.residue,link.hydrogen] = last
            rows.append(f'ATOM  {last:5d} {link.hydrogen:>4s} LYS A{link.residue:4d}    '+''.join(f'{v:8.3f}' for v in point)+'  1.00  0.00           H')
        for i,line in enumerate(rows):
            if line.startswith('ATOM  ') and (int(line[22:26]),line[12:16].strip())==(link.residue,link.nitrogen):
                rows[i] = line.ljust(80)[:78]+('1+' if link.hydrogen!='-' else '  ')
    if sdf:
        _, bonds = topology(sdf)
        for link in updated:
            for a,b,order in bonds:
                if (link.ligand,a) in serial and (link.ligand,b) in serial:
                    rows.append('CONECT'+f'{serial[link.ligand,a]:5d}'+f'{serial[link.ligand,b]:5d}'*max(1,int(order)))
    path.write_text('\n'.join(rows)+'\nEND\n')


def failures(path, require_hydrogens=False):
    text = Path(path).read_text()
    adducts = links(text)
    if not adducts:
        return []
    inventory = atoms(text)
    result = []
    for link in adducts:
        try:
            ce = inventory[link.residue,'CE'][1]
            n = inventory[link.residue,link.nitrogen][1]
            c = inventory[link.ligand,link.carbon][1]
            b = inventory[link.ligand,link.neighbor][1]
        except KeyError:
            result.append('Schiff atom inventory is incomplete')
            continue
        distance = math.dist(n,c)
        if not 1.18 <= distance <= 1.38:
            result.append(f'Schiff C=N distance {distance:.3f} A is outside 1.18..1.38')
        for a,center,z in ((ce,n,c),(n,c,b)):
            u,v = [x-y for x,y in zip(a,center)], [x-y for x,y in zip(z,center)]
            angle = math.degrees(math.acos(max(-1,min(1,sum(x*y for x,y in zip(u,v))/(math.dist(a,center)*math.dist(z,center))))))
            if not 110 <= angle <= 130:
                result.append(f'Schiff bond angle {angle:.2f} degrees is outside 110..130')
        from reec_review import dihedral
        torsion = abs(dihedral(ce,n,c,b))
        if min(torsion,abs(180-torsion))>10:
            result.append(f'Schiff planarity error exceeds 10 degrees: {torsion:.2f}')
        if (torsion>90)!=(link.configuration=='trans'):
            result.append('Schiff imine configuration changed')
        carbon_protons = [p for (r,name),(_,p) in inventory.items() if r==link.ligand
                          and name.startswith('H') and .9<=math.dist(p,c)<=1.2]
        if len(carbon_protons)!=1:
            result.append('Schiff carbon hydrogen inventory changed')
        elif min(abs(dihedral(n,c,b,carbon_protons[0])),abs(180-abs(dihedral(n,c,b,carbon_protons[0]))))>10:
            result.append('Schiff carbon hydrogen is outside the imine plane')
        if (link.ligand,link.oxygen) in inventory:
            result.append('Departed Schiff oxygen reappeared')
        protons = [(name,point) for (r,name),(_,point) in inventory.items() if r==link.residue
                   and name.lstrip('0123456789').startswith('H') and (name.startswith('HZ') or math.dist(point,n)<1.25)]
        if len(protons)>(0 if link.hydrogen=='-' else 1):
            result.append('Excess amine hydrogens reappeared on Schiff nitrogen')
        if protons and (protons[0][0]!=link.hydrogen or not .90<=math.dist(protons[0][1],n)<=1.12):
            result.append('Schiff nitrogen hydrogen identity or distance changed')
        if require_hydrogens:
            if len(protons)!=(0 if link.hydrogen=='-' else 1):
                result.append('Schiff nitrogen protonation state changed')
            charge = inventory[link.residue,link.nitrogen][0][78:80].strip()
            if charge!=('1+' if link.hydrogen!='-' else ''):
                result.append('Schiff nitrogen formal charge changed')
    return result


def restore(path, source):
    """Transfer environment connectivity by atom identity after rigid postprocessing."""
    adducts = links(source)
    if not adducts:
        return
    path = Path(path)
    source_atoms, target_atoms = atoms(source), atoms(path.read_text())
    original = {int(line[6:11]):key for key,(line,_) in source_atoms.items()}
    current = {key:int(line[6:11]) for key,(line,_) in target_atoms.items()}
    rows = [line for line in path.read_text().splitlines() if line.strip()!='END']
    for line in source.splitlines():
        if not line.startswith('CONECT'):
            continue
        serials = [int(line[i:i+5]) for i in range(6,len(line),5) if line[i:i+5].strip()]
        keys = [original.get(number) for number in serials]
        if any(key not in current for key in keys):
            continue
        if any(source_atoms[key][0].startswith('HETATM') for key in keys):
            rows.append('CONECT'+''.join(f'{current[key]:5d}' for key in keys))
    path.write_text('\n'.join(rows)+'\nEND\n')
    annotate(path, adducts)


def receptor_only(path, adducts):
    """Remove the repaired adduct before docking an unreacted ligand again."""
    if not adducts:
        return
    path = Path(path)
    rows = path.read_text().splitlines()
    residues = {x.ligand for x in adducts}
    removed = {int(line[6:11]) for line in rows if line.startswith('HETATM') and line[17:20]=='LIG' and int(line[22:26]) in residues}
    output = []
    for line in rows:
        if line.startswith(PREFIX) or line.strip()=='END':
            continue
        if line.startswith('HETATM') and int(line[6:11]) in removed:
            continue
        if line.startswith('CONECT'):
            serials = [int(line[i:i+5]) for i in range(6,len(line),5) if line[i:i+5].strip()]
            if any(n in removed for n in serials):
                continue
        output.append(line)
    path.write_text('\n'.join(output)+'\nEND\n')
