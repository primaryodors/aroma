"""Read saved docking poses and measure the requested site-specific contacts."""
import math
from pathlib import Path

from rdkit import Chem


def angle(a, b, c):
    u = [x-y for x, y in zip(a, b)]
    v = [x-y for x, y in zip(c, b)]
    scale = math.sqrt(sum(x*x for x in u) * sum(x*x for x in v))
    if scale == 0:
        return None
    return math.degrees(math.acos(max(-1., min(1., sum(x*y for x, y in zip(u, v))/scale))))


def hydroxy_atom_labels(molecule):
    """Locate the 3-hydroxy oxygen through the carboxyl-C2-C3 molecular graph."""
    labels = set()
    for carboxyl in molecule.GetAtoms():
        if carboxyl.GetAtomicNum() != 6:
            continue
        oxygens = [atom for atom in carboxyl.GetNeighbors() if atom.GetAtomicNum() == 8]
        if len(oxygens) != 2:
            continue
        if not any(molecule.GetBondBetweenAtoms(carboxyl.GetIdx(), atom.GetIdx()).GetBondType() == Chem.BondType.DOUBLE for atom in oxygens):
            continue
        for alpha in carboxyl.GetNeighbors():
            if alpha.GetAtomicNum() != 6:
                continue
            for beta in alpha.GetNeighbors():
                if beta.GetAtomicNum() != 6 or beta.GetIdx() == carboxyl.GetIdx():
                    continue
                for oxygen in beta.GetNeighbors():
                    if oxygen.GetAtomicNum() != 8:
                        continue
                    heavy = [atom for atom in oxygen.GetNeighbors() if atom.GetAtomicNum() > 1]
                    if len(heavy) == 1 and molecule.GetBondBetweenAtoms(beta.GetIdx(), oxygen.GetIdx()).GetBondType() == Chem.BondType.SINGLE:
                        labels.add(f'O{oxygen.GetIdx()+1}')
    return sorted(labels)


def read_pose(path):
    protein, ligand, copper = {}, {}, []
    for line in Path(path).read_text().splitlines():
        if not line.startswith(('ATOM  ', 'HETATM')):
            continue
        number = int(line[22:26].strip() or '0')
        name, residue = line[12:16].strip(), line[17:20].strip()
        element = line[76:78].strip() if len(line) >= 78 else name.lstrip('0123456789')[0]
        point = tuple(float(line[i:i+8]) for i in (30, 38, 46))
        if not all(math.isfinite(value) for value in point):
            raise ValueError('Nonfinite pose coordinate')
        record = {'name': name, 'element': element, 'point': point}
        if element.lower() == 'cu':
            copper.append(record)
        elif line.startswith('ATOM  '):
            protein[number, name] = record
        elif residue == 'LIG':
            ligand[name] = record
    return protein, ligand, copper


def copper_check(protein, ligand, copper):
    if not copper or (199, 'SG') not in protein or (200, 'SG') not in protein:
        return {'coordination_pass': False, 'reason': 'Copper or a required cysteine sulfur is absent'}
    possibilities = []
    for metal in copper:
        location = metal['point']
        sulfur = [protein[r, 'SG']['point'] for r in (199, 200)]
        distances = [math.dist(location, point) for point in sulfur]
        coordination_angle = angle(sulfur[0], location, sulfur[1])
        ligand_distances = [math.dist(location, atom['point']) for atom in ligand.values() if atom['element'] == 'S']
        nearest = min(ligand_distances) if ligand_distances else None
        possibilities.append({'cys199_cu_angstrom': distances[0], 'cys200_cu_angstrom': distances[1],
                              'cys199_cu_cys200_degrees': coordination_angle,
                              'ligand_sulfur_cu_angstrom': nearest,
                              'coordination_pass': all(1.8 <= value <= 2.8 for value in distances)
                                  and coordination_angle is not None and 80 <= coordination_angle <= 180,
                              'ligand_copper_contact': nearest is not None and nearest <= 3.0})
    return min(possibilities, key=lambda value: abs(value['cys199_cu_angstrom']-2.25)+abs(value['cys200_cu_angstrom']-2.25))


def asparagine_check(protein, ligand, labels):
    pairs = []
    for label in labels:
        if label not in ligand:
            continue
        oxygen = ligand[label]['point']
        for name in ('OD1', 'ND2'):
            if (208, name) not in protein:
                continue
            site = protein[208, name]['point']
            distance = math.dist(site, oxygen)
            hydrogen_bonds = []
            if name == 'OD1':
                donors = [(oxygen, atom['point'], site) for atom in ligand.values()
                          if atom['element'] == 'H' and math.dist(atom['point'], oxygen) <= 1.25]
            else:
                donors = [(site, atom['point'], oxygen) for (residue, atom_name), atom in protein.items()
                          if residue == 208 and atom['element'] == 'H' and math.dist(atom['point'], site) <= 1.25]
            for donor, hydrogen, acceptor in donors:
                theta = angle(donor, hydrogen, acceptor)
                if 2.2 <= distance <= 3.5 and math.dist(hydrogen, acceptor) <= 2.6 and theta is not None and theta >= 120:
                    hydrogen_bonds.append({'hydrogen_acceptor_angstrom': math.dist(hydrogen, acceptor),
                                           'donor_hydrogen_acceptor_degrees': theta})
            pairs.append({'ligand_3_hydroxy_atom': label, 'asn208_atom': name, 'distance_angstrom': distance,
                          'contact_pass': 2.2 <= distance <= 3.5, 'hydrogen_bond_pass': bool(hydrogen_bonds),
                          'hydrogen_bonds': hydrogen_bonds})
    if not pairs:
        return {'contact_pass': False, 'hydrogen_bond_pass': False,
                'reason': 'The 3-hydroxy oxygen or Asn208 side-chain atoms are absent'}
    best = min(pairs, key=lambda value: value['distance_angstrom'])
    return dict(best, all_contacts=pairs, contact_pass=any(p['contact_pass'] for p in pairs),
                hydrogen_bond_pass=any(p['hydrogen_bond_pass'] for p in pairs))


def check_pose(path, receptor, sdf_files):
    protein, ligand, copper = read_pose(path)
    result = {'pose': str(path), 'ligand_atoms': len(ligand)}
    if receptor == 'OR2T11':
        result['copper'] = copper_check(protein, ligand, copper)
    if receptor == 'OR52E8':
        label_sets = []
        for path in sdf_files:
            molecule = Chem.SDMolSupplier(str(path), removeHs=False)[0]
            if molecule is None:
                raise ValueError('Cannot parse a frozen hydroxy-acid SDF')
            labels = hydroxy_atom_labels(molecule)
            if len(labels) != 1:
                raise ValueError('Expected one chemically identified 3-hydroxy oxygen')
            label_sets.append(labels)
        if any(labels != label_sets[0] for labels in label_sets):
            raise ValueError('Stereoisomer SDFs use different 3-hydroxy atom indices')
        result['asn208'] = asparagine_check(protein, ligand, label_sets[0])
    return result
