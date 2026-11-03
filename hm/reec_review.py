#!/usr/bin/env python3
"""Independent checks of serialized REEC PDB geometry and contact reports."""
import argparse
import json
import math
from pathlib import Path

HEAVY = {
    'ALA': 'CB', 'ARG': 'CB CG CD NE CZ NH1 NH2', 'ASN': 'CB CG OD1 ND2',
    'ASP': 'CB CG OD1 OD2', 'CYS': 'CB SG', 'GLN': 'CB CG CD OE1 NE2',
    'GLU': 'CB CG CD OE1 OE2', 'GLY': '', 'HIS': 'CB CG ND1 CD2 CE1 NE2',
    'ILE': 'CB CG1 CG2 CD1', 'LEU': 'CB CG CD1 CD2', 'LYS': 'CB CG CD CE NZ',
    'MET': 'CB CG SD CE', 'PHE': 'CB CG CD1 CD2 CE1 CE2 CZ', 'PRO': 'CB CG CD',
    'SER': 'CB OG', 'THR': 'CB OG1 CG2', 'TRP': 'CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2',
    'TYR': 'CB CG CD1 CD2 CE1 CE2 CZ OH', 'VAL': 'CB CG1 CG2',
}


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]


def angle(a, b, c):
    u, v = sub(a, b), sub(c, b)
    scale = math.sqrt(dot(u, u) * dot(v, v))
    if scale < 1e-12:
        return 0
    return math.degrees(math.acos(max(-1, min(1, dot(u, v) / scale))))


def dihedral(a, b, c, d):
    u, v, w = sub(b, a), sub(c, b), sub(d, c)
    n, m = cross(u, v), cross(v, w)
    norm = math.sqrt(dot(v, v))
    if norm < 1e-12:
        return 0
    return math.degrees(math.atan2(dot(cross(n, m), v) / norm, dot(n, m)))


def inspect(path):
    path = Path(path)
    residues, serials, links, all_atoms = {}, {}, set(), {}
    for line in path.read_text().splitlines():
        if line.startswith(('ATOM  ', 'HETATM')):
            point = tuple(float(line[i:i+8]) for i in (30, 38, 46))
            all_atoms[int(line[22:26]), line[12:16].strip()] = point
        if line.startswith('ATOM  '):
            no = int(line[22:26])
            name = line[12:16].strip()
            residue = residues.setdefault(no, {'type': line[17:20], 'atoms': {}})
            if name in residue['atoms']:
                raise ValueError(f'Duplicate atom {no}:{name}')
            point = tuple(float(line[i:i+8]) for i in (30, 38, 46))
            if not all(math.isfinite(value) for value in point):
                raise ValueError('Nonfinite coordinates')
            residue['atoms'][name] = point
            serials[int(line[6:11])] = no, name
        elif line.startswith('CONECT'):
            fields = [int(line[i:i+5]) for i in range(6, len(line), 5) if line[i:i+5].strip()]
            if fields:
                links.update(tuple(sorted((fields[0], other))) for other in fields[1:])
    if not residues:
        raise ValueError('No protein atoms')
    failures, distances, angles, omegas = [], [], [], []
    for no, residue in sorted(residues.items()):
        atoms = residue['atoms']
        required = set(('N CA C O ' + HEAVY[residue['type']]).split())
        if no == max(residues):
            required.add('OXT')
        if set(atoms) != required:
            failures.append(f'{no}: atom inventory differs: {sorted(set(atoms) ^ required)}')
        if not required.issubset(atoms):
            continue
        for a, b, ideal in [('N', 'CA', 1.458), ('CA', 'C', 1.525), ('C', 'O', 1.231)]:
            error = abs(math.dist(atoms[a], atoms[b]) - ideal)
            distances.append(error)
            if error > .10:
                failures.append(f'{no}:{a}-{b} bond error {error:.3f}')
        if 'OXT' in atoms:
            error = abs(math.dist(atoms['C'], atoms['OXT']) - 1.25)
            distances.append(error)
            if error > .10:
                failures.append(f'{no}: terminal carboxyl bond error {error:.3f}')
        if 'CB' in atoms:
            volume = dot(sub(atoms['N'], atoms['CA']),
                         cross(sub(atoms['C'], atoms['CA']), sub(atoms['CB'], atoms['CA'])))
            if volume <= 0:
                failures.append(f'{no}: alpha chirality inversion')
        if residue['type'] in ('ILE', 'THR'):
            branch = 'CG1' if residue['type'] == 'ILE' else 'OG1'
            volume = dot(sub(atoms['CA'], atoms['CB']),
                         cross(sub(atoms[branch], atoms['CB']), sub(atoms['CG2'], atoms['CB'])))
            if volume <= 0:
                failures.append(f'{no}: beta chirality inversion')
        if no + 1 not in residues:
            continue
        following = residues[no+1]['atoms']
        error = abs(math.dist(atoms['C'], following['N']) - 1.329)
        distances.append(error)
        if error > .08:
            failures.append(f'{no}-{no+1}: peptide bond error {error:.3f}')
        for actual, ideal in [(angle(atoms['CA'], atoms['C'], following['N']), 116.2),
                              (angle(atoms['C'], following['N'], following['CA']), 121.7)]:
            angles.append(abs(actual - ideal))
            if abs(actual - ideal) > 10:
                failures.append(f'{no}: peptide angle {actual:.2f}')
        omega = dihedral(atoms['CA'], atoms['C'], following['N'], following['CA'])
        omegas.append(abs(abs(omega) - 180))
        if abs(abs(omega) - 180) > 20:
            failures.append(f'{no}: peptide omega {omega:.2f}')
    sulfur_links = []
    for first, second in sorted(links):
        if first not in serials or second not in serials:
            continue
        (a, atom1), (b, atom2) = serials[first], serials[second]
        if atom1 != 'SG' or atom2 != 'SG':
            continue
        aa, bb = residues[a]['atoms'], residues[b]['atoms']
        distance = math.dist(aa['SG'], bb['SG'])
        torsion = dihedral(aa['CB'], aa['SG'], bb['SG'], bb['CB'])
        css = [angle(aa['CB'], aa['SG'], bb['SG']), angle(bb['CB'], bb['SG'], aa['SG'])]
        sulfur_links.append({'residues': [a, b], 'distance': distance,
                             'angles': css, 'torsion': torsion})
        if not 1.90 <= distance <= 2.20 or abs(abs(torsion) - 90) > 30 or any(abs(v - 104) > 15 for v in css):
            failures.append(f'{a}-{b}: disulfide geometry outside review bounds')
    report = json.loads(Path(str(path) + '.json').read_text())
    for contact in report['contacts']:
        a = all_atoms[contact['res1'], contact['atom1']]
        b = all_atoms[contact['res2'], contact['atom2']]
        distance = math.dist(a, b)
        if abs(distance - contact['distance']) > .002:
            failures.append('Serialized contact distance differs from report')
        if contact['hard'] and not contact['min'] - .052 <= distance <= contact['max'] + .052:
            failures.append('Serialized mandatory contact outside bounds')
    return {'model': str(path), 'residues': len(residues), 'valid': not failures,
            'max_backbone_bond_error': max(distances, default=0),
            'max_peptide_angle_error': max(angles, default=0),
            'max_omega_error': max(omegas, default=0), 'disulfides': sulfur_links,
            'failures': failures, 'optimizer_report': report}


def plot(reviews, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    figure, axes = plt.subplots(1, 3, figsize=(13, 4.3), constrained_layout=True)
    labels = [Path(r['model']).name.split('.')[0] for r in reviews]
    positions = list(range(len(labels)))
    initial = [r['optimizer_report']['initial_severe_clashes'] for r in reviews]
    final = [r['optimizer_report']['severe_clashes'] for r in reviews]
    axes[0].bar([v-.18 for v in positions], initial, .36, label='Initial conformer', color='#dc7757')
    axes[0].bar([v+.18 for v in positions], final, .36, label='Relaxed model', color='#238d88')
    for i, count in enumerate(final):
        axes[0].text(i+.18, count + max(initial, default=1)*.025, str(count), ha='center')
    axes[0].set_ylabel('Severe heavy-atom overlaps')
    axes[0].legend(frameon=False)
    axes[1].bar(positions, [r['max_backbone_bond_error'] for r in reviews], color='#238d88')
    axes[1].axhline(.08, linestyle='--', color='#dc7757', label='REEC peptide bond tolerance')
    axes[1].set_ylabel('Largest backbone bond error (Å)')
    axes[1].legend(frameon=False, fontsize=8)
    axes[2].bar(positions, [r['optimizer_report']['template_ca_rmsd'] for r in reviews], color='#465c89')
    axes[2].set_ylabel('Template CA displacement RMS (Å)')
    for axis in axes:
        axis.set_xticks(positions, labels)
        axis.spines[['top', 'right']].set_visible(False)
    figure.suptitle('REEC structural validation — relaxation with contact restraints', fontsize=14)
    figure.savefig(output, dpi=180)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('models', type=Path, nargs='+')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--plot', type=Path)
    args = parser.parse_args()
    reviews = [inspect(path) for path in args.models]
    args.output.write_text(json.dumps(reviews, indent=2) + '\n')
    if args.plot:
        plot(reviews, args.plot)
    print(json.dumps([{k: v for k, v in r.items() if k not in ('optimizer_report', 'disulfides')}
                      for r in reviews], indent=2))
    return 0 if all(r['valid'] for r in reviews) else 2


if __name__ == '__main__':
    raise SystemExit(main())
