"""Relax an existing receptor model while restraining its active backbone."""

import argparse
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from data.disulfides import pdb_disulfides


def relax_model(source, destination, iterations=600, backbone_sigma=0.2):
    from modeller import Environ, Selection, features, forms, log, physical
    from modeller.optimizers import ConjugateGradients
    from modeller.scripts import complete_pdb

    source = Path(source).resolve()
    destination = Path(destination).resolve()
    if source == destination:
        raise ValueError("Choose a separate output file for the candidate model.")
    if iterations < 1 or backbone_sigma <= 0:
        raise ValueError("Iterations and backbone sigma must be positive.")
    links = pdb_disulfides(str(source))
    log.minimal()
    env = Environ()
    env.io.hetatm = False
    env.io.hydrogen = False
    env.libs.topology.read(file="$(LIB)/top_heav.lib")
    env.libs.parameters.read(file="$(LIB)/par.lib")
    env.edat.dynamic_sphere = True
    model = complete_pdb(env, str(source))
    for first, second in links:
        model.patch(residue_type="DISU", residues=(model.residues[f"{first}:A"],
                                                   model.residues[f"{second}:A"]))
    atoms = Selection(model)
    model.restraints.make(atoms, restraint_type="stereo", spline_on_site=False)
    for atom in model.atoms:
        if atom.name not in ("N", "CA", "C", "O"):
            continue
        for axis, coordinate in (("x", atom.x), ("y", atom.y), ("z", atom.z)):
            model.restraints.add(forms.Gaussian(
                group=physical.xy_distance,
                feature=getattr(features, axis + "_coordinate")(atom),
                mean=coordinate, stdev=backbone_sigma))
    initial = atoms.energy()[0]
    ConjugateGradients(output="REPORT").optimize(
        atoms, max_iterations=iterations, min_atom_shift=0.0001)
    final = atoms.energy()[0]
    if final >= initial:
        raise RuntimeError("Relaxation did not improve the model objective.")

    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="aroma-relax-") as directory:
        work = Path(directory)
        heavy, output = work / "heavy.pdb", work / "candidate.pdb"
        model.write(file=str(heavy))
        code = [f'LOAD "{source}" A I', f'LOAD "{heavy}" A A',
                "BWCOPY I A", "UNCHAIN I", "STRAND A"]
        code.extend(line for line in source.read_text().splitlines()
                    if line.startswith("REMARK"))
        code.append(f"REMARK 265 RESTRAINED_RELAXATION: objective {initial:.3f} -> {final:.3f}; "
                    f"backbone sigma {backbone_sigma:g} A; {iterations} iterations")
        code.extend(f"CONECT {a} SG {b} SG" for a, b in links)
        code.extend(("HYDRO", f'SAVE "{output}"', "EXIT"))
        script = work / "export.phew"
        script.write_text("\n".join(code) + "\n")
        subprocess.run([str(ROOT / "bin/phew"), str(script)], cwd=ROOT, check=True)
        if not output.is_file() or sorted(pdb_disulfides(str(output))) != sorted(links):
            raise RuntimeError("Candidate export did not preserve declared disulfides.")
        with tempfile.NamedTemporaryFile(dir=destination.parent, delete=False) as stream:
            temporary = Path(stream.name)
            stream.write(output.read_bytes())
        os.replace(temporary, destination)
    return initial, final


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--iterations", type=int, default=600)
    parser.add_argument("--backbone-sigma", type=float, default=0.2)
    args = parser.parse_args()
    initial, final = relax_model(args.source, args.destination,
                                args.iterations, args.backbone_sigma)
    print(f"Model objective: {initial:.3f} -> {final:.3f}")
