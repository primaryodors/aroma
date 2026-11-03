"""Generate a private viewer structure with RDKit, without database changes."""
import argparse
import sys


def generate(smiles, output):
    from rdkit import Chem
    from rdkit.Chem import AllChem

    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError("Invalid SMILES string.")
    molecule = Chem.AddHs(molecule)
    if not 1 <= molecule.GetNumAtoms() <= 512:
        raise ValueError("The viewer accepts structures with 1 to 512 atoms.")
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = 42
    if AllChem.EmbedMolecule(molecule, parameters) < 0:
        parameters.useRandomCoords = True
        if AllChem.EmbedMolecule(molecule, parameters) < 0:
            raise ValueError("RDKit could not generate 3D coordinates.")
    if AllChem.UFFHasAllMoleculeParams(molecule):
        AllChem.UFFOptimizeMolecule(molecule, maxIters=500)
    with Chem.SDWriter(output) as writer:
        writer.write(molecule)
    print("Generated 3D coordinates with RDKit.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("smiles")
    parser.add_argument("output")
    arguments = parser.parse_args()
    try:
        generate(arguments.smiles, arguments.output)
    except Exception as error:
        print(str(error), file=sys.stderr)
        sys.exit(2)
