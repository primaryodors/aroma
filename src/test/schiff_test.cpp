#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "../classes/protein.h"

using namespace std;

Molecule* mols[4];

int main(int argc, char** argv)
{
    mols[0] = new Molecule("benzaldehyde", "c1ccccc1C=O");
    mols[1] = new Molecule("butylamine", "CCCC[NH3+]");
    mols[2] = mols[3] = nullptr;

    mols[1]->create_Schiff_base(mols[0]);

    FILE* fp = fopen("oh_schiff.sdf", "w");
    if (!fp)
    {
        cerr << "FAILED to open SDF output file for writing." << endl;
        return 1;
    }
    mols[0]->save_sdf(fp, &mols[1]);
    fclose(fp);
    fp = fopen("oh_schiff.pdb", "w");
    if (!fp)
    {
        cerr << "FAILED to open PDB output file for writing." << endl;
        return 1;
    }
    mols[0]->save_pdb(fp, 0, false);
    mols[1]->save_pdb(fp, mols[0]->get_atom_count(), true);
    fclose(fp);
    cout << "Wrote output files." << endl;

    Protein p("test");
    p.add_sequence("TIVKTCALTILTA");
    p.make_helix(1, p.get_end_resno(), ALPHA_PHI, ALPHA_PSI);
    mols[2] = new Molecule("aldehydeMNA", "O=CC(C)CCCCCCCCC");
    AminoAcid* aa = p.get_residue(4);
    aa->create_Schiff_base(mols[2]);

    Atom *NZ = aa->get_atom("NZ");
    Bond *b = NZ->get_bond_by_idx(2);
    if (!b)
    {
        cerr << "Schiff bond is missing." << endl;
        return 2;
    }
    Atom *C = b->atom2;

    cout << "Schiff atoms are: " << NZ->residue << ":" << NZ->name << ", " << C->residue << ":" << C->name << endl;
    cout << "Bond distance: " << NZ->distance_to(C) << endl;

    Bond **bb = ((Molecule*)aa)->get_rotatable_bonds();
    int i;
    if (!bb)
    {
        cerr << "No rotatable bonds." << endl;
        return 2;
    }
    for (i=0; bb[i]; i++)
    {
        cout << "Rotable: " << bb[i]->atom1->name << "-" << bb[i]->atom2->name << endl;
    }

    fp = fopen("well_schiff.pdb", "w");
    if (!fp)
    {
        cerr << "FAILED to open PDB output file for writing." << endl;
        return 1;
    }
    p.save_pdb(fp, mols[2]);
    p.end_pdb(fp);
    fclose(fp);
    cout << "Wrote output file." << endl;

    return 0;
}