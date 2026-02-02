
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "../classes/protein.h"

using namespace std;

int main(int argc, char** argv)
{
    Molecule amine("amine");
    Molecule aldehyde("aldehyde");

    amine.from_smiles("CCCC[NH3+]");        // Phee-yoo!
    aldehyde.from_smiles("C=CC=O");         // Yuck!

    Atom* N5 = amine.get_atom("N5");
    Bond* b = N5->get_bond_by_idx(0);
    b->rotate(M_PI);                        // makes things easier to see

    Molecule* water = amine.create_Schiff_base(&aldehyde);

    FILE* fp = fopen("tmp/schiff.sdf", "w");
    Molecule* ligs[3];
    ligs[0] = &aldehyde;
    ligs[1] = water;
    ligs[2] = nullptr;
    amine.save_sdf(fp, ligs);
    fclose(fp);
    cout << "Wrote output SDF file." << endl;

    Protein p("peptide");
    Molecule ligand("ligand");

    p.add_sequence("LLILFKNLQLVLIST");
    p.make_helix(1, p.get_end_resno(), ALPHA_PHI, ALPHA_PSI);
    AminoAcid* aaK = p.get_residue(6);
    ligand.from_smiles("c1ccccc1C=O");
    ligand.create_Schiff_base(aaK);

    // This should fail.
    Point lwas = ligand.get_barycenter();
    ligand.move(Point(100,0,0));
    Point lmovd = ligand.get_barycenter();
    ligand.recenter(Point(53,81,126));
    Point lrcnd = ligand.get_barycenter();

    if (lwas.get_3d_distance(lmovd) > 0.1) cout << "ERROR! LIGAND MOVED!" << endl;
    if (lrcnd.get_3d_distance(lmovd) > 0.1) cout << "ERROR! LIGAND RECENTERED!" << endl;

    Atom *CA = aaK->get_atom("CA"), *CB = aaK->get_atom("CB");
    b = CA->get_bond_between(CB);
    b->rotate(square);
    Point lflxd = ligand.get_barycenter();

    if (lflxd.get_3d_distance(lrcnd) < 0.1) cout << "ERROR! LIGAND DID NOT FLEX!" << endl;

    Atom *NZ = aaK->get_atom("NZ");
    Atom *aldC = NZ->is_bonded_to("C", 2);
    if (!aldC)
    {
        cout << "ERROR! FAILED TO BOND ALDEHYDE CARBON!" << endl;
    }
    else
    {
        float r = aldC->distance_to(NZ);
        if (r < 1.2 || r > 1.35) cout << "ERROR! NZ AND ALDEHYDE C ARE " << r << "A APART! SHOULD BE 1.28!" << endl;
    }

    fp = fopen("tmp/schiff.pdb", "w");
    p.save_pdb(fp, &ligand);
    fclose(fp);
    cout << "Wrote output PDB file." << endl;

    return 0;
}