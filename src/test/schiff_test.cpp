
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
    cout << "Wrote output file." << endl;
    return 0;
}