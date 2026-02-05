#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "../classes/protein.h"

using namespace std;

Molecule* mols[3];

int main(int argc, char** argv)
{
    mols[0] = new Molecule("aldehyde", "c1ccccc1C=O");
    mols[1] = new Molecule("amine", "CCCC[NH3+]");
    mols[2] = nullptr;

    Atom *O = mols[0]->get_atom("O8");
    Atom *N = mols[1]->get_atom("N5");
    Atom *H1 = N->is_bonded_to("H");
    H1->unbond_all();
    Atom *H2 = N->is_bonded_to("H");
    O->increment_charge(1);
    Vector v = O->get_next_free_geometry(1);
    H1->bond_to(O, 1);
    H1->move(O->loc.add(v));

    FILE* fp = fopen("oh_schiff.sdf", "w");
    mols[0]->save_sdf(fp, &mols[1]);
    fclose(fp);

    return 0;
}