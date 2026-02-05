#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "../classes/molecule.h"

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
    H1->bond_to(N, 1);

    return 0;
}