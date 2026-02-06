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

    mols[0]->create_Schiff_base(mols[1]);

    FILE* fp = fopen("oh_schiff.sdf", "w");
    mols[0]->save_sdf(fp, &mols[1]);
    fclose(fp);

    return 0;
}