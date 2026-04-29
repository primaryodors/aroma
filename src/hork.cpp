
// Heuristic Odor-Receptor Kinetics

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <regex>
#include "classes/protein.h"
#include "classes/kinetic.h"

using namespace std;

int main(int argc, char** argv)
{
    FILE* fp;
    Molecule m("jlgsux");
    m.from_smiles("c1ccccc1CCOC(=O)C");

    MolecularKinetics mk;
    mk.atoms.add(m._atoms);
    mk.set_Boltzmann_momenta();
    mk.dump();
    fp = fopen("tmp/frame1.sdf", "w");
    m.save_sdf(fp);
    fclose(fp);

    mk.advance_clock(1);
    mk.dump();
    fp = fopen("tmp/frame2.sdf", "w");
    m.save_sdf(fp);
    fclose(fp);

    mk.advance_clock(1);
    mk.dump();
    fp = fopen("tmp/frame3.sdf", "w");
    m.save_sdf(fp);
    fclose(fp);

    return 0;
}