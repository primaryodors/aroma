
// Heuristic Odor-Receptor Kinetics

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <regex>
#include "classes/protein.h"

using namespace std;

int main(int argc, char** argv)
{
    Molecule m("jlgsux");
    m.from_smiles("c1ccccc1CCOC(=O)C");

    AtomCollection ac;
    ac.add(m._atoms);
    cout << ac.size() << endl;

    return 0;
}