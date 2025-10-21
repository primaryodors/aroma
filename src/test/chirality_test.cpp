
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "../classes/molecule.h"

using namespace std;

int main(int argc, char** argv)
{
    char tstname[1024] = "BZN.sdf";

    if (argc > 1) strcpy(tstname, argv[1]);

    char buffer[65536];
    Molecule m1(tstname);
    cout << "# Created molecule named " << m1.get_name() << ".\n";
    int nloaded;


    FILE* pf = fopen(tstname, "rb");
    if (!pf)
    {
        m1.from_smiles(argv[1]);
    }
    else
    {
        int got = fread(buffer, 1, 65535, pf);
        fclose(pf);
        cout << "Read data.\n";
        nloaded = m1.from_sdf(buffer);
    }

    cout << "# Loaded " << nloaded << " atoms of " << tstname << " into molecule.\n";

    if (m1.is_chiral()) cout << "Molecule is chiral." << endl;
    else cout << "Molecule is NOT chiral." << endl;

    return 0;
}