
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstring>

#include "../classes/molecule.h"

using namespace std;

int main(int argc, char** argv)
{
    FILE* fp;
    char buffer[65536];
    Molecule m("TheMolecule");
    std:;string ofname = "test.sdf";

    if (argc > 1)
    {
        if (file_exists(argv[1]))
        {
            fp = fopen(argv[1], "rb");
            int fku = fread(buffer, sizeof(char), 65530, fp);
            fclose(fp);
            m.from_sdf(buffer);
        }
        else m.from_smiles(argv[1]);
    }
    else return -1;

    int i, n = m.get_atom_count();
    float strain = 0;
    for (i=0; i<n; i++)
    {
        Atom* a = m.get_atom(i);
        if (!a) continue;
        strain += m.get_atom_bond_length_anomaly(a);
        strain += m.get_atom_bond_angle_anomaly(a);
    }
    cout << "Initial bond strain: " << strain << " kJ/mol." << endl;

    m.mangle();
    m.refine_structure();

    strain = 0;
    for (i=0; i<n; i++)
    {
        Atom* a = m.get_atom(i);
        if (!a) continue;
        strain += m.get_atom_bond_length_anomaly(a);
        strain += m.get_atom_bond_angle_anomaly(a);
    }
    cout << "Post-refine bond strain: " << strain << " kJ/mol." << endl;

    fp = fopen(ofname.c_str(), "wb");
    if (fp)
    {
        m.save_sdf(fp);
        fclose(fp);

        cout << "Wrote " << ofname << endl;
    }

    return 0;
}