
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstring>
#include "../classes/molecule.h"

using namespace std;

int main(int argc, char** argv)
{
    Molecule m("Test");
    char buffer[65536];

    if (argc > 1)
    {
        FILE* pf = fopen(argv[1], "rb");
        if (!pf)
        {
            m.from_smiles(argv[1]);
        }
        else
        {
            int got = fread(buffer, 1, 65535, pf);
            fclose(pf);
            m.from_sdf(buffer);
        }
    }
    else
    {
        cout << "Usage:" << endl << endl << "test/solvent_test {path/to/structure.sdf or SMILES string}" << endl << endl;
        return -1;
    }

    m.hydrogenate();
    m.minimize_internal_clashes();

    cout << "Surface area: " << m.get_surface_area() << "Å²." << endl;
    cout << "Polar surface area: " << m.get_surface_area(true) << "Å²." << endl;

    float DeltaGsolv = m.solvent_free_energy();
    cout << DeltaGsolv << " kJ/mol" << endl << (DeltaGsolv * _kcal_per_kJ) << " kcal/mol" << endl;

    return 0;
}