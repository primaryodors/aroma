
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
    float strain = m.total_bond_strain();
    cout << "Initial bond strain: " << strain << " kJ/mol." << endl;

    Progressbar pgb;
    pgb.set_color(128, 224, 64);

    // m.mangle();
    m.refine_structure(_evolution_default_generations, _default_mutation_rate, _default_population_size,
        nullptr, &pgb);

    strain = m.total_bond_strain();
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