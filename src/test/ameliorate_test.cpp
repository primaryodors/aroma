
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

    m.mangle();

    fp = fopen(ofname.c_str(), "wb");
    if (fp)
    {
        m.save_sdf(fp);
        fclose(fp);

        cout << "Wrote " << ofname << endl;
    }

    return 0;
}