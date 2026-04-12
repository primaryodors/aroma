
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstring>
#include "classes/molecule.h"

using namespace std;

int main(int argc, char** argv)
{
    Molecule mol("molecule");

    std::string fname = "sdf/geraniol.sdf";
    if (argc > 1)
    {
        if (file_exists(argv[1]))
            fname = argv[1];
    }

    char buffer[65536];
    FILE* fp = fopen(fname.c_str(), "rb");
    if (!fp)
    {
        cerr << "FAILED to open " << fname << " for reading." << endl;
        return -1;
    }
    fread(buffer, sizeof(char), 65530, fp);
    fclose(fp);
    mol.from_sdf(buffer);

    mol.optimize();

    fp = fopen(fname.c_str(), "wb");
    if (!fp)
    {
        cerr << "FAILED to open " << fname << " for writing." << endl;
        return -1;
    }

    mol.save_sdf(fp);
    fclose(fp);
    cout << "Wrote " << fname << endl;
}