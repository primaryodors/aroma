
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstring>
#include "classes/protein.h"

using namespace std;

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        cout << "No input file." << endl;
        return -1;
    }

    char* inpfile = NULL;
    char strand = 'A';
    int i, n;

    for (i=1; i<argc; i++)
    {
        if 		(!strcmp(argv[i], "-p")) inpfile = argv[++i];
        if 		(!strcmp(argv[i], "-s")) strand = argv[++i][0];
        else if (file_exists(argv[i])) inpfile = argv[i];
    }

    if (!inpfile)
    {
        cout << "No input file. Please specify a PDB file to extract the protein sequence from." << endl;
        cout << endl;
        return -1;
    }

    Protein p(inpfile);
    FILE* pf = fopen(inpfile, "r");
    if (!pf)
    {
        cout << "Error trying to read " << inpfile << endl;
        return 0xbadf12e;
    }
    p.load_pdb(pf, 0, strand);
    fclose(pf);
    p.set_name_from_pdb_name(inpfile);
    n = p.get_end_resno();

    bool first = true;
    for (i=1; i<=n; i++)
    {
        AminoAcid* aa = p.get_residue(i);
        if (!aa) continue;
        if (first)
        {
            cout << i << " - " << n << endl;
            first = false;
        }
        cout << aa->get_letter();
    }
    cout << endl;

    return 0;
}
