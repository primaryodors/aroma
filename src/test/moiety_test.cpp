
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <sys/stat.h>
#include "../classes/moiety.h"

using namespace std;

int main(int argc, char** argv)
{
    Molecule m("Testing");

    if (argc > 1)
    {
        struct stat buffer;   
        if (stat(argv[1], &buffer) == 0)
        {
            FILE* fp = fopen(argv[1], "rb");
            if (!fp) { cerr << "Could not read ligand file." << endl; return 1; }

            fseek(fp, 0L, SEEK_END);
            int fsize = ftell(fp);
            rewind(fp);

            if (fsize<=0 || fsize>65535) { fclose(fp); cerr << "Invalid ligand file size." << endl; return 1; }
            std::string contents(fsize,'\0');
            size_t got = fread(&contents[0], 1, contents.size(), fp);
            fclose(fp);
            if (got!=contents.size() || !m.from_sdf(contents.c_str()))
            { cerr << "Invalid or incomplete ligand SDF." << endl; return 1; }
        }
        else m.from_smiles(argv[1]);
    }
    else m.from_smiles("CCO");

    Moiety y;
    if (argc > 2) y.pattern = argv[2];
    else y.pattern = "COH";

    Atom* matches[256];
    int result = y.contained_by(&m, matches);

    cout << "Pattern " << y.pattern << " occurs " << result << " times in molecule." << endl;
    int i;
    for (i=0; matches[i]; i++)
    {
        cout << *(matches[i]) << " ";
    }
    cout << endl;
}
