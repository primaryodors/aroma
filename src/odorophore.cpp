#include <iostream>
#include <stdio.h>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include "classes/scoring.h"

using namespace std;

void show_usage()
{
    cout << "Usage: ummm... lemme get back to you on that." << endl << endl;
    cout << "" << endl;
}

int main(int argc, char** argv)
{
    Molecule existing("existing"), added("added");

    int i;
    FILE* fp;

    for (i=1; i<argc; i++)
    {
        char* ext = strstr(argv[i], ".sdf");
        if (ext)
        {
            fp = fopen(argv[i], "r");
            if (!fp)
            {
                cerr << "Failed to open " << argv[i] << " for reading." << endl;
            }

            fseek(fp, 0, SEEK_END);
            int filesize = ftell(fp);
            fseek(fp, 0, 0);

            char buffer[filesize+16];
            int idgaf = fread(buffer, 1, filesize, fp);         // piece of junk throws a warning if we just eat the return value.
            fclose(fp);

            if (!existing.get_atom_count()) existing.from_sdf(buffer);
            else added.from_sdf(buffer);
        }
    }
}