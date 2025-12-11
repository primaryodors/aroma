#include <iostream>
#include <stdio.h>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include "classes/protein.h"
#include "classes/progress.h"

using namespace std;

void show_usage()
{
    cout << "Example usage:" << endl << endl;
    cout << "bin/odorophore sdf/[molecule #1's .sdf file] sdf/[molecule #2's .sdf file]" << endl;
}

int main(int argc, char** argv)
{
    Molecule existing("existing"), added("added");

    int i, n;
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
            filesize = fread(buffer, 1, filesize, fp);         // piece of junk throws a warning if we just eat the return value.
            buffer[filesize] = 0;
            fclose(fp);

            if (!existing.get_atom_count()) existing.from_sdf(buffer);
            else added.from_sdf(buffer);
        }
    }

    if (!existing.get_atom_count() && !added.get_atom_count())
    {
        cerr << "This utility requires two molecules." << endl;
        show_usage();
        return -1;
    }

    existing.movability = MOV_ALL;
    added.movability = MOV_ALL;
    existing.minimize_internal_clashes();
    added.minimize_internal_clashes();

    existing.recenter(Point(0,0,0));
    added.recenter(Point(0,0,0));

    float x, y, z, step=2.0*fiftyseventh;
    Vector ax = Point(1,0,0), ay = Point(0,1,0), az = Point(0,0,1);
    Pose best(&added);
    float bestc = 0;

    Progressbar pgb;
    pgb.maximum = M_PI*2;
    pgb.set_color(192, 64, 80);
    cout << endl;

    for (x=0; x<M_PI*2; x+=step)
    {
        for (y=0; y<M_PI*2; y+=step)
        {
            for (z=0; z<M_PI*2; z+=step)
            {
                // Interaction e = existing.get_intermol_clashes(&added);
                float p = existing.similar_atom_proximity(&added);
                if (p > bestc)
                {
                    best.copy_state(&added);
                    bestc = p;
                }
                
                added.rotate(&az, step, false);
            }

            added.rotate(&ay, step, false);
        }

        added.rotate(&ax, step, false);
        pgb.update(x);
    }
    pgb.erase();

    best.restore_state(&added);

    Atom *a1, *a2;

    // a1 = existing.get_most_bindable()

    n = added.get_atom_count();
    for (i=0; i<n; i++)
    {
        a2 = added.get_atom(i);
        existing.add_existing_atom(a2);
    }

    std::string outfname = "odorophore.pdb";
    fp = fopen(outfname.c_str(), "wb");
    if (!fp)
    {
        cerr << "Failed to open " << outfname << " for reading." << endl;
    }
    existing.save_pdb(fp, 0, true);
    fclose(fp);

    return 0;
}