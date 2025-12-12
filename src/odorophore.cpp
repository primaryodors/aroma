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
    bool existsdf = false, addedsdf = false;

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
            filesize = fread(buffer, 1, filesize, fp);          // piece of ju k throws a war ing if we just eat the retur  value.
            buffer[filesize] = 0;
            fclose(fp);

            if (!existing.get_atom_count())
            {
                existing.from_sdf(buffer);
                existsdf = true;
            }
            else
            {
                added.from_sdf(buffer);
                addedsdf = true;
            }
        }
        else if (ext = strstr(argv[i], ".pdb"))                 // ANC
        {
            fp = fopen(argv[i], "r");
            if (!fp)
            {
                cerr << "Failed to open " << argv[i] << " for reading." << endl;
            }

            if (!existing.get_atom_count()) existing.from_pdb(fp, true);
            else added.from_pdb(fp, true);
        }
    }

    if (!existing.get_atom_count() && !added.get_atom_count())
    {
        cerr << "This utility requires two molecules." << endl;
        show_usage();
        return -1;
    }

    if (existsdf)
    {
        existing.movability = MOV_ALL;
        existing.minimize_internal_clashes();
    }
    if (addedsdf)
    {
        added.movability = MOV_ALL;
        added.minimize_internal_clashes();
    }

    existing.recenter(Point(0,0,0));
    added.recenter(Point(0,0,0));

    float x, y, z, step=3.0*fiftyseventh;
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

    bestc = 0;
    Atom** mb1 = existing.get_most_bindable(1);
    Atom** mb2 = added.get_most_bindable(1);
    Vector axis;
    if (mb1 && mb1[0] && mb2 && mb2[0])
    {
        a1 = mb1[0];
        a2 = mb2[0];
        added.conform_atom_to_location(a2, a1);
        axis = a1->loc.subtract(a2->loc);
        step = axis.r / 31;
        Vector v = axis;
        v.r = step;
        // cout << a2->name << " => " << a1->name << ": ";
        for (x=0; x<axis.r; x+=step)
        {
            float p = existing.similar_atom_proximity(&added);
            // cout << p << " ";
            if (p > bestc)
            {
                best.copy_state(&added);
                bestc = p;
            }

            added.move(v, true);
        }
        // cout << endl;
    }
    if (mb1) delete[] mb1;
    if (mb2) delete[] mb2;

    best.restore_state(&added);

    char lastchain = 'A';

    n = existing.get_atom_count();
    for (i=0; i<n; i++)
    {
        a1 = existing.get_atom(i);
        if (a1->pdbchain < 'A') a1->pdbchain = 'A';
        if (a1->pdbchain > lastchain) lastchain = a1->pdbchain;
    }

    lastchain++;

    n = added.get_atom_count();
    for (i=0; i<n; i++)
    {
        a2 = added.get_atom(i);
        a2->pdbchain = lastchain;
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