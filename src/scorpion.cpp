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
    cout << "Usage:" << endl << endl;
    cout << "scorpion [options] path/to/filename.pdb [strand]" << endl;
}

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        show_usage();
        return -1;
    }

    bool show_splash = true;
    char* pdbfn = nullptr;
    char strand = 'A';
    int i;
    for (i=1; i<argc; i++)
    {
        if (!strcmp(argv[i], "--nosplash") || !strcmp(argv[i], "-n"))
        {
            show_splash = false;
        }
        else if (argv[i][0] == '-')
        {
            cout << "Unsupported option " << argv[i] << endl;
            return -1;
        }
        else if (!pdbfn)
        {
            pdbfn = argv[i];
        }
        else if (pdbfn && strlen(argv[i]) == 1)
        {
            strand = argv[i][0] & 0xdf;
        }
        else
        {
            show_usage();
            return -1;
        }
    }

    if (show_splash) cout << "                 _______\n                /\\ _|_/ \\       \n               /__/   \\__\\\n              /__/    (   )\n             |__/      \\ /       @@@                            @\n             |__|       V       @   @\n             |__|       |       @       @@    @@   @ @@  @ @@  @@@   @@   @ @@\n             |___\\               @@@   @  @  @  @  @@  @ @@  @  @   @  @  @@  @\n             |____\\                 @ @     @    @ @     @   @  @  @    @ @   @\n          ,===\\____\\===.        @   @  @  @  @  @  @     @  @   @   @  @  @   @\n         / ,===\\____\\===.        @@@    @@    @@   @     @@@   @@@   @@   @   @\n        / / ,===\\____\\===.                               @\n         / / ,===\\____\\===.                              @\n          / /    __.o..o. __                             @\n        _--------_ `-VV-'   \\_-------_\n       /  ------  \\/     \\_ /  -----  \\\n      /  _,_,_,_,_,\\       /_,_,_,_,_  \\                 model scoring utility\n      |  \\___               \\          |\n       \\__   ---___          `--------'\n          ----_____\\\n" << endl << endl;

    Protein p("TheProtein");
    Molecule m("TheLigand");
    FILE* fp = fopen(pdbfn, "rb");
    if (!fp)
    {
        cout << "Please check input file exists and is readable." << endl;
        return -2;
    }

    p.load_pdb(fp, 0, strand);
    fseek(fp, 0, 0);
    m.from_pdb(fp, true);
    fclose(fp);

    cout << "Ligand net charge: " << m.get_charge() << endl << endl;

    AminoAcid* aa = nullptr;
    if (p.has_bwnos())
    {
        aa = p.get_residue_bw(3, 50);
        if (!aa) aa = p.get_residue_bw(6, 48);
    }
    if (!aa) aa = p.get_residue(100);
    if (!aa) aa = p.get_residue(p.get_end_resno()/2);

    if (aa && aa->get_hydrogen_count() < 3)
    {
        int i, n = p.get_end_resno();
        for (i=0; i<n; i++)
        {
            aa = p.get_residue(i);
            if (aa) aa->hydrogenate();
        }
    }

    m.identify_acidbase();
    m.hydrogenate();

    fp = fopen("tmp/scored_ligand.sdf", "w");
    m.save_sdf(fp);
    fclose(fp);

    DockResult dr(&p, &m, Point(10000,10000,10000));
    dr.include_pdb_data = false;
    dr.display_binding_atoms = true;
    dr.display_clash_atoms = true;
    cout << dr;
}

