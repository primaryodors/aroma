#include <cstring>
#include <iostream>
#include <stdio.h>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sstream>
#include <unistd.h>
#include <fstream>
#include <sys/stat.h>
#include <sys/wait.h>
#include "classes/protein.h"
#include "classes/dynamic.h"

using namespace std;

void usage()
{
    cout << "Usage:" << endl << "bin/poslig [path/to/input.pdb] [path/to/output.pdb] [ligand_atom] [residue_num]:[residue_atom] ..." << endl;
}

int main(int argc, char** argv)
{
    if (argc < 4)
    {
        usage();
        return 1;
    }

    std::string infname = argv[1], outfname = argv[2];

    Protein p("prot");
    Molecule m("lig");

    if (!file_exists(infname.c_str()))
    {
        cerr << "File not found " << infname << endl;
        return -1;
    }

    FILE* fp = fopen(infname.c_str(), "rb");
    if (!fp)
    {
        cerr << "FAILED to open " << infname << " for writing." << endl;
        return -1;
    }

    int resloaded = p.load_pdb(fp);
    cout << "Loaded " << resloaded << " residues." << endl;
    rewind(fp);
    int atsloaded = m.from_pdb(fp, true);
    cout << "Loaded " << atsloaded << " ligand atoms." << endl;

    int i;
    Point alnloc[3];
    for (i=3; i<argc-1; i++)
    {
        Atom *a = m.get_atom(argv[i]);
        if (!a)
        {
            cerr << "No atom named " << argv[i] << " found in ligand." << endl;
        }
        i++;
        char* colon = strchr(argv[i], ':');
        if (!colon)
        {
            usage();
            return 2;
        }
        char *resno = argv[i], *aname = colon+1;
        *colon = 0;
        while (*resno && (*resno < '0' || *resno > '9')) resno++;
        if (!*resno)
        {
            usage();
            return 2;
        }

        AminoAcid *aa = p.get_residue(atoi(resno));
        if (!aa)
        {
            cerr << "Residue " << resno << " not found in protein." << endl;
            return 2;
        }
        Atom *b = aa->get_atom(aname);
        if (!b)
        {
            cerr << "No atom named " << aname << " found in " << aa->get_name() << endl;
            return 2;
        }

        if (!alnloc[0].magnitude())
        {
            alnloc[0] = b->loc;
            Vector mov = alnloc[0].subtract(a->loc);
            mov.r -= 3;
            if (mov.r > 0) m.move(mov);
        }
        else if (!alnloc[1].magnitude())
        {
            alnloc[1] = b->loc;
            Rotation rot = align_points_3d(a->loc, alnloc[1], alnloc[0]);
            m.rotate(rot, alnloc[0]);
        }
        else if (!alnloc[2].magnitude())
        {
            alnloc[2] = b->loc;
            LocatedVector axis = (Vector)(alnloc[1].subtract(alnloc[0]));
            float theta = find_angle_along_vector(alnloc[2], a->loc, alnloc[0], axis);
            axis.origin = alnloc[0];
            Point plus = rotate3D(a->loc, axis.origin, axis, theta),
                  minus = rotate3D(a->loc, axis.origin, axis, -theta);
            if (minus.get_3d_distance(alnloc[2]) < plus.get_3d_distance(alnloc[2])) theta = -theta;
            m.rotate(axis, theta);
        }
        else
        {
            cerr << "WARNING: Maximum number of contact points exceeded; only the first 3 will be used." << endl;
            break;
        }
    }

    fp = fopen(outfname.c_str(), "wb");
    if (!fp)
    {
        cerr << "FAILED to open " << outfname << " for writing." << endl;
        return -1;
    }

    p.save_pdb(fp, &m);
    p.end_pdb(fp);
    fclose(fp);
    cout << "Wrote " << outfname << endl;
}
