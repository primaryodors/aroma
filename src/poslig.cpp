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
    cout << "Usage:" << endl << "bin/poslig [path/to/input.pdb] [path/to/output.pdb] [--receptor distance] [ligand_atom] [residue_num]:[residue_atom] ..." << endl;
}

int main(int argc, char** argv)
{
    if (argc < 4)
    {
        usage();
        return 1;
    }

    std::string infname = argv[1], outfname = argv[2];
    bool turn_receptor=argc>4 && !strcmp(argv[3],"--receptor");
    float contact_distance=turn_receptor ? atof(argv[4]) : 3;
    int first_pair=turn_receptor?5:3;
    if (argc<=first_pair || (argc-first_pair)%2 || contact_distance<1.5 || contact_distance>4)
    { usage(); return 2; }

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

    char input_chain='A', line[1024];
    while (fgets(line,sizeof(line),fp)) if (!strncmp(line,"ATOM  ",6)) { input_chain=line[21]; break; }
    rewind(fp);
    int resloaded = p.load_pdb(fp);
    p.set_pdb_chain(input_chain);
    cout << "Loaded " << resloaded << " residues." << endl;
    rewind(fp);
    int atsloaded = m.from_pdb(fp, true);
    cout << "Loaded " << atsloaded << " ligand atoms." << endl;

    int i;
    Point alnloc[3];
    for (i=first_pair; i<argc-1; i++)
    {
        Atom *a = m.get_atom(argv[i]);
        if (!a)
        {
            cerr << "No atom named " << argv[i] << " found in ligand." << endl;
            return 2;
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

        if (turn_receptor)
        {
            aa->movability=MOV_FLEXONLY;
            aa->conform_atom_to_location(b,a,10,contact_distance,true);
            if (fabs(b->distance_to(a)-contact_distance)>.2)
            { cerr << "Side chain cannot reach the requested contact." << endl; return 2; }
            cout << "Inward contact " << aa->get_name() << ":" << b->name << " to " << a->name
                 << " at " << b->distance_to(a) << " A; surrounding clashes are left for repair." << endl;
            cout << "Contact restraint " << a->name << " " << aa->get_residue_no() << ":" << b->name
                 << " " << contact_distance << endl;
            Atom* h=a->is_bonded_to("H");
            float radius=InteratomicForce::hydrogen_bond_distance(h,b);
            bool ligand_donor=radius>0;
            if (!ligand_donor) { h=b->is_bonded_to("H"); radius=InteratomicForce::hydrogen_bond_distance(h,a); }
            if (radius>0)
            {
                Bond* axis=nullptr;
                for (int j=0;j<a->get_bonded_atoms_count();j++)
                    if (a->get_bond_by_idx(j)->atom2->Z>1) axis=a->get_bond_by_idx(j)->get_reversed();
                Atom* best_h=h;
                if (!ligand_donor || (axis && axis->can_rotate))
                {
                    Pose origin(&m), best(&m), receptor_origin(aa), receptor_best(aa); float error=INFINITY;
                    std::vector<Bond*> axes;
                    Bond** bonds=aa->get_rotatable_bonds();
                    if (bonds) for (int j=0;bonds[j] && axes.size()<2;j++)
                        if (bonds[j]->can_rotate && !bonds[j]->is_ring_bond() && !bonds[j]->atom2->is_backbone
                            && (!bonds[j]->atom1->is_backbone || !strcmp(bonds[j]->atom1->name,"CA")))
                            axes.push_back(bonds[j]);
                    for (int x=0;x<(axes.empty()?1:24);x++) for (int y=0;y<(axes.size()<2?1:24);y++)
                    {
                        receptor_origin.restore_state(aa);
                        if (axes.size()) axes[0]->rotate(x*M_PI/12,false,true);
                        if (axes.size()>1) axes[1]->rotate(y*M_PI/12,false,true);
                        for (int j=0;j<(ligand_donor?72:1);j++)
                        {
                            origin.restore_state(&m); if (ligand_donor) axis->rotate(j*M_PI/36,false,true);
                            for (int k=0;k<(ligand_donor?1:b->get_bonded_atoms_count());k++)
                            {
                                Atom* hydrogen=ligand_donor?h:b->get_bond_by_idx(k)->atom2;
                                if (hydrogen->Z!=1) continue;
                                Atom* acceptor=ligand_donor?b:a;
                                Atom* donor=ligand_donor?a:b;
                                float trial=fabs(hydrogen->distance_to(acceptor)-radius)+fabs(a->distance_to(b)-contact_distance);
                                if (find_3d_angle(donor->loc,acceptor->loc,hydrogen->loc)*fiftyseven<135) trial+=1;
                                if (trial<error) {error=trial; best.copy_state(&m); receptor_best.copy_state(aa); best_h=hydrogen;}
                            }
                        }
                    }
                    best.restore_state(&m);
                    receptor_best.restore_state(aa);
                }
                cout << (ligand_donor?"Contact restraint ":"Contact donor ") << (ligand_donor?best_h->name:a->name) << " " << aa->get_residue_no() << ":" << (ligand_donor?b->name:best_h->name)
                     << " " << radius << endl;
            }
        }
        else if (!alnloc[0].magnitude())
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
