#include <cstring>
#include <iostream>
#include <stdio.h>
#include <string>
#include <math.h>
#include "../classes/molecule.h"

using namespace std;

void set_color(float r, float g, float b)
{
    int lr = max(0, min(255, (int)(r*255)));
    int lg = max(0, min(255, (int)(g*255)));
    int lb = max(0, min(255, (int)(b*255)));

    cout << "\x1b[48;2;" << lr << ";" << lg << ";" << lb << "m";
}

void clear_color()
{
    cout << "\x1b[0m";
}

int main(int argc, char** argv)
{
    Molecule m("Test"), p("Probe");
    cout << "Created empty molecule named " << m.get_name() << ".\n";

    Atom* anisoa;
    bool colors = false;
    bool hfwid = false;

    if (argc > 1 && argv[1][0] != '-')
    {
        m.from_smiles(argv[1]);
    }
    else
    {
        Atom* C1 = m.add_atom("C", "C1", 0, 0);
        cout << "Added a carbon atom. Its location is " << C1->loc.printable() << "." << endl;

        Atom* C2 = m.add_atom("C", "C2", C1, 1);
        cout << "Added another carbon atom. Its location is " << C2->loc.printable() << "." << endl;

        Atom* O3 = m.add_atom("O", "O3", C2, 2);
        cout << "Added an oxygen atom. Its location is " << O3->loc.printable() << "." << endl;

        Atom* C4 = m.add_atom("C", "C4", C2, 1);
        cout << "Added another carbon atom. Its location is " << C4->loc.printable() << "." << endl;

        anisoa = O3;
    }

    m.hydrogenate();
    cout << "Hydrogenated." << endl;

    if (argc > 2 && argv[2][0] != '-')
    {
        anisoa = m.get_atom(argv[2]);
        if (!anisoa)
        {
            cout << argv[2] << " not found in molecule. Exiting." << endl;
            return -1;
        }
    }

    int i;
    for (i=0; i<argc; i++)
    {
        if (!strcmp(argv[i], "--colors")) colors = true;
        else if (!strcmp(argv[i], "--asciiart")) colors = false;
        else if (!strcmp(argv[i], "--hfwid")) hfwid = true;
    }

    const int search_size=22;
    float ar = 2.1;		// Aspect ratio.
    int x, y;
    const char* asciiart = " .':+=inm@";
    int asciilen = strlen(asciiart);

    if (hfwid) ar /= 2;

    Atom* probe;
    if (argc > 4 && argv[3][0] != '-')
    {
        p.from_smiles(argv[3]);
        probe = p.get_atom(argv[4]);
    }
    else
    {
        p.from_smiles("O");
        probe = p.get_atom(1);
    }

    Atom* oxy = probe->get_bonded_atoms_count() ? probe->get_bond_by_idx(0)->atom2 : nullptr;

    InteratomicForce* ifs[32];
    only_closest_conj_charge = false;
    InteratomicForce::fetch_applicable(probe, anisoa, ifs);
    if (!ifs || !ifs[0])
    {
        cout << "No forces to measure; check bindings.dat." << endl;
        return -1;
    }

    InteratomicForce* hb=0;
    float strongest = 0;
    for (x=0; ifs[x]; x++)
    {
        cout << ifs[x]->get_type() << endl;
        if (ifs[x]->get_kJmol() > strongest)
        {
            hb = ifs[x];
            strongest = ifs[x]->get_kJmol();
        }
    }

    if (!hb)
    {
        cout << "No hbond force to measure; check bindings.dat." << endl;
        return -1;
    }

    Point paim(0,10000000,0);

    Bond* bb[16];
    anisoa->fetch_bonds(bb);
    Point aloc = anisoa->loc;
    Point bloc;
    Point bblocs[16];
    if (bb)
    {
        for (i=0; bb[i]; i++)
        {
            if (!bb[i]->atom2) break;
            Point pt = bb[i]->atom2->loc;
            pt = pt.subtract(aloc);
            pt.scale(1);
            pt = pt.add(aloc);
            pt.weight = 1;
            bblocs[i] = pt;
            // cout << bb[i]->atom2->name << endl;
        }

        // cout << i << " points for average." << endl;
        if (i) bloc = average_of_points(bblocs, i);
        // mp.add_atom("He", "He1", &bloc, NULL, 0);
    }

    paim = paim.add(bloc);

    // cout << "Align " << aloc << " to " << paim << " centered at " << bloc << endl;
    Rotation rot = align_points_3d(&aloc, &paim, &bloc);

    // cout << rot << endl;

    m.rotate(&rot.v, rot.a);

    aloc = anisoa->loc;
    minimum_searching_aniso = 0;

    float best_energy = 0;
    for (y=-search_size; y<=search_size; y++)
    {
        for (x=-search_size*ar; x<=search_size*ar; x++)
        {
            float lx = (float)x/ar;
            float r = sqrt(lx*lx + y*y);
            if (r > search_size) cout << " ";
            else
            {
                float phi = find_angle(lx, y);
                float theta = (1.0-r/(search_size))*(M_PI/2);

                Vector v(hb->get_distance(), theta, phi);
                Point loc = anisoa->loc.add(&v);

                Vector motion = loc.subtract(probe->loc);
                p.move(motion);
                
                if (oxy)
                {
                    loc = loc.add(v);
                    Rotation rot = align_points_3d(oxy->loc, loc, probe->loc);
                    p.rotate(rot, probe->loc);
                }

                if (probe->is_pi())
                {
                    Point norm = probe->loc.add(probe->get_pi_normal(&(anisoa->loc)));
                    Rotation rot = align_points_3d(norm, anisoa->loc, probe->loc);
                    // if (!x && !y) cout << endl << (rot.a*fiftyseven) << " degrees" << endl;
                    p.rotate(rot, probe->loc);
                }

                probe->clear_geometry_cache();
                oxy->clear_geometry_cache();

                if (!x && !y)
                {
                    p.been_flexed = false;              // put debug breakpoints here
                }
                float tb = InteratomicForce::total_binding(anisoa, probe).summed();

                if (tb < best_energy) best_energy = tb;

                tb /= hb->get_kJmol() * -1.5;
                if (tb<0) tb=0;

                if (!colors)
                {
                    tb *= asciilen;
                    if (tb > asciilen-1) tb = asciilen-1;

                    cout << asciiart[(int)tb];
                }
                else
                {
                    set_color(tb*tb,tb,sqrt(tb));
                    cout << " ";
                    clear_color();
                }

                if (!x && !y)
                {
                    int anisg = anisoa->get_geometry();
                    Vector* anisgeo = anisoa->get_geometry_aligned_to_bonds();
                    Molecule mptemp("Very temporary");

                    int n = p.get_atom_count();
                    for (i=0; i<n; i++) mptemp.add_existing_atom(p.get_atom(i));

                    if (anisgeo)
                        for (i=0; i<anisg; i++)
                        {
                            Point pt(&anisgeo[i]);
                            pt.scale(0.4);
                            pt = pt.add(aloc);
                            char buffer[16];
                            sprintf(buffer, "He%d", i);
                            mptemp.add_atom("He", buffer, &pt, NULL, 0);
                        }

                    FILE* pf = fopen("aniso.sdf", "wb");
                    Molecule* ligands[3];
                    ligands[0] = &mptemp;
                    ligands[1] = NULL;
                    m.save_sdf(pf, ligands);
                    fclose(pf);
                }
            }
        }
        cout << endl;
    }

    cout << "Binding energy: " << best_energy << " out of an optimal " << -hb->get_kJmol() << " " << hb->get_type() << " kJ/mol." << endl;

}









