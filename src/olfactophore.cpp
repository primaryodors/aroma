#include <iostream>
#include <stdio.h>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include "classes/search.h"
#include "classes/reshape.h"
#include "classes/progress.h"

#define phore_polar_maxr 3.5
#define phore_aliphatic_maxr 1.0
#define phore_pi_maxr 2.5

#define feature_type_hbacc 0x04b0
#define feature_type_hbdon 0x04b1
#define feature_type_sulph 0x05b0
#define feature_type_aliph 0x0a20
#define feature_type_pi    0x0b11
#define feature_type_excl  0xf000

using namespace std;

void show_usage()
{
    cout << "Example usage:" << endl << endl;
    cout << "bin/olfactophore sdf/[molecule #1's .sdf file] sdf/[molecule #2's .sdf file]" << endl;
}

Pose phore_tumble_3d(Molecule& existing, Molecule& added)
{
    float x, y, z, step=3.6*fiftyseventh;
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

    return best;
}

Pose phore_bindable_adjustment(Molecule& existing, Molecule& added)
{
    Pose best(&added);
    Atom *a1, *a2, *a3;
    float bestc = 0;
    Atom** mb1 = existing.get_most_bindable(1);
    Atom** mb2 = added.get_most_bindable(1);
    Vector axis;
    if (mb1 && mb1[0] && mb2 && mb2[0])
    {
        a1 = mb1[0];
        a2 = mb2[0];
        Pose putitback(&added);
        float p = existing.similar_atom_proximity(&added);
        added.conform_atom_to_location(a2, a1);
        if (existing.similar_atom_proximity(&added) < 0.75*p) putitback.restore_state(&added);
        axis = a1->loc.subtract(a2->loc);
        float x, step = axis.r / 31;
        Vector v = axis;
        v.r = step;
        // cout << a2->name << " => " << a1->name << ": ";
        for (x=0; x<axis.r; x+=step)
        {
            p = existing.similar_atom_proximity(&added);
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
    return best;
}

int find_features(Molecule& existing, Molecule& added, Point* feature, float* featurer, int* feattype)
{
    int i, j, l, m, n, nfeature=0;
    Atom *a1, *a2, *a3;

    n = existing.get_atom_count();
    bool dirty[n+4];
    for (i=0; i<n; i++) dirty[i] = false;
    for (i=0; i<n; i++)
    {
        if (dirty[i]) continue;
        a1 = existing.get_atom(i);
        int fam1 = a1->get_family();
        float pol1 = a1->is_polar();
        if (fabs(pol1) < hydrophilicity_cutoff) pol1 = 0;
        if (fam1 == TETREL) pol1 = 0;
        if (fam1 == TETREL && a1->is_bonded_to(CHALCOGEN, 2)) continue;
        bool pi1 = a1->is_pi();
        if (fam1 == PNICTOGEN && (a1->get_bonded_atoms_count() > 2 || a1->is_bonded_to("H"))) pol1 = 0;
        for (j=i+1; j<n; j++)
        {
            if (dirty[j]) continue;
            a2 = existing.get_atom(j);
            int fam2 = a2->get_family();
            if (a1->pdbchain == a2->pdbchain && a2->is_bonded_to(a1)) continue;
            float pol2 = a2->is_polar();
            if (fabs(pol2) < hydrophilicity_cutoff) pol2 = 0;
            if (fam2 == TETREL) pol2 = 0;
            if (fam2 == TETREL && a2->is_bonded_to(CHALCOGEN, 2)) continue;
            bool pi2 = a2->is_pi();
            if (fam2 == PNICTOGEN && (a2->get_bonded_atoms_count() > 2 || a2->is_bonded_to("H"))) pol2 = 0;

            float rij = a1->distance_to(a2);
            m = -1;
            if ((pol1 && sgn(pol2) == sgn(pol1) && rij <= phore_polar_maxr)
                ||
                (!pol1 && !pol2 && pi1 && pi2 && rij <= phore_pi_maxr)
                ||
                (!pol1 && !pol2 && rij <= phore_aliphatic_maxr)
                )
            {
                float pol3;
                int fam3;
                bool pi3;
                float best3 = 0;
                for (l=j+1; l<n; l++)
                {
                    if (dirty[l]) continue;
                    a3 = existing.get_atom(l);
                    fam3 = a3->get_family();
                    if (a3->pdbchain == a1->pdbchain && a3->is_bonded_to(a1)) continue;
                    if (a3->pdbchain == a2->pdbchain && a3->is_bonded_to(a2)) continue;
                    pol3 = a3->is_polar();
                    if (fabs(pol3) < hydrophilicity_cutoff) pol3 = 0;
                    if (fam3 == TETREL) pol3 = 0;
                    if (fam3 == TETREL && a3->is_bonded_to(CHALCOGEN, 2)) continue;
                    pi3 = a3->is_pi();
                    if (fam3 == PNICTOGEN && (a3->get_bonded_atoms_count() > 2 || a3->is_bonded_to("H"))) pol3 = 0;

                    float ril = a1->distance_to(a3);
                    float rjl = a2->distance_to(a3);

                    if ((fam1 == CHALCOGEN && fam2 == CHALCOGEN && fam3 == CHALCOGEN && ril <= phore_polar_maxr && rjl <= phore_polar_maxr)
                        ||
                        (pol1 && sgn(pol2) == sgn(pol1) && sgn(pol2) == sgn(pol3) && ril <= phore_polar_maxr && rjl <= phore_polar_maxr)
                        ||
                        (!pol1 && !pol2 && !pol3 && pi1 && pi2 && pi3 && ril <= phore_pi_maxr && rjl <= phore_pi_maxr)
                        ||
                        (!pol1 && !pol2 && !pol3 && ril <= phore_aliphatic_maxr && rjl <= phore_aliphatic_maxr)
                        )
                    {
                        float f = 0;
                        if (pol1 && pol2 && pol3) f += 35;
                        if (pi1 && pi2 && pi3) f += 5;
                        if (!pol1 && !pol2 && !pol3) f += 15;
                        if (f > best3)
                        {
                            m = l;
                            best3 = f;
                        }
                    }
                }

                if (m>=0)
                {
                    a3 = existing.get_atom(m);
                    pol3 = a3->is_polar();
                    fam3 = a3->get_family();
                    if (fabs(pol3) < hydrophilicity_cutoff) pol3 = 0;
                    if (fam3 == TETREL) pol3 = 0;
                    pi3 = a3->is_pi();
                    if (fam3 == PNICTOGEN && (a3->get_bonded_atoms_count() > 2 || a3->is_bonded_to("H"))) pol3 = 0;
                }

                Point featcen = a1->loc.add(a2->loc);
                if (m>=0) featcen = featcen.add(a3->loc);
                featcen.multiply(1.0 / ((m>=0)?3:2));
                float featr = a1->loc.get_3d_distance(featcen) + a1->vdW_radius;
                float r = a2->loc.get_3d_distance(featcen) + a2->vdW_radius;
                if (r > featr) featr = r;
                if (m>=0)
                {
                    r = a3->loc.get_3d_distance(featcen) + a3->vdW_radius;
                    if (r > featr) featr = r;
                }

                feature[nfeature] = featcen;
                featurer[nfeature] = featr;

                #if 0
                cout << a1->pdbchain << ":" << a1->name << " "
                    << a2->pdbchain << ":" << a2->name << " "
                    << (m<0 ? ' ' : a3->pdbchain) << ":" << (m<0 ? "" : a3->name) << endl;
                #endif
                if (fam1 == CHALCOGEN && fam2 == CHALCOGEN && (m<0 || fam3 == CHALCOGEN) && a1->Z > 10 && a2->Z > 10 && (m<0 || a3->Z > 10))
                    feattype[nfeature] = feature_type_sulph;
                else if (pol1 < 0 && pol2 < 0 && (m<0 || pol3 < 0))
                    feattype[nfeature] = (a1->Z < 10 && a2->Z < 10 && (m<0 || a3->Z < 10)) ? feature_type_hbacc : feature_type_sulph;
                else if (pol1 > 0 && pol2 > 0 && (m<0 || pol3 > 0))
                    feattype[nfeature] = (a1->get_heavy_atom()->Z < 10 && a2->get_heavy_atom()->Z < 10 
                        && (m<0 || a3->get_heavy_atom()->Z < 10)) ? feature_type_hbdon : feature_type_sulph;
                else if (m>=0 && a1->Z > 1 && a2->Z > 1 && a3->Z > 1 && !pol1 && !pol2 && !pol3 && pi1 && pi2 && pi3)
                    feattype[nfeature] = feature_type_pi;
                else if (m>=0 && a1->Z > 1 && a2->Z > 1 && a3->Z > 1 && !pol1 && !pol2 && !pol3)
                    feattype[nfeature] = feature_type_aliph;
                else continue;

                // cout << hex << feattype[nfeature] << dec << endl;

                nfeature++;

                dirty[i] = true;
                dirty[j] = true;
                if (m>=0) dirty[m] = true;
            }
            if (dirty[j]) continue;
        }
        if (dirty[i]) continue;
    }

    return nfeature;
}

int main(int argc, char** argv)
{
    Protein prot;
    Molecule existing("existing"), added("added");

    int i, j, l, m, n;
    FILE* fp;
    bool existsdf = false, addedsdf = false;
    Point feature[256];
    float featurer[256];
    int feattype[256];
    int nfeature = 0;
    char *infname1 = nullptr, *infname2 = nullptr;
    std::string outfname = "olfactophore.pdb";

    for (i=1; i<argc; i++)
    {
        if (!strcmp(argv[i], "-o"))
        {
            outfname = argv[++i];
            continue;
        }

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
                infname1 = argv[i];
                existing.from_sdf(buffer);
                existsdf = true;
            }
            else
            {
                infname2 = argv[i];
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

            if (!existing.get_atom_count())
            {
                if (!prot.get_seq_length())
                {
                    prot.load_pdb(fp);
                    if (prot.get_seq_length())
                    {
                        fclose(fp);
                        continue;
                    }
                }

                fseek(fp, 0, 0);
                existing.from_pdb(fp, true);
                infname1 = argv[i];
            }
            else
            {
                added.from_pdb(fp, true);
                infname2 = argv[i];
            }
            fclose(fp);
        }
    }

    if (infname1 && infname2)
    {
        if (!strcmp(outfname.c_str(), "1")) outfname = infname1;
        if (!strcmp(outfname.c_str(), "2")) outfname = infname2;

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

        Pose best;
        best = phore_tumble_3d(existing, added);
        best = phore_bindable_adjustment(existing, added);

        char lastchain = 'A';

        Atom *a1, *a2;
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

        nfeature = find_features(existing, added, feature, featurer, feattype);
        // cout << "Found " << nfeature << " features." << endl;

        fp = fopen(outfname.c_str(), "wb");
        if (!fp)
        {
            cerr << "Failed to open " << outfname << " for reading." << endl;
        }

        for (i=0; i<nfeature; i++)
        {
            fprintf(fp, "REMARK 821 0 %.3f %.3f %.3f %.3f ...%c%c%c%c \n",
                feature[i].x, feature[i].y, feature[i].z, featurer[i],
                feattype[i] == feature_type_hbacc ? 'A' : (feattype[i] == feature_type_hbdon ? 'D' : '.'),
                feattype[i] == feature_type_sulph ? 'S' : '.',
                feattype[i] == feature_type_pi ? 'P' : '.',
                feattype[i] == feature_type_excl ? 'X' : '.'
            );
        }

        existing.save_pdb(fp, 0, true);
        fclose(fp);
        cout << "Wrote " << outfname << endl;
    }
    else if (prot.get_seq_length() && infname1)
    {
        Progressbar pgb;
        pgb.set_color(224, 32, 144);
        cout << "Performing pocket search..." << endl << flush;
        Point sz = prot.get_region_bounds(1, 99999);
        sz.multiply(0.666);
        Point cen = prot.find_loneliest_point(Point(0,5,0), sz);                // Do this BEFORE processing rshpms and deleting residues!
        existing.recenter(cen);

        Reshape rshp;
        if (prot.get_residue_bw(7, 50))
        {
            rshp.load_rshpm_file(rshp_GPCR, &existing);
            cout << "Loaded reshape file." << endl << flush;
            rshp.apply(&prot, false);
            cout << "Applied reshape motions." << endl << flush;
        }

        AminoAcid* bsr[SPHREACH_MAX+4];
        int nbsr = prot.get_residues_can_clash_ligand(bsr, &existing, cen, Point(7,7,7));
        #if 1
        // get barycenter of hydrophilicity
        Point hbary(0,0,0);
        j=0;
        for (i=0; i<nbsr; i++)
        {
            if (fabs(bsr[i]->hydrophilicity()) >= hydrophilicity_cutoff)
            {
                Atom* a = bsr[i]->get_reach_atom();
                if (a)
                {
                    Atom* la = existing.get_nearest_atom(a->loc);
                    if (!la) continue;
                    if (la->distance_to(a) > _INTERA_R_CUTOFF) continue;
                    hbary = hbary.add(a->loc);
                    j++;
                }
            }
        }
        if (j)
        {
            hbary.multiply(1.0/j);
            Point hnew(0,0,0);
            j=0;
            for (i=0; i<nbsr; i++)
            {
                if (fabs(bsr[i]->hydrophilicity()) >= hydrophilicity_cutoff)
                {
                    Atom* a = bsr[i]->get_reach_atom();
                    if (a)
                    {
                        float r = a->loc.get_3d_distance(hbary);
                        if (r > _INTERA_R_CUTOFF+1) continue;
                        hnew = hnew.add(a->loc);
                        j++;
                        cout << bsr[i]->get_name() << ":" << a->name << " is " << r << "A from wet barycenter." << endl;
                    }
                }
            }
            if (j)
            {
                hnew.multiply(1.0 / j);
                cen = hnew;
                Atom* a = prot.get_nearest_atom(hnew);
                if (a) cout << "Wet spot is " << a->loc.get_3d_distance(hnew) << " A from "
                    << a->residue << ":" << a->name << endl;
            }
        }

        Point lhc = existing.polar_barycenter();
        if (j)
        {
            Rotation rot = align_points_3d(lhc, cen, existing.get_barycenter());
            existing.rotate(rot);
            lhc = existing.polar_barycenter();
        }
        else prot.tumble_ligand_inside_pocket(&existing, cen, 1, &pgb, &lhc);
        #else
        existing.recenter(cen);
        LigandTarget lt[256];
        int nlt = Search::identify_ligand_pairing_targets(&existing, lt, 256);
        BestBindingResult bbr;
        pgb.set_color(128, 64, 255);
        cout << endl;
        Search::pair_targets(&existing, lt, bsr, cen, &bbr, nullptr, true, &pgb);
        Search::align_targets(&existing, cen, &bbr, 1);
        #endif
        if (prot.get_residue_bw(7, 50))
        {
            rshp.apply(&prot, true);
            cout << "Applied activation motions." << endl << flush;
        }

        fp = fopen(outfname.c_str(), "wb");
        if (!fp)
        {
            cerr << "Failed to open " << outfname << " for reading." << endl;
        }

        prot.set_pdb_chain('A');
        prot.save_pdb(fp, &existing);
        fclose(fp);
        cout << "Wrote " << outfname << endl;
    }
    else
    {
        cerr << "This utility requires two molecules." << endl;
        show_usage();
        return -1;
    }


    return 0;
}