#include <iostream>
#include <stdio.h>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include "classes/search.h"
#include "classes/reshape.h"
#include "classes/progress.h"

using namespace std;

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        cout << "Usage:" << endl << endl;
        cout << "bin/ic path/to/structure.pdb [strand ID] [save]" << endl;
        cout << endl << "Optional strand ID is a letter between A and Z specifying the target strand of the source PDB file." << endl;
        cout << "This program will optimize the positions of certain atoms. If the optional save param is specified, it will overwrite the original PDB with the optimized structure." << endl;
        return -1;
    }
    char* fname;
    Protein p("testing");
    bool dosave = false, dominc = false, polaronly = false, dohyd = false, totalall = false, dorshps = false;
    ResiduePlaceholder rotres;
    char* rota1 = nullptr;
    char* rota2 = nullptr;
    float rota = 0;
    Reshape rshp;

    int i;
    FILE* fp;
    float threshold = -3;
    for (i=1; i<argc; i++)
    {
        if (!strcmp(argv[i], "save"))
        {
            dosave = true;
        }
        if (!strcmp(argv[i], "minc"))
        {
            dominc = true;
        }
        if (!strcmp(argv[i], "hydro"))
        {
            dohyd = true;
        }
        else if (!strcmp(argv[i], "nooil"))
        {
            polaronly = true;
        }
        else if (!strcmp(argv[i], "tall"))
        {
            totalall = true;
        }
        else if (!strcmp(argv[i], "rot"))
        {
            rotres.set(argv[++i]);
            rota1 = argv[++i];
            rota2 = argv[++i];
            rota = atof(argv[++i]) * fiftyseventh;
        }
        else if (atof(argv[i])) threshold = atof(argv[i]);
        else if (strstr(argv[i], ".pdb"))
        {
            fname = argv[i];
            p.set_name_from_pdb_name(fname);
            fp = fopen(fname, "rb");
            if (!fp)
            {
                cout << "Failed to open " << fname << " for reading." << endl;
                return -1;
            }
            if (i<argc-1 && strlen(argv[i+1]) == 1) p.load_pdb(fp, 0, argv[i+1][0]);
            else p.load_pdb(fp);
            fclose(fp);

            cout << "Loaded " << p.get_seq_length() << " residues." << endl;
        }
        else if (strstr(argv[i], ".rshpm"))
        {
            rshp.load_rshpm_file(argv[i], nullptr);
            dorshps = true;
        }
    }

    if (dohyd)
    {
        int resno, endres = p.get_end_resno();
        cout << "Hydrogenating...";
        for (resno=1; resno<=endres; resno++)
        {
            AminoAcid* res = p.get_residue(resno);
            if (res)
            {
                res->hydrogenate();
            }
            cout << "." << flush;
        }
        cout << endl;
    }

    if (dorshps) rshp.apply(&p, false);

    if (rota1 && rota2 && rota)
    {
        rotres.resolve_resno(&p);
        AminoAcid* aa = p.get_residue(rotres.resno);
        if (aa)
        {
            Atom* a = aa->get_atom(rota1);
            if (a)
            {
                Bond* b = a->get_bond_between(rota2);
                if (b)
                {
                    if (!b->rotate(rota))
                    {
                        cout << "Failed to rotate " << aa->get_name() << ":" << rota1 << "-" << aa->get_name() << ":" << rota2 
                            << " reason " << b->last_fail << "." << endl;
                    }
                }
                else cout << "No bond between " << aa->get_name() << ":" << rota1 << " and " << aa->get_name() << ":" << rota2 << "." << endl;
            }
            else cout << "Atom " << aa->get_name() << ":" << rota1 << " not found." << endl;
        }
        else cout << "Residue " << rotres.resno << " not found." << endl;
    }

    int j, n = p.get_end_resno();
    if (dominc)
    {
        for (i=1; i<=n; i++)
        {
            AminoAcid* aa = p.get_residue(i);
            if (!aa) continue;
            AminoAcid* cl[SPHREACH_MAX+4];
            p.get_residues_can_clash_ligand(cl, aa, aa->get_barycenter(), Point(8,8,8), nullptr);
            Interaction e = aa->get_intermol_binding(cl, false);
            if (e.clash > clash_limit_per_aa) p.minimize_residue_clashes(i);
            cout << "." << flush;
        }
        cout << endl;
    }

    p.optimize_hydrogens();

    float ttl = 0;
    for (i=1; i<=n; i++)
    {
        AminoAcid* aa = p.get_residue(i);
        if (!aa) continue;
        AminoAcid* cl[SPHREACH_MAX+4];
        p.get_residues_can_clash_ligand(cl, aa, aa->get_barycenter(), Point(8,8,8), nullptr);
        BallesterosWeinstein bw1 = p.get_bw_from_resno(i);
        float bestf = 0;
        Pose best(aa);

        for (j=0; cl[j]; j++)
        {
            BallesterosWeinstein bw2 = p.get_bw_from_resno(cl[j]->get_residue_no());
            if (bw1.helix_no && bw2.helix_no && bw1.helix_no >= bw2.helix_no) continue;
            // if (p.get_bw_from_resno(cl[j]->get_residue_no()).helix_no < bw1.helix_no) continue;
            if (fabs(cl[j]->get_residue_no() - i) < 5) continue;

            Atom *a, *b;
            aa->Molecule::mutual_closest_atoms(cl[j], &a, &b);
            if (!a || !b) continue;
            float r = a->distance_to(b);
            if (r > _INTERA_R_CUTOFF) continue;

            bool contact_polar = (fabs(a->is_polar()) >= hydrophilicity_cutoff && fabs(a->get_heavy_atom()->is_polar()) >= hydrophilicity_cutoff)
                && (fabs(b->is_polar()) >= hydrophilicity_cutoff && fabs(b->get_heavy_atom()->is_polar()) >= hydrophilicity_cutoff);

            if (polaronly && !totalall)
            {
                if (!contact_polar) continue;
            }

            Interaction e = aa->get_intermol_binding(cl[j], false);
            if (a->distance_to(b) > (a->vdW_radius+b->vdW_radius)) e.clash = 0;
            float f = e.summed();
            ttl += f;
            if (polaronly)
            {
                if (!contact_polar) continue;
            }

            a = a->get_heavy_atom();
            b = b->get_heavy_atom();
            r = a->distance_to(b);

            if ((threshold < 0 && f <= threshold) || (threshold > 0 && f >= threshold))
            {
                cout << *aa;
                if (bw1.helix_no && bw1.member_no) cout << "(" << bw1 << ")";
                if (a) cout << "." << a->name;
                cout << "-" << *cl[j];
                if (bw2.helix_no && bw2.member_no) cout << "(" << bw2 << ")";
                if (b) cout << "." << b->name;
                cout << ": " << r << " Å; " << f << " kJ/mol." << endl;
            }
            if (f < bestf)
            {
                best.copy_state(aa);
                bestf = f;
            }
        }

        if (bestf) best.restore_state(aa);
    }

    cout << "Total: " << ttl << " kJ/mol." << endl;

    if (dosave)
    {
        // Save optimized PDB structure.
        fp = fopen(fname, "wb");
        if (fp)
        {
            p.save_pdb(fp);
            p.end_pdb(fp);
            fclose(fp);
        }
        else cerr << "Failed to write " << fname << ", check permissions." << endl;
    }
}
