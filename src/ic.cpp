#include <iostream>
#include <stdio.h>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include "classes/search.h"
#include "classes/reshape.h"
#include "classes/progress.h"
#include "classes/dynamic.h"

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
    bool dosave = false, dominc = false, polaronly = false, dohyd = false, totalall = false,
        dorshps = false, dohg = false, dowind = false, doat = false, doprot = false;
    ResiduePlaceholder rotres;
    char* rota1 = nullptr;
    char* rota2 = nullptr;
    float rota = 0;
    Reshape rshp;
    ICHelixGroup hg;
    DynamicMotion winding(&p);
    ResidueAtomPlaceholder atmov, attgt;
    ResiduePlaceholder protr;

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
        else if (!strcmp(argv[i], "verbose"))
        {
            rshp_verbose = true;
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
        else if (!strcmp(argv[i], "wind"))
        {
            winding.type = dyn_wind;
            winding.start_resno.from_string(argv[++i]);
            winding.end_resno.from_string(argv[++i]);
            winding.bias = atof(argv[++i]);
            dowind = true;
        }
        else if (!strcmp(argv[i], "atomto"))
        {
            atmov.set(argv[++i]);
            attgt.set(argv[++i]);
            doat = true;
        }
        else if (!strcmp(argv[i], "prot"))
        {
            protr.set(argv[++i]);
            doprot = true;
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
        else if (strstr(argv[i], ".ic"))
        {
            cout << "Loading internal contacts file " << argv[i] << "..." << endl;
            hg.load_ic_file(argv[i]);
            dohg = true;
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

    if (doprot)
    {
        protr.resolve_resno(&p);
        if (protr.resno)
        {
            AminoAcid *aa = p.get_residue(protr.resno);
            if (aa)
            {
                if (aa->protonate()) cout << "Protonated " << aa->get_name() << endl;
                else cerr << "Failed to protonate " << aa->get_name() << endl;
            }
            else cerr << "Failed to protonate: specified residue lacking from model." << endl;
        }
        else cerr << "Failed to protonate: specified residue not found in model." << endl;
    }

    if (dowind)
    {
        winding.apply_absolute(1);
        cout << (winding.bias < 0 ? "Unwound " : "Wound ") << winding.start_resno.to_string() << "~" << winding.end_resno.to_string() << endl;
    }

    if (dorshps)
    {
        rshp.apply(&p, false);
        cout << "Applied reshaping file." << endl;
    }

    if (dohg)
    {
        Progressbar pbr;
        pbr.maximum = 0;
        pbr.set_color(192, 24, 156);
        cout << "Optimizing helix groups..." << endl;
        float anomaly = hg.optimize_helices(&p, 20, &pbr);
        cout << "Post-optimization anomaly: " << anomaly << endl;
    }

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
                    if (b->rotate(rota))
                        cout << "Rotated " << aa->get_name() << ":" << rota1 << "-" << rota2 << "." << endl;
                    else cout << "Failed to rotate " << aa->get_name() << ":" << rota1 << "-" << rota2
                            << " reason " << b->last_fail << "." << endl;
                }
                else cout << "No bond between " << aa->get_name() << ":" << rota1 << " and " << aa->get_name() << ":" << rota2 << "." << endl;
            }
            else cout << "Atom " << aa->get_name() << ":" << rota1 << " not found." << endl;
        }
        else cout << "Residue " << rotres.resno << " not found." << endl;
    }

    if (doat)
    {
        atmov.resolve_resno(&p);
        attgt.resolve_resno(&p);

        if (atmov.resno && attgt.resno)
        {
            AminoAcid *aamov = p.get_residue(atmov.resno), *aatgt = p.get_residue(attgt.resno);
            if (aamov && aatgt)
            {
                Atom *amov = aamov->get_atom(atmov.get_aname().c_str()), *atgt = aatgt->get_atom(attgt.get_aname().c_str());
                if (amov && atgt)
                {
                    aamov->conform_atom_to_location(amov, atgt, 50, 3);
                    cout << "Moved " << aamov->get_name() << ":" << amov->name
                        << " to " << amov->distance_to(atgt) << "A from "
                        << aatgt->get_name() << ":" << atgt->name
                        << "." << endl;
                }
                else cerr << "Failed to perform atomto: one or both atoms not present in model." << endl;
            }
            else cerr << "Failed to perform atomto: one or both residues not present in model." << endl;
        }
        else cerr << "Failed to perform atomto: one or both residues not found in protein." << endl;
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
    cout << "Optimized hydrogens." << endl;

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
            p.set_pdb_chain('A');           // Important for fixfail.py.
            p.save_pdb(fp);
            p.end_pdb(fp);
            fclose(fp);
            cout << "Wrote " << fname << endl;
        }
        else cerr << "Failed to write " << fname << ", check permissions." << endl;
    }
}
