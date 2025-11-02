
#include <string.h>
#include "../classes/search.h"
#include "../classes/cavity.h"

int main(int argc, char** argv)
{
    Protein p("Test Receptor");
    Molecule m("Test Ligand");
    std::vector<BallesterosWeinstein> pocketcen_res;
    std::vector<std::string> atoms_of_interest;
    bool priorities[256];
    bool fullscan = false;
    int maxbbr = 10;
    char* cvtyfname = nullptr;
    Cavity cvtys[256];
    int ncvtys = 0;
    int showtopn = 0;

    bool save_tmp_pdbs = false;

    int i, j;
    for (i=0; i<256; i++) priorities[i] = false;

    FILE* fp;
    char buffer[49152];
    int got;
    for (i=1; i<argc; i++)
    {
        char* dot;
        dot = strrchr(argv[i], '.');
        if (dot && !strcasecmp(dot, ".pdb"))
        {
            fp = fopen(argv[i], "rb");
            if (!fp)
            {
                cout << "Failed to open " << argv[i] << " for reading." << endl;
                return -1;
            }
            p.load_pdb(fp);
            fclose(fp);
        }
        else if (dot && !strcasecmp(dot, ".sdf"))
        {
            fp = fopen(argv[i], "rb");
            if (!fp)
            {
                cout << "Failed to open " << argv[i] << " for reading." << endl;
                return -1;
            }
            got = fread(buffer, 1, 49150, fp);
            fclose(fp);
            m.from_sdf(buffer);
        }
        else if (dot && !strcasecmp(dot, ".cvty"))
        {
            cvtyfname = argv[i];
        }
        else if (!strcasecmp(argv[i], "aoi"))
        {
            for (i++; i<argc; i++)
            {
                if (!strcasecmp(argv[i], "ioa")) break;
                else if (!strcasecmp(argv[i], "*")) break;
                else atoms_of_interest.push_back(argv[i]);
            }
        }
        else if (!strcasecmp(argv[i], "cen"))
        {
            for (i++; i<argc; i++)
            {
                if (!strcasecmp(argv[i], "nec")) break;
                else if (!strcasecmp(argv[i], "*")) break;
                else
                {
                    if (argv[i][strlen(argv[i])-1] == '!') priorities[pocketcen_res.size()] = true;
                    pocketcen_res.push_back(argv[i]);
                }
            }
        }
        else if (!strcasecmp(argv[i], "full"))
        {
            fullscan = true;
            if (i+1 < argc && argv[i+1][0] >= '0' && argv[i+1][0] <= '9') maxbbr = atoi(argv[++i]);
        }
        else if (!strcasecmp(argv[i], "savtmp"))
        {
            save_tmp_pdbs = true;
        }
        else if (!strcasecmp(argv[i], "top"))
        {
            showtopn = atoi(argv[++i]);
        }
        else cout << "Warning: unknown command argument " << argv[i] << endl;
    }

    if (!p.get_end_resno())
    {
        cout << "ERROR: No protein." << endl;
        return -2;
    }

    if (!m.get_atom_count())
    {
        cout << "ERROR: No ligand." << endl;
        return -2;
    }

    if (cvtyfname && strlen(cvtyfname))
    {
        if (!file_exists(cvtyfname))
        {
            cerr << "ERROR file not found: " << cvtyfname << endl;
            return -7;
        }
        FILE* fp = fopen(cvtyfname, "rb");
        if (!fp)
        {
            cerr << "FAILED to open " << cvtyfname << " for reading." << endl;
            return -7;
        }

        cout << "Reading " << cvtyfname << "..." << endl;
        char buffer[1024];
        while (!feof(fp))
        {
            char* got = fgets(buffer, 1022, fp);
            if (!got) break;
            CPartial cp;
            int cno = cp.from_cvty_line(buffer);
            cvtys[cno].add_partial(cp);
            if (cno+1 > ncvtys) ncvtys = cno+1;
        }
        fclose(fp);
        cout << "Read " << ncvtys << " cavities." << endl;
    }

    Point pocketcen = p.get_region_center(1, p.get_end_resno());
    size = Point(100,100,100);

    int maxlt = m.get_heavy_atom_count()+8;
    LigandTarget lt[maxlt];
    int ltargs = Search::identify_ligand_pairing_targets(&m, lt, maxlt);
    cout << ltargs << " ligand target(s) found:" << endl;
    for (i=0; i<ltargs; i++)
    {
        if (lt[i].single_atom) cout << "Atom " << lt[i].single_atom->name << endl;
        else if (lt[i].conjgrp) cout << "Group " << *lt[i].conjgrp << endl;
        else cout << "(null)" << endl;
    }
    cout << endl;
    lt[ltargs].conjgrp = nullptr;
    lt[ltargs].single_atom = nullptr;

    if (fullscan)
    {
        BestBindingResult bbrs[maxbbr+2];

        AminoAcid* aafloor = p.get_residue_bw(3, 40);
        Point ptfloor = aafloor ? aafloor->get_CA_location() : p.get_region_center(1, p.get_end_resno());
        Box b(-Avogadro, Avogadro, ptfloor.y, Avogadro, -Avogadro, Avogadro);

        Search::scan_protein(&p, &m, lt, bbrs, maxbbr, b);

        for (i=0; i<maxbbr; i++)
        {
            if (!bbrs[i].probability) continue;
            if (bbrs[i].pri_res && bbrs[i].pri_tgt)
            {
                cout << "Paired: " << bbrs[i].pri_res->get_name();
                BallesterosWeinstein bw = p.get_bw_from_resno(bbrs[i].pri_res->get_residue_no());
                if (bw.helix_no) cout << "(" << bw << ")";
                cout << "...";
                if (bbrs[i].pri_tgt->single_atom) cout << bbrs[i].pri_tgt->single_atom->name;
                else if (bbrs[i].pri_tgt->conjgrp) cout << *(bbrs[i].pri_tgt->conjgrp);
                cout << endl;
            }
            if (bbrs[i].sec_res && bbrs[i].sec_tgt)
            {
                cout << "Paired: " << bbrs[i].sec_res->get_name();
                BallesterosWeinstein bw = p.get_bw_from_resno(bbrs[i].sec_res->get_residue_no());
                if (bw.helix_no) cout << "(" << bw << ")";
                cout << "...";
                if (bbrs[i].sec_tgt->single_atom) cout << bbrs[i].sec_tgt->single_atom->name;
                else if (bbrs[i].sec_tgt->conjgrp) cout << *(bbrs[i].sec_tgt->conjgrp);
                cout << endl;
            }
            if (bbrs[i].tert_res && bbrs[i].tert_tgt)
            {
                cout << "Paired: " << bbrs[i].tert_res->get_name();
                BallesterosWeinstein bw = p.get_bw_from_resno(bbrs[i].tert_res->get_residue_no());
                if (bw.helix_no) cout << "(" << bw << ")";
                cout << "...";
                if (bbrs[i].tert_tgt->single_atom) cout << bbrs[i].tert_tgt->single_atom->name;
                else if (bbrs[i].tert_tgt->conjgrp) cout << *(bbrs[i].tert_tgt->conjgrp);
                cout << endl;
            }
            else if (bbrs[i].tert_res)
            {
                cout << "Extra: " << bbrs[i].tert_res->get_name();
                BallesterosWeinstein bw = p.get_bw_from_resno(bbrs[i].tert_res->get_residue_no());
                if (bw.helix_no) cout << "(" << bw << ")";
                cout << endl;
            }
            cout << "Probability: " << bbrs[i].probability << endl;
            cout << endl;
        }
    }
    else
    {
        if (p.get_bw50(6))
        {
            int pcn = pocketcen_res.size();
            if (!pcn)
            {
                pocketcen_res.push_back("3.33");
                pocketcen_res.push_back("3.37");
                pocketcen_res.push_back("5.43");
                pocketcen_res.push_back("6.48");
                pocketcen_res.push_back("7.38");
                pcn = pocketcen_res.size();
            }

            AminoAcid* aa;
            Point pt4avg[pcn+2];
            
            int j=0;
            for (i=0; i<pcn; i++)
            {
                if (pocketcen_res[i].member_no) aa = p.get_residue(pocketcen_res[i]);
                else aa = p.get_residue(pocketcen_res[i].helix_no);         // a trick to allow plain integer resnos on the command line.
                if (aa)
                {
                    aa->priority = priorities[i];
                    pt4avg[j++] = aa->get_CA_location();
                }
            }

            pocketcen = average_of_points(pt4avg, j);
            size = size_of_point_space(pt4avg, j);
        }

        loneliest = p.find_loneliest_point(pocketcen, size);
        cout << "Loneliest point = " << loneliest << endl;

        AminoAcid* reaches_spheroid[256];

        Cavity* cv = nullptr;
        float cvr = Avogadro;
        if (ncvtys) for (i=0; i<ncvtys; i++)
        {
            float r = cvtys[i].get_center().get_3d_distance(pocketcen);
            if (r < cvr)
            {
                cv = &cvtys[i];
                cvr = r;
            }
        }

        m.recenter(pocketcen);
        p.get_residues_can_clash_ligand(reaches_spheroid, &m, pocketcen, Point(7,7,7), nullptr);

        BestBindingResult bbr;
        Search::pair_targets(&m, lt, reaches_spheroid, loneliest, &bbr);
        if (bbr.pri_res && bbr.pri_tgt)
        {
            cout << "Paired: " << bbr.pri_res->get_name() << "...";
            if (bbr.pri_tgt->single_atom) cout << bbr.pri_tgt->single_atom->name;
            else if (bbr.pri_tgt->conjgrp) cout << *(bbr.pri_tgt->conjgrp);
            cout << endl;
        }
        if (bbr.sec_res && bbr.sec_tgt)
        {
            cout << "Paired: " << bbr.sec_res->get_name() << "...";
            if (bbr.sec_tgt->single_atom) cout << bbr.sec_tgt->single_atom->name;
            else if (bbr.sec_tgt->conjgrp) cout << *(bbr.sec_tgt->conjgrp);
            cout << endl;
        }
        if (bbr.tert_res && bbr.tert_tgt)
        {
            cout << "Paired: " << bbr.tert_res->get_name() << "...";
            if (bbr.tert_tgt->single_atom) cout << bbr.tert_tgt->single_atom->name;
            else if (bbr.tert_tgt->conjgrp) cout << *(bbr.tert_tgt->conjgrp);
            cout << endl;
        }
    }

    if (showtopn)
    {
        cout << endl << "Top candidates:" << endl;
        for (i=0; i<showtopn; i++)
        {
            if (!bbr_candidates[i].pri_res || !bbr_candidates[i].pri_tgt)
            {
                cout << "End of list." << endl;
                break;
            }
            cout << bbr_candidates[i] << "score: " << bbr_candidates[i].cached_score << endl << endl;
        }
        cout << endl;
    }

    return 0;
}