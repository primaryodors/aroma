
#include "reshape.h"
#include "dynamic.h"

bool rshp_verbose = false;

void Reshape::load_rshpm_file(ReshapeType rshpt, Molecule* ligand)
{
    char infname[256];
    switch (rshpt)
    {
        case rshp_GPCR:
        strcpy(infname, "data/gpcr.rshpm");
        load_rshpm_file(infname, ligand);
        break;

        default:
        cerr << "Unknown reshape type." << endl;
        throw 0xbadc0de;
    }
}

void Reshape::load_rshpm_file(const char* infname, Molecule* ligand)
{
    FILE* fp = fopen(infname, "rb");
    if (!fp)
    {
        cerr << "Failed to open " << infname << " for reading." << endl;
        throw 0xbadf12e;
    }

    char buffer[1024];
    m_rshpm = new ReshapeMotion[1024];
    while (!feof(fp))
    {
        char* ln = fgets(buffer, 1020, fp);
        if (!ln || feof(fp)) break;
        char* hash = strchr(buffer, '#');
        if (hash) *hash = 0;
        char** fields = chop_spaced_words(buffer);
        if (!fields || !fields[0] || !strlen(fields[0])) continue;

        if (!strcmp(fields[0], "PIVOT"))                    // https://www.youtube.com/watch?v=n67RYI_0sc0
        {
            m_rshpm[nrshpm].rshpmt = rshpm_pivot;
            m_rshpm[nrshpm].rap_start.set(fields[1]);
            m_rshpm[nrshpm].rap_end.set(fields[2]);
            m_rshpm[nrshpm].rap_fulcrum.set(fields[3]);
            if (!strcmp(fields[4], "entire"))
                m_rshpm[nrshpm].entire = true;
            else m_rshpm[nrshpm].rap_index.set(fields[4]);
            m_rshpm[nrshpm].tgtdist = atof(fields[5]);
            if (!strcmp(fields[5], "noclash"))
            {
                m_rshpm[nrshpm].fixclash = true;
            }
            else if (fields[5][0] == '>')
            {
                m_rshpm[nrshpm].morethan = true;
                m_rshpm[nrshpm].tgtdist = atof(&(fields[5][1]));
            }
            if (!strcmp(fields[6], "ligand"))
            {
                m_rshpm[nrshpm].tgtligand = true;
                m_rshpm[nrshpm].ligand = ligand;
            }
            else m_rshpm[nrshpm].rap_target.set(fields[6]);
            if (fields[7])
            {
                if (!strcmp(fields[7], "MEAS")) m_rshpm[nrshpm].measure = true;
            }
            nrshpm++;
        }
        else if (!strcmp(fields[0], "BEND"))
        {
            m_rshpm[nrshpm].rshpmt = rshpm_bend;
            m_rshpm[nrshpm].rap_start.set(fields[1]);
            m_rshpm[nrshpm].rap_end.set(fields[2]);
            m_rshpm[nrshpm].rap_index.set(fields[3]);
            m_rshpm[nrshpm].tgtdist = atof(fields[4]);
            /*if (!strcmp(fields[4], "noclash"))
            {
                m_rshpm[nrshpm].fixclash = true;
            }
            else*/ if (fields[4][0] == '>')
            {
                m_rshpm[nrshpm].morethan = true;
                m_rshpm[nrshpm].tgtdist = atof(&(fields[4][1]));
            }
            m_rshpm[nrshpm].rap_target.set(fields[5]);
            if (fields[6])
            {
                if (!strcmp(fields[6], "MEAS")) m_rshpm[nrshpm].measure = true;
            }
            nrshpm++;
        }
        else if (!strcmp(fields[0], "PROX"))
        {
            m_rshpm[nrshpm].rshpmt = rshpm_prox;
            m_rshpm[nrshpm].rap_index.set(fields[1]);
            m_rshpm[nrshpm].tgtdist = atof(fields[2]);
            m_rshpm[nrshpm].rap_target.set(fields[3]);
            if (fields[4])
            {
                if (!strcmp(fields[4], "MEAS")) m_rshpm[nrshpm].measure = true;
            }
            nrshpm++;
        }
        else if (!strcmp(fields[0], "TRNL"))
        {
            m_rshpm[nrshpm].rshpmt = rshpm_xlate;
            m_rshpm[nrshpm].rap_start.set(fields[1]);
            m_rshpm[nrshpm].rap_end.set(fields[2]);
            m_rshpm[nrshpm].rap_index.set(fields[3]);
            m_rshpm[nrshpm].tgtdist = atof(fields[4]);
            if (fields[4][0] == '>')
            {
                m_rshpm[nrshpm].morethan = true;
                m_rshpm[nrshpm].tgtdist = atof(&(fields[4][1]));
            }
            if (!strcmp(fields[4], "noclash"))
            {
                m_rshpm[nrshpm].fixclash = true;
            }
            else if (fields[4][0] == '>')
            {
                m_rshpm[nrshpm].morethan = true;
                m_rshpm[nrshpm].tgtdist = atof(&(fields[4][1]));
            }
            if (!strcmp(fields[5], "ligand"))
            {
                m_rshpm[nrshpm].tgtligand = true;
                m_rshpm[nrshpm].ligand = ligand;
            }
            else m_rshpm[nrshpm].rap_target.set(fields[5]);
            if (fields[6])
            {
                if (!strcmp(fields[6], "MEAS")) m_rshpm[nrshpm].measure = true;
            }
            nrshpm++;
        }
        else if (!strcmp(fields[0], "DEL"))
        {
            m_rshpm[nrshpm].rshpmt = rshpm_delete;
            m_rshpm[nrshpm].rap_start.set(fields[1]);
            m_rshpm[nrshpm].rap_end.set(fields[2]);
            if (fields[3])
            {
                if (!strcmp(fields[3], "MEAS")) m_rshpm[nrshpm].measure = true;
            }
            nrshpm++;
        }
        else if (!strcmp(fields[0], "WIND"))
        {
            m_rshpm[nrshpm].rshpmt = rshpm_wind;
            m_rshpm[nrshpm].rap_start.set(fields[1]);
            m_rshpm[nrshpm].rap_end.set(fields[2]);
            m_rshpm[nrshpm].rap_index.set(fields[3]);
            m_rshpm[nrshpm].tgtdist = atof(fields[4]);
            if (fields[4][0] == '>')
            {
                m_rshpm[nrshpm].morethan = true;
                m_rshpm[nrshpm].tgtdist = atof(&(fields[4][1]));
            }
            m_rshpm[nrshpm].rap_target.set(fields[5]);
            if (fields[6])
            {
                if (!strcmp(fields[6], "MEAS")) m_rshpm[nrshpm].measure = true;
            }
            nrshpm++;
        }
        else if (!strcmp(fields[0], "FLEX"))
        {
            m_rshpm[nrshpm].rshpmt = rshpm_flex;
            m_rshpm[nrshpm].rap_start.set(fields[1]);
            if (strlen(fields[2]) > 13) fields[2][13] = 0;
            strcpy(m_rshpm[nrshpm].ba1, fields[2]);
            if (strlen(fields[3]) > 13) fields[3][13] = 0;
            strcpy(m_rshpm[nrshpm].ba2, fields[3]);
            m_rshpm[nrshpm].rap_index.set(fields[4]);
            m_rshpm[nrshpm].tgtdist = atof(fields[5]);
            if (fields[5][0] == '>')
            {
                m_rshpm[nrshpm].morethan = true;
                m_rshpm[nrshpm].tgtdist = atof(&(fields[5][1]));
            }
            m_rshpm[nrshpm].rap_target.set(fields[6]);
            if (fields[7])
            {
                if (!strcmp(fields[7], "MEAS")) m_rshpm[nrshpm].measure = true;
            }
            nrshpm++;
        }
        else if (!strcmp(fields[0], "EXIT"))
        {
            break;
        }
    }
    fclose(fp);
}

bool ReshapeMotion::get_pt_index_and_tgt(Protein* p, Point* index, Point* target, Atom** a_index, Atom** a_target)
{
    if (!entire)
    {
        rap_index.resolve_resno(p);
        if (rap_index.get_aname().c_str()[0] == 'n')
            rap_index.resolve_special_atom(p, rap_target.loc());
        if (!rap_index.resno) return false;
        #if _dbg_rshpm_apply
        cout << "Resolved index " << rap_index.resno << endl;
        #endif
        *index = rap_index.loc();
        if (a_index) *a_index = rap_index.atom();
    }
    if (!tgtligand)
    {
        rap_target.resolve_resno(p);
        if (rap_target.get_aname().c_str()[0] == 'n')
            rap_target.resolve_special_atom(p, rap_index.loc());
        if (!rap_target.resno) return false;
        #if _dbg_rshpm_apply
        cout << "Resolved index " << rap_target.resno << endl;
        #endif
        *target = rap_target.loc();
        if (a_target) *a_target = rap_target.atom();
    }
    if (tgtligand && !ligand)
    {
        cerr << "Called ligand-target reshape motion without a ligand set." << endl;
        throw 0xbadc0de;                // bad code monkey no banana!
    }
    if (tgtligand)
    {
        Atom* liga = ligand->get_nearest_atom_to_line(rap_start.loc(), rap_end.loc());
        if (!liga) return false;
        *target = liga->loc;
        if (a_target) *a_target = liga;
        #if _dbg_rshpm_apply
        cout << "Target is ligand." << endl;
        #endif
    }
    if (entire)
    {
        closest_to_ligand = p->get_nearest_atom(*target, min(rap_start.resno, rap_end.resno), max(rap_start.resno, rap_end.resno));
        if (!closest_to_ligand)
        {
            cout << "No atom near " << *target << " between " << rap_start.resno << " and " << rap_end.resno << endl;
            return false;
        }
        *index = closest_to_ligand->loc;
        if (a_index) *a_index = closest_to_ligand;
        #if _dbg_rshpm_apply
        cout << "Index is entire." << endl;
        #endif
    }

    return true;
}

float ReshapeMotion::measure_index_tgt_clashes(Protein *p)
{
    Molecule* ltarget;
    if (tgtligand) ltarget = ligand;
    else
    {
        if (!rap_target.resno) rap_target.resolve_resno(p);
        ltarget = (Molecule*)(p->get_residue(rap_target.resno));
    }
    if (!ltarget) return 0.0f;

    if (entire)
        return p->get_intermol_clashes(ltarget, min(rap_start.resno, rap_end.resno), max(rap_start.resno, rap_end.resno));

    if (!rap_index.resno) rap_index.resolve_resno(p);
    AminoAcid* lindex = p->get_residue(rap_index.resno);

    return ltarget->get_intermol_clashes(lindex);
}

bool ReshapeMotion::fix_clash(Protein* p, int sr, int er, Point pt_fulcrum, int iters, float step)
{
    #if _dbg_rshpm_apply
    // cout << "Motion fixes clash." << endl;
    #endif
    int i;

    float oldclash = measure_index_tgt_clashes(p), clash = oldclash;
    if (rshp_verbose) cout << "Rotating " << sr << "->" << er << " about " << rap_fulcrum.resno
        << " by " << (step*fiftyseven) << " deg "
        << " to avoid clash of " << oldclash << "kJ/mol..." << flush;
    for (i=0; i<iters; i++)
    {
        Point ptidx, pttgt;
        if (!get_pt_index_and_tgt(p, &ptidx, &pttgt)) return false;
        Rotation rot = align_points_3d(pttgt, ptidx, pt_fulcrum);
        p->rotate_piece(sr, er, pt_fulcrum, rot.v, step);
        if (rshp_verbose) cout << ".";
        clash = measure_index_tgt_clashes(p);
        step *= 0.97;
        if (clash <= clash_limit_per_aa) break;
    }
    if (rshp_verbose) cout << endl;
    return (clash < oldclash);
}

bool ReshapeMotion::set_distance(Protein* p, int sr, int er, Point pt_fulcrum, Point pt_index, Point pt_target, int moreorless, float amount)
{
    float r = pt_index.get_3d_distance(pt_target);

    if ((moreorless > 0 && r < tgtdist)
        || (moreorless < 0 && r > tgtdist)
        || !moreorless
        )
    {
        Vector tolerance = pt_index.subtract(pt_target);
        Rotation rot;
        tolerance.r = tgtdist;
        pt_target = pt_target.add(tolerance);
        rot = align_points_3d(pt_index, pt_target, pt_fulcrum);
        rot.a *= amount;
        if (rshp_verbose) cout << "Rotating " << sr << "->" << er << " " 
            << (rot.a*fiftyseven) << " deg about " << rap_fulcrum.resno
            << " to minimum residue distance " << tolerance.r << "A"
            << " from previous " << r
            << endl;
        p->rotate_piece(sr, er, pt_fulcrum, rot.v, rot.a);
    }

    return true;
}

void Reshape::apply(Protein *p, bool ool)            // This is our ool. There is no p in it. 🌊🏊🌊
{
    int i;
    bool post = false;
    for (i=0; i<nrshpm; i++)
    {
        if (m_rshpm[i].tgtligand) post = true;
        if (post == ool)
        {
            m_rshpm[i].apply(p);
            if (m_rshpm[i].measure)
                m_rshpm[i].do_measurement(p);
        }
    }
}

void ReshapeMotion::apply(Protein *p)
{
    if (rshpmt == rshpm_pivot)
    {
        #if _dbg_rshpm_apply
        cout << "Pivot motion." << endl;
        #endif
        rap_start.resolve_resno(p);
        if (!rap_start.resno) return;
        rap_end.resolve_resno(p);
        if (!rap_end.resno) return;
        #if _dbg_rshpm_apply
        cout << "Resolved start resno " << rap_start.resno << " and end " << rap_end.resno << endl;
        #endif
        rap_fulcrum.resolve_resno(p);
        if (!rap_fulcrum.resno) return;
        #if _dbg_rshpm_apply
        cout << "Resolved fulcrum " << rap_fulcrum.resno << endl;
        #endif

        Point pt_fulcrum = rap_fulcrum.loc(), pt_index, pt_target;
        if (!get_pt_index_and_tgt(p, &pt_index, &pt_target)) return;
        #if _dbg_rshpm_apply
        cout << "Resolved target and index." << endl;
        #endif

        if (fixclash)
        {
            #if _dbg_rshpm_apply
            cout << "Fixing clash." << endl;
            #endif
            fix_clash(p, rap_start.resno, rap_end.resno, pt_fulcrum);
        }
        else
        {
            #if _dbg_rshpm_apply
            cout << "Aiming for distance." << endl;
            #endif
            set_distance(p, rap_start.resno, rap_end.resno, pt_fulcrum, pt_index, pt_target, morethan?1:-1, 1);
        }
    }
    else if (rshpmt == rshpm_bend)
    {
        #if _dbg_rshpm_apply
        cout << "Bend motion." << endl;
        #endif
        rap_start.resolve_resno(p);
        if (!rap_start.resno) return;
        rap_end.resolve_resno(p);
        if (!rap_end.resno) return;
        #if _dbg_rshpm_apply
        cout << "Resolved start resno " << rap_start.resno << " and end " << rap_end.resno << endl;
        #endif
        rap_index.resolve_resno(p);
        if (!rap_index.resno && !entire) return;
        #if _dbg_rshpm_apply
        if (entire) cout << "Entire." << endl;
        else cout << "Resolved index." << endl;
        #endif

        Point pt_fulcrum, pt_index, pt_target;

        int i, j, sign, m, n;
        m = n = abs(rap_end.resno - rap_start.resno);
        sign = sgn(rap_end.resno - rap_start.resno);
        if (n < 2) return;
        #if _dbg_rshpm_apply
        cout << "Sufficient residues for getting bent." << endl;
        #endif

        // Atom* pa = p->get_nearest_atom(ligand->get_barycenter(), min(rap_start.resno, rap_end.resno), max(rap_start.resno, rap_end.resno));
        for (j=0; j<(fixclash?6:1); j++)
        {
            for (i = rap_start.resno; i != rap_end.resno; i += sign)
            {
                if (rshp_verbose) cout << "Bending at " << i << endl;
                // if (pa && i == pa->residue) break;
                if (!get_pt_index_and_tgt(p, &pt_index, &pt_target)) return;
                #if _dbg_rshpm_apply
                cout << "Have index and target points." << endl;
                #endif
                AminoAcid* aa_fulcrum = p->get_residue(i);
                if (!aa_fulcrum) continue;
                #if _dbg_rshpm_apply
                cout << "Fulcrum is " << aa_fulcrum->get_name() << endl;
                #endif
                rap_fulcrum.resno = i;
                pt_fulcrum = aa_fulcrum->get_CA_location();

                float factor = 1.0 / (n-1);
                if (fixclash)
                {
                    bool result = fix_clash(p, min(i, rap_end.resno), max(i, rap_end.resno), pt_fulcrum, 10, (2.5/m)*fiftyseventh);
                    if (!result && abs(i-rap_start.resno) > 5) break;
                }
                else set_distance(p, min(i, rap_end.resno), max(i, rap_end.resno), pt_fulcrum, pt_index, pt_target, morethan?1:-1, factor);
                n--;
            }
            if (measure_index_tgt_clashes(p) < clash_limit_per_aa) break;
        }
        if (rshp_verbose) cout << "Bent " << rap_start.resno << "->" << rap_end.resno;
        if (rshp_verbose && fixclash)
        {
            cout << " to fix clash between ";
            if (entire) cout << "nearest side chain";
            else cout << rap_index.resno;
            cout << " and ";
            if (tgtligand) cout << "ligand";
            else cout << rap_target.resno;
        }
        else if (rshp_verbose)
        {
            cout << " to bring ";
            if (entire) cout << "nearest side chain";
            else cout << rap_index.resno;
            cout << (morethan ? " at least " : " within ") << tgtdist << " A from ";
            if (tgtligand) cout << "ligand";
            else cout << rap_target.resno;
        }
        if (rshp_verbose) cout << endl;
    }
    else if (rshpmt == rshpm_prox)
    {
        rap_index.resolve_resno(p);
        if (!rap_index.resno) return;
        rap_target.resolve_resno(p);
        if (!rap_target.resno) return;

        if (rap_index.get_aname().c_str()[0] == 'n')
            if (!rap_index.resolve_special_atom(p, rap_target.loc())) return;
        if (rap_target.get_aname().c_str()[0] == 'n')
            if (!rap_target.resolve_special_atom(p, rap_index.loc())) return;

        AminoAcid *aa1, *aa2;
        Atom *a1, *a2;

        aa1 = p->get_residue(rap_index.resno);
        aa2 = p->get_residue(rap_target.resno);

        /* cout << "Want to point " << aa1->get_name() << ":" << rap_index.aname
            << " and " << aa2->get_name() << ":" << rap_target.aname
            << " toward each other." << endl; */

        if (aa1 && aa2)
        {
            a1 = aa1->get_atom(rap_index.get_aname().c_str());
            a2 = aa2->get_atom(rap_target.get_aname().c_str());

            if (a1 && a2)
            {
                aa1->conform_atom_to_location(a1, a2, 20, tgtdist);
                aa2->conform_atom_to_location(a2, a1, 20, tgtdist);
                aa1->conform_atom_to_location(a1, a2, 20, tgtdist);
                aa2->conform_atom_to_location(a2, a1, 20, tgtdist);
                if (rshp_verbose) cout << "Pointing " << aa1->get_name() << ":" << a1->name
                    << " and " << aa2->get_name() << ":" << a2->name
                    << " toward each other." << endl;
            }
        }
    }
    else if (rshpmt == rshpm_delete)
    {
        rap_start.resolve_resno(p);
        if (!rap_start.resno) return;
        rap_end.resolve_resno(p);
        if (!rap_end.resno) return;

        p->delete_residues(min(rap_start.resno, rap_end.resno), max(rap_start.resno, rap_end.resno));
        if (rshp_verbose) cout << "Deleting residues " << rap_start.resno << " - " << rap_end.resno << endl;
    }
    else if (rshpmt == rshpm_xlate)
    {
        rap_start.resolve_resno(p);
        if (!rap_start.resno) return;
        rap_end.resolve_resno(p);
        if (!rap_end.resno) return;

        if (!entire)
        {
            rap_index.resolve_resno(p);
            if (rap_index.get_aname().c_str()[0] == 'n')
                rap_index.resolve_special_atom(p, rap_target.loc());
            if (!rap_index.resno) return;
            if (rap_index.get_aname().c_str()[0] == 'n')
                if (!rap_index.resolve_special_atom(p, rap_target.loc())) return;
        }
        if (!tgtligand)
        {
            rap_target.resolve_resno(p);
            if (rap_target.get_aname().c_str()[0] == 'n')
                rap_target.resolve_special_atom(p, rap_index.loc());
            if (!rap_target.resno) return;
            if (rap_target.get_aname().c_str()[0] == 'n')
                if (!rap_target.resolve_special_atom(p, rap_index.loc())) return;
        }

        AminoAcid *aa1, *aa2;
        Atom *a1, *a2;

        // TODO: fixclash and tgtligand behavior.
        aa1 = p->get_residue(rap_index.resno);
        aa2 = p->get_residue(rap_target.resno);
        // cout << "Translate ";

        if (aa1 && aa2)
        {
            // cout << rap_start.resno << " - " << rap_end.resno << " to bring " << aa1->get_name() << " and " << aa2->get_name() << " " << flush;
            a1 = aa1->get_atom(rap_index.get_aname().c_str());
            a2 = aa2->get_atom(rap_target.get_aname().c_str());

            if (a1 && a2)
            {
                // TODO: morethan behavior.
                Vector v = a2->loc.subtract(a1->loc);
                v.r = fmax(0, v.r - tgtdist);
                if (v.r) p->move_piece(min(rap_start.resno, rap_end.resno), max(rap_start.resno, rap_end.resno), v);
                if (rshp_verbose) cout << "Translated " << rap_start.resno << " - " << rap_end.resno 
                    << " by " << v.r << " A to bring "
                    << aa1->get_name() << ":" << a1->name
                    << " " << a1->distance_to(a2) << "A from "
                    << aa2->get_name() << ":" << a2->name
                    << endl;
            }
        }
    }
    else if (rshpmt == rshpm_wind)
    {
        rap_start.resolve_resno(p);
        if (!rap_start.resno) return;
        rap_end.resolve_resno(p);
        if (!rap_end.resno) return;
        rap_index.resolve_resno(p);
        if (!rap_index.resno) return;
        rap_target.resolve_resno(p);
        if (!rap_target.resno) return;

        if (rap_index.get_aname().c_str()[0] == 'n')
            if (!rap_index.resolve_special_atom(p, rap_target.loc())) return;
        if (rap_target.get_aname().c_str()[0] == 'n')
            if (!rap_target.resolve_special_atom(p, rap_index.loc())) return;

        AminoAcid *aa1, *aa2;
        Atom *a1, *a2;

        aa1 = p->get_residue(rap_index.resno);
        aa2 = p->get_residue(rap_target.resno);
        if (rshp_verbose) cout << "Wind ";

        if (aa1 && aa2)
        {
            if (rshp_verbose) cout << rap_start.resno << " - " << rap_end.resno << " to bring " << aa1->get_name() << " and " << aa2->get_name() << " " << flush;
            a1 = aa1->get_atom(rap_index.get_aname().c_str());
            a2 = aa2->get_atom(rap_target.get_aname().c_str());

            if (a1 && a2)
            {
                if (rshp_verbose) cout << a1->name << " and " << a2->name << " ";
                DynamicMotion d(p);
                d.start_resno.from_string(rap_start.bw.c_str());
                d.end_resno.from_string(rap_end.bw.c_str());
                d.type = dyn_wind;
                d.bias = 15;

                float rold, rnew;
                int i;
                for (i=0; i<100; i++)
                {
                    rold = a1->distance_to(a2);
                    d.apply_incremental(0.1);
                    rnew = a1->distance_to(a2);
                    float anomaly = fabs(rnew-tgtdist);
                    if (anomaly > fabs(rold-tgtdist)) d.bias *= -0.666;
                    else if (anomaly < 0.1) break;
                    if (rshp_verbose) cout << "." << flush;
                }
                if (rshp_verbose) cout << "Wound " << rap_start.resno << " - " << rap_end.resno 
                    << " by " << (d.get_total_applied()*15) << " helical units to bring "
                    << aa1->get_name() << ":" << a1->name
                    << " " << a1->distance_to(a2) << "A from "
                    << aa2->get_name() << ":" << a2->name
                    << endl;
            }
        }
    }
    else if (rshpmt == rshpm_flex)
    {
        rap_start.resolve_resno(p);
        if (!rap_start.resno) return;
        rap_index.resolve_resno(p);
        if (!rap_index.resno) return;
        rap_target.resolve_resno(p);
        if (!rap_target.resno) return;

        if (rap_index.get_aname().c_str()[0] == 'n')
            if (!rap_index.resolve_special_atom(p, rap_target.loc())) return;
        if (rap_target.get_aname().c_str()[0] == 'n')
            if (!rap_target.resolve_special_atom(p, rap_index.loc())) return;

        AminoAcid *aa1, *aa2, *aa3;
        aa1 = p->get_residue(rap_start.resno);
        aa2 = p->get_residue(rap_index.resno);
        aa3 = p->get_residue(rap_target.resno);

        if (aa1 && aa2 && aa3)
        {
            Atom *a1, *a2;
            a1 = aa1->get_atom(ba1);
            a2 = aa1->get_atom(ba2);
            if (a1 && a2)
            {
                Bond* b = a1->get_bond_between(a2);
                if (b && b->atom2)
                {
                    a1 = aa2->get_atom(rap_index.get_aname().c_str());
                    a2 = aa3->get_atom(rap_target.get_aname().c_str());
                    if (a1 && a2)
                    {
                        float rold, rnew, total_rot = 0;
                        if (b->can_rotate)
                        {
                            float theta = fiftyseventh*5;
                            int i;
                            for (i=0; i<1000; i++)
                            {
                                rold = a1->distance_to(a2);
                                b->rotate(theta);
                                rnew = a1->distance_to(a2);
                                if (morethan ? (rold > rnew) : (rnew > rold)) theta *= -0.666;
                                else if (morethan ? (rnew > tgtdist) : (rnew < tgtdist+0.1)) break;
                                total_rot += theta;
                            }
                        }
                        else if (b->can_flip)
                        {
                            rold = a1->distance_to(a2);
                            b->rotate(b->flip_angle);
                            total_rot = b->flip_angle;
                            rnew = a1->distance_to(a2);
                            if (morethan ? (rold > rnew && rnew < tgtdist) : (rnew > rold))
                            {
                                b->rotate(b->flip_angle);
                                total_rot = 0;
                            }
                        }
                        else
                        {
                            cerr << "Bond " << aa1->get_name() << ":" << b->atom1->name << "-" << b->atom2->name
                                << " cannot rotate or flip." << endl;
                            return;
                        }

                        if (rshp_verbose) cout << "Rotated " << aa1->get_name() << ":" << b->atom1->name << "-" << b->atom2->name
                            << " by " << (total_rot*fiftyseven) << " deg to bring "
                            << aa2->get_name() << ":" << a1->name
                            << " " << a1->distance_to(a2) << "A from "
                            << aa3->get_name() << ":" << a2->name
                            << endl;
                    }
                    else
                    {
                        if (!a1) cerr << aa2->get_name() << ":" << rap_index.get_orig_aname() << " not found." << endl;
                        if (!a2) cerr << aa3->get_name() << ":" << rap_target.get_orig_aname() << " not found." << endl;
                    }
                }
                else
                {
                    cerr << aa1->get_name() << ":" << ba1 << " and " << ba2 << " are not bonded." << endl;
                }
            }
            else
            {
                if (!a1) cerr << aa1->get_name() << ":" << ba1 << " not found." << endl;
                if (!a2) cerr << aa1->get_name() << ":" << ba2 << " not found." << endl;
            }
        }
    }
}

void ReshapeMotion::do_measurement(Protein *p)
{
    Point pt_index, pt_target;
    Atom *aidx, *atgt;
    if (get_pt_index_and_tgt(p, &pt_index, &pt_target, &aidx, &atgt))
    {
        if (entire) cout << "Region " << rap_start.resno << "->" << rap_end.resno;
        else cout << rap_index.resno << ":" << (aidx ? aidx->name : "cen");
        cout << " is " << pt_index.get_3d_distance(pt_target) << "A from ";
        if (tgtligand) cout << "ligand" << ":" << (atgt ? atgt->name : "cen");
        else cout << rap_target.resno << ":" << (atgt ? atgt->name : "cen");
        cout << "." << endl;
    }
    // else cout << "Dismal failure can't measure a simple point distance." << endl;
}

LocRotation ICHelixGroup::get_motion(InternalContact *ic, ICHelix *ich, Protein *prot)
{
    LocRotation result;

    if (prot) m_prot = prot;
    else if (!m_prot)
    {
        cerr << "No protein." << endl;
        throw 0xbadca22;
    }

    int i, j, l;
    if (!ich)
    {
        for (i=0; i<n_helix; i++)
        {
            if (helices[i].contains(ic))
            {
                ich = &helices[i];
                break;
            }
        }
        if (!ich)
        {
            cerr << "Internal contact not found in helix group." << endl;
            throw 0xbadca22;
        }
    }

    if (ich->hxstatic) return result;

    ic->res1.resolve_resno(m_prot);
    ic->res2.resolve_resno(m_prot);
    ich->start.resolve_resno(m_prot);
    ich->end.resolve_resno(m_prot);

    Vector av = ic->atom_distance(m_prot);
    if (!av.r) return result;

    if (ich->n_ic == 1)
    {
        // A little trick: if the return value has no vector r and no rotation angle,
        // then the origin represents linear translation coordinates.
        result.a = 0;
        result.v.r = 0;
        if (av.r > ic->r_optimal + ic->tolerance) av.r = av.r - (ic->r_optimal + ic->tolerance);
        else if (av.r < ic->r_optimal - ic->tolerance) av.r = (ic->r_optimal - ic->tolerance) - av.r;
        else av.r = 0;           // already in range? no translation necessary.
        result.origin = av;
        return result;
    }

    int farthest, fdist = 0;
    j = ic->res1.resno;
    for (i=0; i<ich->n_ic; i++)
    {
        ich->ic[i].res1.resolve_resno(prot);
        ich->ic[i].res2.resolve_resno(prot);

        l = ich->ic[i].res1.resno;
        int r = abs(j-l);
        if (r > fdist)
        {
            farthest = i;
            fdist = r;
        }
    }

    // Determine how far to rotate the helix about the farthest contact in order to bring ic within tolerance.
    if (av.r > ic->r_optimal + ic->tolerance) av.r = av.r - (ic->r_optimal + ic->tolerance);
    else if (av.r < ic->r_optimal - ic->tolerance) av.r = (ic->r_optimal - ic->tolerance) - av.r;
    else return result;           // already in range? no rotation necessary.

    AminoAcid* sueruon_rocenon = m_prot->get_residue(ich->ic[i].res1.resno);
    if (!sueruon_rocenon) return result;
    result.origin = sueruon_rocenon->get_CA_location();

    AminoAcid* aa1 = prot->get_residue(ic->res1.resno);
    if (!aa1) return result;
    Atom* a1 = aa1->get_reach_atom(hbond);
    if (!a1) return result;
    Point target = a1->loc.add(av);
    Rotation rot = align_points_3d(a1->loc, target, result.origin);
    result.v = rot.v;
    result.a = rot.a;

    return result;
}

void ICHelixGroup::load_ic_file(const char *filename)
{
    FILE* fp = fopen(filename, "rb");
    if (!fp) return;

    int i, j;
    char buffer[1024];
    char icbw[10];
    while (fgets(buffer, 1022, fp))
    {
        char** words = chop_spaced_words(buffer);

        if (!strcmp(words[0], "DEL"))
        {
            deletion_start[n_deletion].set(words[1]);
            deletion_end[n_deletion].set(words[2]);
            n_deletion++;
        }
        if (!strcmp(words[0], "STATIC"))
        {
            j = atoi(words[1]);
            for (i=0; i<n_helix; i++)
            {
                if (helices[i].hxno == j)
                {
                    helices[i].hxstatic = true;
                    break;
                }
            }
        }
        if (!strcmp(words[0], "CNTCT"))
        {
            InternalContact ictmp;
            ictmp.res1.set(words[1]);
            ictmp.res2.set(words[2]);
            if (words[3])
            {
                ictmp.r_optimal = atof(words[3]);
                if (words[4]) ictmp.tolerance = atof(words[4]);
            }

            j = n_helix;
            for (i=0; i<n_helix; i++)
            {
                if (helices[i].hxno == ictmp.res1.hxno)
                {
                    j = i;
                    break;
                }
            }
            if (j < MAX_ICHX)
            {
                if (helices[j].n_ic < MAX_HXIC)
                {
                    helices[j].ic[helices[j].n_ic++] = ictmp;
                }
                if (j == n_helix)
                {
                    helices[j].hxno = ictmp.res1.hxno;
                    sprintf(icbw, "%d.s", ictmp.res1.hxno);
                    helices[j].start.set(icbw);
                    sprintf(icbw, "%d.e", ictmp.res1.hxno);
                    helices[j].end.set(icbw);
                    n_helix++;
                }
            }
            ResiduePlaceholder swap = ictmp.res1;
            ictmp.res1 = ictmp.res2;
            ictmp.res2 = swap;
            j = n_helix;
            for (i=0; i<n_helix; i++)
            {
                if (helices[i].hxno == ictmp.res1.hxno)
                {
                    j = i;
                    break;
                }
            }
            if (j < MAX_ICHX)
            {
                if (helices[j].n_ic < MAX_HXIC)
                {
                    helices[j].ic[helices[j].n_ic++] = ictmp;
                }
                if (j == n_helix)
                {
                    helices[j].hxno = ictmp.res1.hxno;
                    sprintf(icbw, "%d.s", ictmp.res1.hxno);
                    helices[j].start.set(icbw);
                    sprintf(icbw, "%d.e", ictmp.res1.hxno);
                    helices[j].end.set(icbw);
                    n_helix++;
                }
            }
        }
    }

    fclose(fp);
}

float ICHelixGroup::optimize_helices(Protein *prot, int iters)
{
    if (!prot)
    {
        cerr << "Null protein specified for helix optimization." << endl;
        throw 0xbadca22;
    }
    int iter, i, j;
    m_prot = prot;
    float amount = 1.0 / iters;

    for (i=0; i<n_deletion; i++)
    {
        deletion_start[i].resolve_resno(prot);
        deletion_end[i].resolve_resno(prot);

        if (!deletion_start[i].resno || !deletion_end[i].resno) continue;

        prot->delete_residues(deletion_start[i].resno, deletion_end[i].resno);
    }

    for (iter=0; iter<iters; iter++)
    {
        for (i=0; i<n_helix; i++)
        {
            helices[i].start.resolve_resno(prot);
            helices[i].end.resolve_resno(prot);
            for (j=0; j<helices[i].n_ic; j++)
            {
                LocRotation lr = get_motion(&(helices[i].ic[j]), &(helices[i]), prot);

                if (!lr.a && !lr.v.r && lr.origin.magnitude())
                {
                    Vector mov = lr.origin;
                    mov.r *= amount;
                    prot->move_piece(helices[i].start.resno, helices[i].end.resno, mov);
                }
                else
                {
                    lr.a *= amount;
                    prot->rotate_piece(helices[i].start.resno, helices[i].end.resno, lr);
                }
            }
        }
    }

    float result = 0;
    for (i=0; i<n_helix; i++)
    {
        for (j=0; j<helices[i].n_ic; j++)
        {
            InternalContact* lic = &(helices[i].ic[j]);
            if (lic->res1.resno < lic->res2.resno)
            {
                Vector av = lic->atom_distance(prot);
                if (av.r < lic->r_optimal - lic->tolerance) result += (lic->r_optimal - lic->tolerance - av.r);
                else if (av.r > lic->r_optimal + lic->tolerance) result += (av.r - lic->r_optimal - lic->tolerance);
            }
        }
    }

    return result;
}

bool ICHelix::contains(InternalContact *lic)
{
    int i;
    for (i=0; i<n_ic; i++)
    {
        if (&(ic[i]) == lic) return true;

        if (ic[i].res1.hxno && ic[i].res2.hxno
            && ic[i].res1.hxno == lic->res1.hxno
            && ic[i].res1.bwpos == lic->res1.bwpos
            && ic[i].res2.hxno == lic->res2.hxno
            && ic[i].res2.bwpos == lic->res2.bwpos
            ) return true;
    }
    return false;
}

Vector InternalContact::atom_distance(Protein *prot)
{
    AminoAcid* aa1 = prot->get_residue(res1.resno);
    AminoAcid* aa2 = prot->get_residue(res2.resno);
    Atom* a1 = aa1->get_reach_atom(hbond);
    Atom* a2 = aa2->get_reach_atom(hbond);
    return a2->loc.subtract(a1->loc);
}
