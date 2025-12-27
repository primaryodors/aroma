
#include "reshape.h"
#include "dynamic.h"

void Reshape::load_rshpm_file(ReshapeType rshpt, Molecule* ligand)
{
    char infname[256];
    switch (rshpt)
    {
        case rshp_GPCR:
        strcpy(infname, "data/gpcr.rshpm");
        break;

        default:
        cerr << "Unknown reshape type." << endl;
        throw 0xbadc0de;
    }

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
            nrshpm++;
        }
        else if (!strcmp(fields[0], "PROX"))
        {
            m_rshpm[nrshpm].rshpmt = rshpm_prox;
            m_rshpm[nrshpm].rap_index.set(fields[1]);
            m_rshpm[nrshpm].tgtdist = atof(fields[2]);
            m_rshpm[nrshpm].rap_target.set(fields[3]);
            nrshpm++;
        }
        else if (!strcmp(fields[0], "TRNL"))
        {
            m_rshpm[nrshpm].rshpmt = rshpm_xlate;
            m_rshpm[nrshpm].rap_start.set(fields[1]);
            m_rshpm[nrshpm].rap_end.set(fields[2]);
            m_rshpm[nrshpm].rap_index.set(fields[3]);
            m_rshpm[nrshpm].tgtdist = atof(fields[4]);
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
            nrshpm++;
        }
        else if (!strcmp(fields[0], "DEL"))
        {
            m_rshpm[nrshpm].rshpmt = rshpm_delete;
            m_rshpm[nrshpm].rap_start.set(fields[1]);
            m_rshpm[nrshpm].rap_end.set(fields[2]);
            nrshpm++;
        }
        else if (!strcmp(fields[0], "WIND"))
        {
            m_rshpm[nrshpm].rshpmt = rshpm_wind;
            m_rshpm[nrshpm].rap_start.set(fields[1]);
            m_rshpm[nrshpm].rap_end.set(fields[2]);
            m_rshpm[nrshpm].rap_index.set(fields[3]);
            m_rshpm[nrshpm].tgtdist = atof(fields[4]);
            m_rshpm[nrshpm].rap_target.set(fields[5]);
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
            m_rshpm[nrshpm].rap_target.set(fields[6]);
            nrshpm++;
        }
        else if (!strcmp(fields[0], "EXIT"))
        {
            break;
        }
    }
    fclose(fp);
}

bool ReshapeMotion::get_pt_index_and_tgt(Protein* p, Point* index, Point* target)
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
        #if _dbg_rshpm_apply
        cout << "Target is ligand." << endl;
        #endif
    }

    if (entire)
    {
        closest_to_ligand = p->get_nearest_atom(*target, rap_start.resno, rap_end.resno);
        if (!closest_to_ligand) return false;
        *index = closest_to_ligand->loc;
        #if _dbg_rshpm_apply
        cout << "Index is entire." << endl;
        #endif
    }

    return true;
}

bool ReshapeMotion::fix_clash(Protein* p, int sr, int er, Point pt_fulcrum, int iters, float step)
{
    #if _dbg_rshpm_apply
    // cout << "Motion fixes clash." << endl;
    #endif
    int i;
    Molecule *aai = p->get_residue(rap_index.resno), *aat = p->get_residue(rap_target.resno);
    if (entire) aai = (Molecule*)closest_to_ligand->mol;
    if (tgtligand) aat = ligand;
    if (!aai || !aat)
    {
        cerr << "Something went wrong in ReshapeMotion::fix_clash()." << endl;
        throw 0xbadc0de;
    }

    float oldclash = aai->get_intermol_clashes(aat), clash = oldclash;
    cout << "Rotating " << sr << "->" << er << " about " << rap_fulcrum.resno 
        << " to avoid clash of " << oldclash << "kJ/mol..." << flush;
    for (i=0; i<iters; i++)
    {
        Rotation rot = align_points_3d(aat->get_barycenter(), aai->get_barycenter(), pt_fulcrum);
        p->rotate_piece(sr, er, pt_fulcrum, rot.v, step);
        cout << ".";
        clash = aai->get_intermol_clashes(aat);
        // if (!i && clash > oldclash) step *= -1;
        if (clash <= clash_limit_per_aa) break;
    }
    cout << endl;
    return (clash < oldclash);
}

bool ReshapeMotion::set_distance(Protein* p, int sr, int er, Point pt_fulcrum, Point pt_index, Point pt_target, int moreorless, float amount)
{
    Vector tolerance = pt_index.subtract(pt_target);
    Rotation rot;

    tolerance.r = tgtdist;
    pt_target = pt_target.add(tolerance);
    float r = pt_index.get_3d_distance(pt_target);

    if ((moreorless > 0 && r < tgtdist)
        || (moreorless < 0 && r > tgtdist)
        || !moreorless
        )
    {
        rot = align_points_3d(pt_index, pt_target, pt_fulcrum);
        rot.a *= amount;
        cout << "Rotating " << sr << "->" << er << " " 
            << (rot.a*fiftyseven) << " deg about " << rap_fulcrum.resno << " to minimum residue distance " << tolerance.r << "A..." << endl;
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
        if (post == ool) m_rshpm[i].apply(p);
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

        if (fixclash) fix_clash(p, rap_start.resno, rap_end.resno, pt_fulcrum);
        else set_distance(p, rap_start.resno, rap_end.resno, pt_fulcrum, pt_index, pt_target, morethan?1:-1, 1);
    }
    else if (rshpmt == rshpm_bend)
    {
        rap_start.resolve_resno(p);
        if (!rap_start.resno) return;
        rap_end.resolve_resno(p);
        if (!rap_end.resno) return;
        rap_index.resolve_resno(p);
        if (!rap_index.resno) return;

        Point pt_fulcrum, pt_index, pt_target;

        int i, j, sign, m, n;
        m = n = abs(rap_end.resno - rap_start.resno);
        sign = sgn(rap_end.resno - rap_start.resno);
        if (n < 2) return;

        for (i = rap_start.resno; i != rap_end.resno; i += sign)
        {
            if (!get_pt_index_and_tgt(p, &pt_index, &pt_target)) return;
            AminoAcid* aa_fulcrum = p->get_residue(i);
            if (!aa_fulcrum) continue;
            rap_fulcrum.resno = i;
            pt_fulcrum = aa_fulcrum->get_CA_location();

            float factor = 1.0 / (n-1);
            if (fixclash)
            {
                bool result = fix_clash(p, min(i, rap_end.resno), max(i, rap_end.resno), pt_fulcrum, 5, (1.0/m)*fiftyseventh);
                if (!result && abs(i-rap_start.resno) > 2) break;
            }
            else set_distance(p, min(i, rap_end.resno), max(i, rap_end.resno), pt_fulcrum, pt_index, pt_target, morethan?1:-1, factor);
            n--;
        }
        cout << "Bent " << rap_start.resno << "->" << rap_end.resno;
        if (fixclash)
        {
            cout << " to fix clash between ";
            if (entire) cout << "nearest side chain";
            else cout << rap_index.resno;
            cout << " and ";
            if (tgtligand) cout << "ligand";
            else cout << rap_target.resno;
        }
        else
        {
            cout << " to bring ";
            if (entire) cout << "nearest side chain";
            else cout << rap_index.resno;
            cout << (morethan ? " at least " : " within ") << tgtdist << " A from ";
            if (tgtligand) cout << "ligand";
            else cout << rap_target.resno;
        }
        cout << endl;
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
                cout << "Pointing " << aa1->get_name() << ":" << a1->name
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

        p->delete_residues(rap_start.resno, rap_end.resno);
        cout << "Deleting residues " << rap_start.resno << " - " << rap_end.resno << endl;
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
                if (v.r) p->move_piece(rap_start.resno, rap_end.resno, v);
                cout << "Translated " << rap_start.resno << " - " << rap_end.resno 
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
        cout << "Wind ";

        if (aa1 && aa2)
        {
            cout << rap_start.resno << " - " << rap_end.resno << " to bring " << aa1->get_name() << " and " << aa2->get_name() << " " << flush;
            a1 = aa1->get_atom(rap_index.get_aname().c_str());
            a2 = aa2->get_atom(rap_target.get_aname().c_str());

            if (a1 && a2)
            {
                cout << a1->name << " and " << a2->name << " ";
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
                    cout << "." << flush;                }
                cout << "Wound " << rap_start.resno << " - " << rap_end.resno 
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
                                if (rnew > rold) theta *= -0.666;
                                else if (rnew < tgtdist+0.1) break;
                                total_rot += theta;
                            }
                        }
                        else if (b->can_flip)
                        {
                            cout << "and can flip." << endl;
                            rold = a1->distance_to(a2);
                            b->rotate(b->flip_angle);
                            total_rot = b->flip_angle;
                            rnew = a1->distance_to(a2);
                            if (rnew > rold)
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

                        cout << "Rotated " << aa1->get_name() << ":" << b->atom1->name << "-" << b->atom2->name
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
