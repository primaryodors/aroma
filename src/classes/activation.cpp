
#include "activation.h"

void Activation::load_acvm_file(AcvType acvt, Molecule* ligand)
{
    char infname[256];
    switch (acvt)
    {
        case acv_GPCR:
        strcpy(infname, "data/gpcr.acvm");
        break;

        default:
        cerr << "Unknown activation type." << endl;
        throw 0xbadc0de;
    }

    FILE* fp = fopen(infname, "rb");
    if (!fp)
    {
        cerr << "Failed to open " << infname << " for reading." << endl;
        throw 0xbadf12e;
    }

    char buffer[1024];
    m_acvm = new ActiveMotion[256];
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
            m_acvm[nacvm].acvmt = acvm_pivot;
            m_acvm[nacvm].rap_start.set(fields[1]);
            m_acvm[nacvm].rap_end.set(fields[2]);
            m_acvm[nacvm].rap_fulcrum.set(fields[3]);
            if (!strcmp(fields[4], "entire"))
                m_acvm[nacvm].entire = true;
            else m_acvm[nacvm].rap_index.set(fields[4]);
            m_acvm[nacvm].tgtdist = atof(fields[5]);
            if (!strcmp(fields[5], "noclash"))
            {
                m_acvm[nacvm].fixclash = true;
            }
            else if (fields[5][0] == '>')
            {
                m_acvm[nacvm].morethan = true;
                m_acvm[nacvm].tgtdist = atof(&(fields[5][1]));
            }
            if (!strcmp(fields[6], "ligand"))
            {
                m_acvm[nacvm].tgtligand = true;
                m_acvm[nacvm].ligand = ligand;
            }
            else m_acvm[nacvm].rap_target.set(fields[6]);
            nacvm++;
        }
        else if (!strcmp(fields[0], "PROX"))
        {
            m_acvm[nacvm].acvmt = acvm_prox;
            m_acvm[nacvm].rap_index.set(fields[1]);
            m_acvm[nacvm].tgtdist = atof(fields[2]);
            m_acvm[nacvm].rap_target.set(fields[3]);
            nacvm++;
        }
    }
    fclose(fp);
}

void Activation::apply(Protein *p, bool ool)            // This is our ool. There is no p in it. 🌊🏊🌊
{
    int i;
    bool post = false;
    for (i=0; i<nacvm; i++)
    {
        if (m_acvm[i].tgtligand) post = true;
        if (post == ool) m_acvm[i].apply(p);
    }
}

void ActiveMotion::apply(Protein *p)
{
    if (acvmt == acvm_pivot)
    {
        #if _dbg_acvm_apply
        cout << "Pivot motion." << endl;
        #endif
        rap_start.resolve_resno(p);
        if (!rap_start.resno) return;
        rap_end.resolve_resno(p);
        if (!rap_end.resno) return;
        #if _dbg_acvm_apply
        cout << "Resolved start resno " << rap_start.resno << " and end " << rap_end.resno << endl;
        #endif
        rap_fulcrum.resolve_resno(p);
        if (!rap_fulcrum.resno) return;
        #if _dbg_acvm_apply
        cout << "Resolved fulcrum " << rap_fulcrum.resno << endl;
        #endif
        if (!entire)
        {
            rap_index.resolve_resno(p);
            if (rap_index.aname.c_str()[0] == 'n')
                rap_index.resolve_special_atom(p, rap_target.loc());
            if (!rap_index.resno) return;
            #if _dbg_acvm_apply
            cout << "Resolved index " << rap_index.resno << endl;
            #endif
        }
        if (!tgtligand)
        {
            rap_target.resolve_resno(p);
            if (rap_target.aname.c_str()[0] == 'n')
                rap_target.resolve_special_atom(p, rap_index.loc());
            if (!rap_target.resno) return;
            #if _dbg_acvm_apply
            cout << "Resolved index " << rap_target.resno << endl;
            #endif
        }

        Point pt_fulcrum = rap_fulcrum.loc(), pt_index = rap_index.loc(), pt_target = rap_target.loc();
        if (tgtligand && !ligand)
        {
            cerr << "Called ligand-target active motion without a ligand set." << endl;
            throw 0xbadc0de;                // bad code monkey no banana!
        }
        if (tgtligand)
        {
            Atom* liga = ligand->get_nearest_atom_to_line(rap_start.loc(), rap_end.loc());
            if (!liga) return;
            pt_target = liga->loc;
            #if _dbg_acvm_apply
            cout << "Target is ligand." << endl;
            #endif
        }
        Atom* closest_to_ligand = nullptr;
        if (entire)
        {
            closest_to_ligand = p->get_nearest_atom(pt_target, rap_start.resno, rap_end.resno);
            pt_index = closest_to_ligand->loc;
            #if _dbg_acvm_apply
            cout << "Index is entire." << endl;
            #endif
        }
        Vector tolerance = pt_index.subtract(pt_target);
        Rotation rot;
        if (fixclash)
        {
            #if _dbg_acvm_apply
            // cout << "Motion fixes clash." << endl;
            #endif
            int i;
            Molecule *aai = p->get_residue(rap_index.resno), *aat = p->get_residue(rap_target.resno);
            if (entire) aai = (Molecule*)closest_to_ligand->mol;
            if (tgtligand) aat = ligand;
            if (!aai || !aat)
            {
                cerr << "Something went wrong." << endl;
                throw 0xbadc0de;
            }
            rot = align_points_3d(aat->get_barycenter(), aai->get_barycenter(), pt_fulcrum);

            cout << "Rotating " << rap_start.resno << "->" << rap_end.resno << " about " << rap_fulcrum.resno 
                << " to avoid clash..." << flush;
            float oldclash = aai->get_intermol_clashes(aat), step = 1.5*fiftyseventh;
            for (i=0; i<60; i++)
            {
                p->rotate_piece(rap_start.resno, rap_end.resno, pt_fulcrum, rot.v, step);
                cout << ".";
                float clash = aai->get_intermol_clashes(aat);
                // if (!i && clash > oldclash) step *= -1;
                if (clash <= clash_limit_per_aa) break;
            }
            cout << endl;
        }
        else if (morethan)
        {
            tolerance.r = tgtdist;
            pt_target = pt_target.add(tolerance);
            rot = align_points_3d(pt_index, pt_target, pt_fulcrum);

            cout << "Rotating " << rap_start.resno << "->" << rap_end.resno << " " 
                << (rot.a*fiftyseven) << " deg about " << rap_fulcrum.resno << " to minimum residue distance " << tolerance.r << "A..." << endl;
            p->rotate_piece(rap_start.resno, rap_end.resno, pt_fulcrum, rot.v, rot.a);
        }
        else
        {
            tolerance.r = tgtdist;
            pt_target = pt_target.add(tolerance);
            rot = align_points_3d(pt_index, pt_target, pt_fulcrum);

            cout << "Rotating " << rap_start.resno << "->" << rap_end.resno << " " 
                << (rot.a*fiftyseven) << " deg about " << rap_fulcrum.resno << " to maximum contact distance " << tolerance.r << "A..." << endl;
            p->rotate_piece(rap_start.resno, rap_end.resno, pt_fulcrum, rot.v, rot.a);
        }
    }
    else if (acvmt == acvm_prox)
    {
        rap_index.resolve_resno(p);
        if (!rap_index.resno) return;
        rap_target.resolve_resno(p);
        if (!rap_target.resno) return;

        if (rap_index.aname.c_str()[0] == 'n')
            if (!rap_index.resolve_special_atom(p, rap_target.loc())) return;
        if (rap_target.aname.c_str()[0] == 'n')
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
            a1 = aa1->get_atom(rap_index.aname.c_str());
            a2 = aa2->get_atom(rap_target.aname.c_str());

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
}
