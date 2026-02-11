
#include "soft.h"

float soft_contact_elasticity = initial_soft_contact_elasticity;

int SoftRegion::num_contacts()
{
    int i;
    for (i=0; i<allocated && contacts[i].local && contacts[i].distant; i++);
    return i;
}

AminoAcid *SoftRegion::get_local_contact(int i, Protein *p)
{
    return p->get_residue(contacts[i].local);
}

AminoAcid *SoftRegion::get_distant_contact(int i, Protein *p)
{
    return p->get_residue(contacts[i].distant);
}

float SoftRegion::get_contact_original_distance(int i)
{
    return contacts[i].CA_distance;
}

Atom* SoftRegion::get_pivot_atom_by_contact_idx(int i, Protein *p)
{
    if (i >= allocated) return nullptr;
    int n = num_contacts();
    if (i >= n) return nullptr;
    Atom *a, *b;
    AminoAcid *aa1 = p->get_residue(contacts[i].local), *aa2 = p->get_residue(contacts[i].distant);
    if (!aa1 || !aa2) return nullptr;
    aa1->mutual_closest_atoms(aa2, &a, &b);
    return a;
}

void SoftRegion::optimize_contact(Protein *p, int cidx)
{
    int i;
    if (i >= allocated) return;
    if (contacts[i].local && contacts[i].distant)
    {
        AminoAcid *aa1 = p->get_residue(contacts[i].local), *aa2 = p->get_residue(contacts[i].distant);
        if (!aa1 || !aa2) return;
        Molecule* mols[3];
        mols[0] = (Molecule*)aa1;
        mols[1] = (Molecule*)aa2;
        mols[2] = nullptr;
        MovabilityType mt1 = aa1->movability, mt2 = aa2->movability;
        if (!((aa1->movability & MOV_PINNED) || (aa2->movability & MOV_PINNED)))
        {
            aa1->movability = aa2->movability = MOV_FLEXONLY;
            Molecule::conform_molecules(mols, 200);
        }
        aa1->movability = mt1;
        aa2->movability = mt2;
    }
}

float SoftRegion::contact_anomaly(Protein *p, int cidx, bool ip)
{
    int i;
    float anomaly = 0;
    if (cidx >= allocated) cidx = -1;

    int ifrom = (cidx >= 0) ? cidx : 0, ito = (cidx >= 0) ? cidx : allocated+1;
    for (i=ifrom; i<=ito; i++)
    {
        if (contacts[i].local && contacts[i].distant)
        {
            if (ip && contacts[i].paired) continue;
            AminoAcid *aa1 = p->get_residue(contacts[i].local), *aa2 = p->get_residue(contacts[i].distant);
            if (!aa1 || !aa2) continue;
            Interaction e = aa1->get_intermol_binding(aa2);
            Atom *lca, *lcb;
            aa1->Molecule::mutual_closest_atoms(aa2, &lca, &lcb);
            if (lca && lcb && lca->distance_to(lcb) > (lca->vdW_radius+lcb->vdW_radius)) e.clash = 0;
            float f = e.summed();
            #if optimize_energies_in_contact_anomaly_check
            if (f < (1.0-contact_energy_allowance_for_optimization)*contacts[i].energy)
            {
                Molecule* mols[3];
                mols[0] = (Molecule*)aa1;
                mols[1] = (Molecule*)aa2;
                mols[2] = nullptr;
                Molecule::conform_molecules(mols, contact_energy_optimization_iterations);
                e = aa1->get_intermol_binding(aa2).summed();
                f = e.summed();
            }
            #endif
            if (f > contacts[i].energy) anomaly += (f - contacts[i].energy);
            /* cout << "Binding energy between " << aa1->get_name() << " and " << aa2->get_name()
                << " = " << f << " originally " << contacts[i].energy << endl; */
        }
        else break;
    }

    return anomaly;
}

float SoftRegion::contact_distance_anomaly(Protein *p, int cidx, bool ip)
{
    int i;
    float anomaly = 0;
    if (cidx >= allocated) cidx = -1;

    for (i=0; i<allocated; i++)
    {
        if (cidx >= 0) i = cidx;
        if (contacts[i].local && contacts[i].distant)
        {
            if (ip && contacts[i].paired) continue;
            AminoAcid *aa1 = p->get_residue(contacts[i].local), *aa2 = p->get_residue(contacts[i].distant);
            if (!aa1 || !aa2) continue;

            float r = aa1->get_CA_location().get_3d_distance(aa2->get_CA_location());
            if (r > contacts[i].CA_distance) anomaly += (r - contacts[i].CA_distance);
        }
        else break;
        if (cidx >= 0) break;
    }
    return anomaly;
}

void SoftRegion::link_region(SoftRegion* prev)
{
    prev_rgn_end = prev->rgn.end;
    prev->next_rgn_start = rgn.start;
}

SoftRegion::SoftRegion()
{
    ;
}

SoftRegion::~SoftRegion()
{
    if (contacts) delete[] contacts;
    contacts = nullptr;
}

SoftRegion::SoftRegion(SoftRegion&& sr) noexcept
{
    allocated = sr.allocated;
    contacts = sr.contacts;
    sr.contacts = nullptr;
}

SoftRegion& SoftRegion::operator=(SoftRegion&& sr) noexcept
{
    if (this != &sr)
    {
        allocated = sr.allocated;
        contacts = sr.contacts;
        sr.contacts = nullptr;
    }
    return *this;
}

SoftRegion& SoftRegion::operator=(const SoftRegion& sr)
{
    if (this != &sr)
    {
        contacts = new SoftContact[allocated = sr.allocated];
        int i;
        for (i=0; i<allocated; i++)
        {
            contacts[i] = sr.contacts[i];
        }
    }
    return *this;
}

SoftRegion::SoftRegion(const SoftRegion& sr)
{
    contacts = new SoftContact[allocated = sr.allocated];
    int i;
    for (i=0; i<allocated; i++)
    {
        contacts[i] = sr.contacts[i];
    }
}

bool SoftRegion::check_chain_constraints(Protein* prot)
{
    AminoAcid *aapre=nullptr, *aanrs=nullptr, *saa, *eaa;

    prev_rgn_violated = next_rgn_violated = false;
    saa = prot->get_residue(rgn.start);
    eaa = prot->get_residue(rgn.end);
    if (!saa || !eaa) return true;

    if (prev_rgn_end) aapre = prot->get_residue(prev_rgn_end);
    if (next_rgn_start) aanrs = prot->get_residue(next_rgn_start);

    if (aapre)
    {
        int Dn = saa->get_residue_no() - aapre->get_residue_no();
        prev_rgn_dist = saa->get_CA_location().get_3d_distance(aapre->get_CA_location());
        if (prev_rgn_dist > chain_length_per_aa * Dn) prev_rgn_violated = true;
    }
    if (aanrs)
    {
        int Dn = aanrs->get_residue_no() - eaa->get_residue_no();
        next_rgn_dist = aanrs->get_CA_location().get_3d_distance(eaa->get_CA_location());
        if (next_rgn_dist > chain_length_per_aa * Dn) next_rgn_violated = true;
    }
    return !prev_rgn_violated && !next_rgn_violated;
}

void SoftRegion::add_contact(int local, int distant, Protein* p, bool paired)
{
    int i;

    if (!allocated)
    {
        allocated = soft_contact_alloc_block;
        contacts = new SoftContact[allocated];
        i=0;
    }
    else
    {
        for (i=0; contacts[i].local; i++)               // get count
        {
            if (contacts[i].local == local && contacts[i].distant == distant) return;               // prevent duplicates
            if (contacts[i].local == distant && contacts[i].distant == local) return;
        }
        if (i >= allocated-1)
        {
            int j;
            SoftContact* tmp = new SoftContact[allocated+soft_contact_alloc_block];
            for (j=0; j<allocated; j++) tmp[j] = contacts[j];
            tmp[j].local = tmp[j].distant = 0;
            delete[] contacts;
            contacts = tmp;
            allocated += soft_contact_alloc_block;
        }
    }

    if (!local || !distant) return;

    contacts[i].local = local;
    contacts[i].distant = distant;
    AminoAcid *aa1 = p->get_residue(local), *aa2 = p->get_residue(distant);
    Interaction e;
    if (aa1 && aa2)
    {
        e = aa1->get_intermol_binding(aa2);
        contacts[i].CA_distance = aa1->get_CA_location().get_3d_distance(aa2->get_CA_location());
    }
    contacts[i].energy = e.summed();
    contacts[i].paired = paired;
    // cout << contacts[i].local->get_name() << "..." << contacts[i].distant->get_name() << " with energy " << contacts[i].energy << endl;
    i++;
    contacts[i].local = contacts[i].distant = 0;
}

void soft_docking_iteration(Protein *protein, Molecule* ligand, int nsoftrgn, SoftRegion* softrgns, float softness)
{
    int i, j, l;

    for (i=0; i<nsoftrgn; i++)
    {
        int srnc = softrgns[i].num_contacts();
        Point softpush(0,0,0);
        Point foravg[softrgns[i].rgn.end-softrgns[i].rgn.start+16];
        float pushmax = 0;
        l = 0;
        float near_pivot = 0, near_terminus = 0;

        Atom* piva;
        int srci = srnc ? (rand() % srnc) : 0;             // soft region contact index.
        piva = srnc ? softrgns[i].get_pivot_atom_by_contact_idx(srci, protein) : protein->region_pivot_atom(softrgns[i].rgn);
        if (!piva) continue;
        if (piva->vanished) continue;

        Point C = piva->loc;
        Point S = softrgns[i].rgn.start_CA_location(protein);
        Point E = softrgns[i].rgn.end_CA_location(protein);
        float rotation_share = 0, translation_share = 0;

        if (srnc)
        {
            softrgns[i].get_local_contact(srci, protein)->been_flexed = true;
            softrgns[i].get_distant_contact(srci, protein)->been_flexed = true;
        }
        else
        {
            AminoAcid* aa = protein->get_residue(piva->residue);
            if (aa) aa->been_flexed = true;
        }

        for (j=softrgns[i].rgn.start; j<=softrgns[i].rgn.end; j++)
        {
            AminoAcid* aa = protein->get_residue(j);
            if (!aa) continue;

            aa->been_flexed = true;

            Point A = aa->get_CA_location();
            float aa_pivot_dist = A.get_3d_distance(C);
            float aa_term_dist = fmin(A.get_3d_distance(S), A.get_3d_distance(E));
            float pivter = aa_pivot_dist + aa_term_dist;

            float c = aa->get_intermol_clashes(ligand);
            if (c > clash_limit_per_aa*10)
            {
                Point A = aa->get_CA_location();
                foravg[l++] = A;
                Vector AB = A.subtract(ligand->get_nearest_atom(A)->loc);
                AB.r = pow(c, 1.0/3);
                if (!isinf(AB.r) && !isnan(AB.r))
                {
                    softpush = softpush.add(AB);
                    if (AB.r > pushmax) pushmax = AB.r;
                    // Clash will have more of an effect on translation or rotation depending on its distance from the pivot.
                    rotation_share += (aa_pivot_dist / pivter) * AB.r;
                    translation_share += (aa_term_dist / pivter) * AB.r;
                }
            }
        }

        float pushmagn = softpush.magnitude();
        float pullmagn = 0;
        softpush.scale(pushmax);
        for (j=0; j<srnc; j++)
        {
            AminoAcid *aa1 = softrgns[i].get_local_contact(j, protein), *aa2 = softrgns[i].get_distant_contact(j, protein);
            if (!aa1 || !aa2) continue;
            Atom *CA1 = aa1->get_atom("CA"), *CA2 = aa2->get_atom("CA");
            float ra = CA1->distance_to(CA2) - softrgns[i].get_contact_original_distance(j);
            Vector ctpull = CA2->loc.subtract(CA1->loc);
            ctpull.r = ra * soft_contact_elasticity * frand(0,1);
            softpush = softpush.subtract(ctpull);
            pullmagn += ctpull.r;
        }
        pushmax = softpush.magnitude();

        pushmax *= soft_push_multiplier;
        if (pushmax > speed_limit) pushmax = speed_limit;

        if (l)
        {
            // Figure how much to attempt rotation and how much to attempt translation.
            float roxl = rotation_share + translation_share;
            rotation_share /= roxl;
            translation_share /= roxl;

            Vector AB = softpush;
            AB.r = pushmax;
            Vector ABr = AB, ABx = AB;
            ABr.r *= rotation_share;
            ABx.r *= translation_share;
            ABr.r = fmin(fmax(0, ABr.r), speed_limit);
            ABx.r = fmin(fmax(0, ABx.r), speed_limit);

            // Do the translation stepwise and stop when contacts anomaly exceeds threshold.
            float translation_accomplished;
            float translation_step = 0.05;
            Vector ABx_step = ABx;
            ABx_step.r *= translation_step;
            float cwaybefore = protein->get_internal_clashes(softrgns[i].rgn.start, softrgns[i].rgn.end, repack_soft_clashes, soft_repack_iterations)
                + protein->get_intermol_clashes(ligand) + softrgns[i].contact_distance_anomaly(protein);
            float cbefore, cafter, clbefore, clafter;
            for (translation_accomplished = 0; translation_accomplished < 1; translation_accomplished += translation_step)
            {
                // Get energy before performing soft motion
                clbefore = protein->get_intermol_clashes(ligand);
                cbefore = protein->get_internal_clashes(softrgns[i].rgn.start, softrgns[i].rgn.end, repack_soft_clashes, soft_repack_iterations)
                    + clbefore + softrgns[i].contact_distance_anomaly(protein, -1, false);

                #if move_ligand_with_soft_motion
                Pose ligand_was(ligand);
                #endif

                // Perform soft motion
                #if _dbg_soft_motions
                int dbgi, dbgn = protein->get_end_resno();
                Point dbgaaloc[dbgn+2];
                for (dbgi = 1; dbgi <= dbgn; dbgi++)
                {
                    AminoAcid *dbgaa = protein->get_residue(dbgi);
                    if (dbgaa) dbgaaloc[dbgi] = dbgaa->get_CA_location();
                }
                #endif
                protein->move_piece(softrgns[i].rgn.start, softrgns[i].rgn.end, ABx_step);
                #if _dbg_soft_motions
                cout << "Moved " << softrgns[i].rgn.start << "-" << softrgns[i].rgn.end << endl;
                for (dbgi = 1; dbgi <= dbgn; dbgi++)
                {
                    AminoAcid *dbgaa = protein->get_residue(dbgi);
                    if (dbgaa)
                    {
                        float dbgmovd = dbgaaloc[dbgi].get_3d_distance(dbgaa->get_CA_location());
                        if (dbgmovd) cout << dbgaa->get_name() << " moved " << dbgaaloc[dbgi].get_3d_distance(dbgaa->get_CA_location()) << "A." << endl;
                    }
                }
                cout << endl;
                #endif

                // If the ligand is "staying" near any part of the soft region, move it too.
                #if move_ligand_with_soft_motion
                if (ligand->stay_close_other 
                    && !ligand->glued_to_mol()                                      // glued ligand will move with side chain already
                    && !ligand->stay_close_other->vanished
                    && ligand->stay_close_other->residue >= softrgns[i].rgn.start
                    && ligand->stay_close_other->residue <= softrgns[i].rgn.end
                   )
                   ligand->move(ABx_step);
                #endif

                // Get energy after performing soft motion
                clafter = protein->get_intermol_clashes(ligand);
                cafter = protein->get_internal_clashes(softrgns[i].rgn.start, softrgns[i].rgn.end, repack_soft_clashes, soft_repack_iterations)
                    + clafter + softrgns[i].contact_distance_anomaly(protein, -1, false);

                #if move_ligand_with_soft_motion
                if (clafter > (1.0-contact_energy_allowance_for_optimization)*clbefore)
                {
                    if (!ligand->glued_to_mol()) ligand_was.restore_state(ligand);
                    clafter = protein->get_intermol_clashes(ligand);
                    cafter = protein->get_internal_clashes(softrgns[i].rgn.start, softrgns[i].rgn.end, repack_soft_clashes, soft_repack_iterations)
                        + clafter + softrgns[i].contact_distance_anomaly(protein, -1, false);
                }
                #endif

                if (cafter > (1.0-contact_energy_allowance_for_optimization)*cbefore 
                    || (pushmagn > pullmagn && !softrgns[i].check_chain_constraints(protein)))
                {
                    protein->undo();
                    #if move_ligand_with_soft_motion
                    if (!ligand->glued_to_mol()) ligand_was.restore_state(ligand);
                    #endif
                    cafter = cbefore;
                    clafter = clbefore;
                    break;
                }
                #if move_ligand_with_soft_motion
                else if (!ligand->glued_to_mol()) ligand_was.copy_state(ligand);
                #endif
            }
            /*cout << "Moved residues " << softrgns[i].rgn.start << "-" << softrgns[i].rgn.end << " "
                << translation_accomplished*ABx.r << "A." << endl << endl;*/

            rotation_share += translation_share * (1.0 - translation_accomplished);
            translation_share *= translation_accomplished;
            ABx.r *= translation_accomplished;
            Point CX = C.add(ABx);
            if (audit) fprintf(audit, "Accepted soft dock translation from %g to %g.\n", cwaybefore, cafter);

            Point A = average_of_points(foravg, l);
            Point B = A.add(ABr);
            Rotation rot = align_points_3d(A, B, CX);
            rot.a = fmin(rot.a, 0.1*softness*fiftyseventh);

            cbefore = protein->get_internal_clashes(softrgns[i].rgn.start, softrgns[i].rgn.end, repack_soft_clashes, soft_repack_iterations)
                + protein->get_intermol_clashes(ligand) + softrgns[i].contact_distance_anomaly(protein, -1, false);
            protein->rotate_piece(softrgns[i].rgn.start, softrgns[i].rgn.end, CX, rot.v, rot.a);
            cafter = protein->get_internal_clashes(softrgns[i].rgn.start, softrgns[i].rgn.end, repack_soft_clashes, soft_repack_iterations)
                + protein->get_intermol_clashes(ligand) + softrgns[i].contact_distance_anomaly(protein, -1, false);
            if (cafter > (1.0-contact_energy_allowance_for_optimization)*cbefore || !softrgns[i].check_chain_constraints(protein))
                protein->undo();
            else if (audit) fprintf(audit, "Accepted soft dock rotation from %g to %g.\n", cbefore, cafter);
        }
    }

    #if progressively_increase_elasticity
    soft_contact_elasticity = 1.0 - ((1.0 - soft_contact_elasticity)*(1.0-soft_contact_elasticity_progress));
    #endif
}
