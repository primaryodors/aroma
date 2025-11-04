
#include "scoring.h"

float init_total_binding_by_type[_INTER_TYPES_LIMIT];
float fin_total_binding_by_type[_INTER_TYPES_LIMIT];

float* initial_binding;
#if compute_vdw_repulsion
float* initial_vdWrepl;
#endif

DockResult::DockResult()
{
    ;
}

DockResult::DockResult(Protein* protein, Molecule* ligand, Point search_size, int* addl_resno, int drcount, Molecule** waters, bool is_movie)
{
    int end1 = SPHREACH_MAX+4;
    AminoAcid* reaches_spheroid[end1];
    int sphres = protein->get_residues_can_clash_ligand(reaches_spheroid, ligand, ligand->get_barycenter(), search_size, addl_resno);
    // cout << "sphres " << sphres << endl;
    mprot = protein;
    mlig = ligand;
    Molecule* met = protein->metals_as_molecule();
    int i, li, j, k, l, n;

    char metrics[end1][20];
    float lmkJmol[end1];
    float lmstab[end1];
    float limkJmol[end1];
    #if compute_vdw_repulsion
    float lmvdWrepl[end1];
    float limvdWrepl[end1];
    #endif
    int nmetrics = 0;
    float btot = 0;
    float btots = 0;
    float pstot = 0;
    char* lmba1n[end1];
    char* lmba2n[end1];
    char* lmca1n[end1];
    char* lmca2n[end1];
    float lmc[end1];
    bool lis_mcr[end1];

    movie_mode = is_movie;
    if (!movie_mode)
    {
        ligpos.copy_state(ligand);
        ligvol = ligand->get_volume();
        stay_close_ligand = ligand->stay_close_mine;
        stay_close_protein = ligand->stay_close_other;
        ligand_solvation_energy = ligand->solvent_free_energy();
        ligand_pocket_wet_energy = ligand_solvation_energy;
        ligand_waters_energy = 0;
        if (!global_water.get_atom_count())
        {
            global_water.from_smiles("O{O1}");
            global_water.obtain_vdW_surface(vdw_surface_density);
        }
        Atom* polarH = global_water.get_atom(1);
        Atom* polarO = global_water.get_atom(0);

        #if compute_lsrb
        const Point* ligsurf = ligand->obtain_vdW_surface(vdw_surface_density);
        Atom** ligva = ligand->get_vdW_vertex_atoms();
        int nlsv = ligand->get_vdW_vertex_count();
        lsrb_points = new Point[nlsv+4];
        float lsrb_rawval[nlsv+4];
        for (i=0; i<nlsv; i++)
        {
            lsrb_rawval[i] = 0;
            lsrb_points[i] = ligsurf[i];
            lsrb_points[i].weight = 0;
        }
        for (i=0; i<nlsv; i++)
        {
            float polarity = fabs(ligva[i]->is_polar());
            Atom* a = protein->get_nearest_atom(ligsurf[i]);
            if (a)
            {
                Point pt = ligsurf[i];
                Point aloc = a->loc;
                float r = pt.get_3d_distance(aloc);
                float vdw = a->vdW_radius;

                float lsrbpart = 0;
                if (r < vdw*lsrb_vdw_multiplier)
                {
                    lsrbpart = 1;
                    #if _dbg_lsrb
                    cout << a->aaletter << a->residue << ":" << a->name << " at " << r << "Å, lsrb +" << lsrbpart << endl;
                    #endif
                }
                else if (r < vdw*lsrb_vdw_multiplier*2)
                {
                    r -= vdw*lsrb_vdw_multiplier;
                    r /= (vdw*lsrb_vdw_multiplier);
                    lsrbpart = 1.0-r;
                    #if _dbg_lsrb
                    cout << a->aaletter << a->residue << ":" << a->name << " at " << r << "Å, lsrb +" << lsrbpart << endl;
                    #endif
                }

                float theta = find_3d_angle(aloc, pt, ligva[i]->loc);
                float thcos = cos(theta);
                if (thcos < 0) continue;
                lsrbpart *= thcos;

                if (lsrbpart)
                {
                    pt.weight = lsrbpart;
                    lsrb_rawval[i] = lsrbpart;
                    lsrb_points[i] = pt;
                }
            }
        }
        for (i=0; i<nlsv; i++)
        {
            Point pt = ligsurf[i];
            for (j=0; j<nlsv; j++)
            {
                if (j==i) continue;
                float r = pt.get_3d_distance(ligsurf[j]);
                if (r > 0.5) continue;
                r = r*4+1;
                float f = 1.0 - lsrb_rawval[j] / (r*r);
                lsrb_points[i].weight = 1.0 - ((1.0 - lsrb_points[i].weight)*f);
            }
            lsrb_points[i].weight = sqrt(lsrb_points[i].weight);
            ligand_surface_receptor_binding += lsrb_points[i].weight;
        }
        if (nlsv)
        {
            ligand_surface_receptor_binding /= nlsv;

            #if _dbg_lsrb
            cout << endl << endl;
            #endif
        }
        nlsrb_points = nlsv;
        #endif

        if (waters) ligand_waters_energy += ligand->get_intermol_binding(waters).summed();
    }

    for (i=0; i<end1; i++)
    {
        #if compute_vdw_repulsion
        lmkJmol[i] = lmstab[i] = limkJmol[i] = lmvdWrepl[i] = limvdWrepl[i] = lmc[i] = 0;
        #else
        lmkJmol[i] = lmstab[i] = limkJmol[i] = lmc[i] = 0;
        #endif
    }

    worst_energy = worst_nrg_aa = 0;
    worst_clash_1 = worst_clash_2 = nullptr;

    // if (debug) *debug << "Pose " << pose << " pathnode " << nodeno /*<< " clashes " << clash*/ << endl;

    ligand->clear_atom_binding_energies();

    float final_binding[end1];
    #if compute_vdw_repulsion
    float final_vdWrepl[end1];
    for (i=0; i<end1; i++) final_binding[i] = final_vdWrepl[i] = 0;
    #else
    for (i=0; i<end1; i++) final_binding[i] = 0;
    #endif

    AminoAcid* allres[protein->get_end_resno()+16]; 
    int qpr = protein->fetch_residues_near(ligand->get_barycenter(), 100000, allres, false);
    Molecule* postaa[qpr];
    postaa[0] = ligand;
    for (i=0; i<qpr; i++)
    {
        postaa[i+1] = reinterpret_cast<Molecule*>(allres[i]);
    }

    for (i=0; i<_INTER_TYPES_LIMIT; i++) fin_total_binding_by_type[i] = total_binding_by_type[i];

    #if _peratom_audit
    interaudit.clear();
    interauditing = true;
    #endif

    for (i=0; i<_INTER_TYPES_LIMIT; i++)
    {
        init_total_binding_by_type[i] = 0;
        fin_total_binding_by_type[i] = 0;
        total_binding_by_type[i] = 0;
    }

    int prot_seq_len = protein->get_end_resno();
    #if compute_clashdirs
    residue_clash = new float[prot_seq_len+8];
    res_clash_dir = new Vector[prot_seq_len+8];

    for (i=0; i<=prot_seq_len; i++)
    {
        residue_clash[i] = 0;
        res_clash_dir[i] = Vector(0,0,0);
    }
    #endif

    int mtlmetric = -1;
    mcoord_charge_repulsion = 0;
    ninterall = 0;
    atomlvl = "";
    reaches_spheroid[sphres] = nullptr;
    for (i=0; i<sphres; i++)
    {
        if (!reaches_spheroid[i]) continue;
        if (!protein->aa_ptr_in_range(reaches_spheroid[i])) continue;
        reaches_spheroid[i]->clear_atom_binding_energies();
        int resno = reaches_spheroid[i]->get_residue_no();

        if ((reaches_spheroid[i]->has_hbond_acceptors() && ligand->has_hbond_donors())
            || (reaches_spheroid[i]->has_hbond_donors() && ligand->has_hbond_acceptors())
           )
        {
            int aai, aan = reaches_spheroid[i]->get_atom_count(), ln = ligand->get_atom_count();
            for (li=0; li<ln; li++)
            {
                Atom* a = ligand->get_atom(li);
                if (!a) continue;
                int aZ = a->Z;

                #if make_thiolates_in_scoring
                // this stupid thing is thiolating here too
                if (aZ > 10 && a->get_family() == CHALCOGEN)
                {
                    if (met && !a->get_charge() && met->get_atom_count())
                    {
                        Atom* mtl = met->get_nearest_atom(a->loc);
                        if (mtl && mtl->is_metal() && mtl->distance_to(a) < _INTERA_R_CUTOFF)
                        {
                            Atom* H = a->is_bonded_to("H");
                            if (H)
                            {
                                ligand->delete_atom(H);
                                a->increment_charge(-1);
                            }
                        }
                    }
                }
                #endif

                float apol = a->is_polar();
                if (fabs(apol) < hydrophilicity_cutoff) continue;

                if (sgn(apol) > 0 && aZ > 1) continue;
                for (aai=0; aai<aan; aai++)
                {
                    Atom* b = reaches_spheroid[i]->get_atom(aai);
                    if (!b) continue;
                    float bpol = b->is_polar();
                    if (fabs(bpol) < hydrophilicity_cutoff) continue;
                    int bZ = b->Z;
                    if (sgn(bpol) > 0 && bZ > 1) continue;
                    if (sgn(apol) == sgn(bpol)) continue;
                    float r = a->distance_to(b);
                    if (r > 6) continue;

                    Atom* H = (aZ==1) ? a : b;
                    Atom* target = (aZ==1) ? b : a;
                    Atom* targetH = target->is_bonded_to("H");
                    Bond* bond = H->get_bond_by_idx(0);
                    if (!bond) continue;
                    Atom* heavy = bond->atom2;
                    if (!heavy) continue;
                    if (heavy->get_family() != CHALCOGEN) continue;
                    if (heavy->is_pi()) continue;

                    bond = heavy->get_bond_by_idx(0);
                    if (!bond || !bond->atom2 || !bond->can_rotate) continue;

                    Vector axis = heavy->loc.subtract(bond->atom2->loc);
                    float theta = find_angle_along_vector(H->loc, target->loc, heavy->loc, axis);
                    if (targetH && heavy->distance_to(targetH) < target->distance_to(H)) theta += M_PI;

                    Point newHloc1 = rotate3D(H->loc, heavy->loc, axis, theta);
                    Point newHloc2 = rotate3D(H->loc, heavy->loc, axis, -theta);

                    if (newHloc1.get_3d_distance(target->loc) < newHloc2.get_3d_distance(target->loc))
                        H->move(newHloc1);
                    else H->move(newHloc2);
                }
            }
        }

        lis_mcr[nmetrics] = false;
        mc_bpotential = 0;
        compute_interall = true;
        Interaction lb = ligand->get_intermol_binding(reaches_spheroid[i], false);
        compute_interall = false;
        float clash = lb.clash;
        if (ligand->get_worst_clash() > worst_energy)
        {
            worst_energy = ligand->get_worst_clash();
            worst_clash_1 = ligand->clash1;
            worst_clash_2 = ligand->clash2;
        }
        if (/* lb.summed() > 0 && */ clash > worst_nrg_aa)
        {
            worst_nrg_aa = clash;
            #if _dbg_worst_energy
            cout << reaches_spheroid[i]->get_name() << " binding strength " << lb << " clashes " << clash << " updates worst_nrg_aa to " << worst_nrg_aa << endl;
            #endif
        }

        if (fabs(lb.summed()) >= 0.1)
        {
            float sfe = reaches_spheroid[i]->sc_hfe();
            float occlusion_by_ligand = reaches_spheroid[i]->surface_occlusion(ligand);
            float occlusion_by_neighbors = reaches_spheroid[i]->surface_occlusion((Molecule**)reaches_spheroid);
            sfe *= fmax(0, 1.0 - occlusion_by_neighbors);
            pocket_wet_solvation_energy += sfe;
            float sbe = sfe * fmax(0, 1.0 - occlusion_by_ligand);
            #if _dbg_pocket_DeltaG_solv
            cout << reaches_spheroid[i]->get_name() << " sfe = " << sfe << " sbe = " << sbe << endl;
            #endif
            pocket_bound_solvation_energy += sbe;
            if (sfe < 0 && reaches_spheroid[i]->is_ic_res) pocket_ic_DeltaG_solvation += (sfe-sbe);     // if sfe is more negative than sbe, ligand stabilizes contact by reducing its solvation effect.
        }

        #if compute_clashdirs
        if (lb > 0 && ligand->clash1 && ligand->clash2)
        {
            residue_clash[resno] += lb;
            Vector clashdir = ligand->clash2->loc.subtract(ligand->clash1->loc);
            clashdir.r = ligand->clash1->vdW_radius + ligand->clash2->vdW_radius - ligand->clash1->distance_to(ligand->clash2);
            res_clash_dir[resno] = res_clash_dir[resno].add(clashdir);
        }
        #endif

        #if include_residue_eclipses
        lb -= fmax(reaches_spheroid[i]->total_eclipses() - reaches_spheroid[i]->initial_eclipses, 0);
        #endif

        #if _dbg_51e2_ionic
        if (resno == 262)
        {
            cout << endl << resno << " charge " << reaches_spheroid[i]->get_charge()
                << " vs. ligand charge " << ligand->get_charge()
                << ": " << lb << endl << endl;
        }
        #endif

        if (lb.summed() < -500) lb = 0;
        lmkJmol[nmetrics] = lb.summed();
        lmstab[nmetrics] = reaches_spheroid[i]->estimated_stability((Molecule**)reaches_spheroid);
        lmc[nmetrics] = -mc_bpotential / missed_connection.r;

        lmba1n[nmetrics] = ligand->best_intera ? ligand->best_intera->name : nullptr;
        lmba2n[nmetrics] = ligand->best_other_intera ? ligand->best_other_intera->name : nullptr;
        lmca1n[nmetrics] = ligand->clash1 ? ligand->clash1->name : nullptr;
        lmca2n[nmetrics] = ligand->clash2 ? ligand->clash2->name : nullptr;

        BallesterosWeinstein bw = protein->get_bw_from_resno(resno);
        if (bw.helix_no)
            sprintf(metrics[nmetrics], "%s%d(%d.%d)", reaches_spheroid[i]->get_3letter(), resno, bw.helix_no, bw.member_no);
        else
            sprintf(metrics[nmetrics], "%s%d", reaches_spheroid[i]->get_3letter(), resno);
        // cout << metrics[metcount] << ": " << lb << " . ";

        #if compute_vdw_repulsion
        lmvdWrepl[metcount] = 0;
        lmvdWrepl[metcount] += ligand->get_vdW_repulsion(reaches_spheroid[i]);
        /*for (j=0; j<sphres; j++)
        {
            if (j == i) continue;
            mvdWrepl[metcount] += reaches_spheroid[i]->get_vdW_repulsion(reaches_spheroid[j]);
        }*/
        limvdWrepl[metcount] = 0;
        #endif
        limkJmol[nmetrics] = 0;
        if (reaches_spheroid[i]->coordmtl) lis_mcr[nmetrics] = true;

        nmetrics++;
        btot += lb.summed();
        btots += lb.summed() * reaches_spheroid[i]->estimated_stability((Molecule**)reaches_spheroid);
        // cout << *(reaches_spheroid[i]) << " adds " << lb << " to btot, making " << btot << endl;

        if (reaches_spheroid[i]->coordmtl && reaches_spheroid[i]->coordmtl->residue == reaches_spheroid[i]->get_residue_no())
        {
            mtlmetric = nmetrics;
            int mi;

            for (mi=0; mi < nmetrics; mi++)
            {
                if (!strcmp(metrics[mi], "Metals"))
                {
                    mtlmetric = mi;
                    break;
                }
            }

            strcpy(metrics[mtlmetric], "Metals");
            lmba1n[mtlmetric] = nullptr;
            lmba2n[mtlmetric] = nullptr;
            lmca1n[mtlmetric] = nullptr;
            lmca2n[mtlmetric] = nullptr;
            float lmb = (mtlmetric >= nmetrics) ? 0 : lmkJmol[mtlmetric];

            Molecule lm("MTL");
            lm.add_existing_atom(reaches_spheroid[i]->coordmtl);
            float lmf = ligand->get_intermol_binding(&lm).summed();
            btot += lmf;
            btots += lmf;
            lmb += lmf;

            lmkJmol[mtlmetric] = lmb;
            if (mtlmetric >= nmetrics) nmetrics++;
        }

        float lf = ligand->get_intermol_polar_sat(reaches_spheroid[i]);
        pstot -= lf;

        #if _dbg_polsat
        cout << *(reaches_spheroid[i]) << " adds " << lf << " to pstot, making " << pstot << endl;
        #endif

        float ichg;
        if (reaches_spheroid[i]->coordmtl && (ichg = reaches_spheroid[i]->get_charge()))
        {
            for (j=0; j<sphres; j++)
            {
                if (j==i) continue;
                if (reaches_spheroid[j]->coordmtl) continue;
                float jchg = reaches_spheroid[j]->get_charge();
                if (!jchg) continue;
                if (sgn(ichg) == sgn(jchg))
                {
                    float er = reaches_spheroid[i]->distance_to(reaches_spheroid[j]);
                    float eadd = 60.0 / pow(fmax(1, er / 2.0), 2);
                    mcoord_charge_repulsion += eadd;
                    #if 0
                    cout << endl << "Charge repulsion between " << reaches_spheroid[i]->get_name()
                        << " and " << reaches_spheroid[j]->get_name()
                        << ": " << eadd << " kJ/mol at " << er << "A separation." << endl;
                    #endif
                }
            }
        }
    }

    #if _dbg_pocket_DeltaG_solv
    cout << endl;
    #endif
    ligand_pocket_wet_energy = ligand->solvent_free_energy(80, 0.316, false);
    ligand_pocket_occlusion = ligand->octant_occlusion((Molecule**)reaches_spheroid);
    ligand_pocket_wet_energy *= 1.0 - fmin(1, ligand_pocket_occlusion);

    for (l=0; l<ninterall; l++)
    {
        if (fabs(interall[l]) < 0.1) continue;
        if (!interall_a1[l] || !interall_a2[l]) continue;

        if (interall_a2[l]->residue)
        {
            if (mprot)
            {
                AminoAcid* aa = mprot->get_residue(interall_a2[l]->residue);
                if (aa) atomlvl += (std::string)(aa->get_name());
                else atomlvl += std::to_string(interall_a2[l]->residue);
            }
            else atomlvl += std::to_string(interall_a2[l]->residue);
        }
        else if (interall_a2[l]->mol) atomlvl += (std::string)(((Molecule*)interall_a2[l]->mol)->get_name());
        else atomlvl += (std::string)"(ligand)";
        atomlvl += (std::string)":" + (std::string)(interall_a2[l]->name);

        atomlvl += (std::string)"~";

        if (interall_a1[l]->residue)
        {
            if (mprot)
            {
                AminoAcid* aa = mprot->get_residue(interall_a1[l]->residue);
                if (aa) atomlvl += (std::string)(aa->get_name());
                else atomlvl += std::to_string(interall_a1[l]->residue);
            }
            else atomlvl += std::to_string(interall_a1[l]->residue);
        }
        else if (interall_a1[l]->mol) atomlvl += (std::string)(((Molecule*)interall_a1[l]->mol)->get_name());
        else atomlvl += (std::string)"(ligand)";
        atomlvl += (std::string)":" + (std::string)interall_a1[l]->name;
        atomlvl += (std::string)": ";
        atomlvl += std::to_string(interall[l]) + (std::string)"\n";
    }

    #if include_mcr_in_dock_energy
    if (mtlmetric>=0)
    {
        // lmkJmol[mtlmetric] += mcoord_charge_repulsion;
        btot += mcoord_charge_repulsion;
        btots += mcoord_charge_repulsion;
    }
    // cout << btot << endl;
    #endif

    #if _peratom_audit
    cout << endl << "Interatomic Audit:" << endl;
    cout << "Total energy: " << btot << endl;
    int ian = interaudit.size(), iai;
    for (iai=0; iai<ian; iai++) cout << interaudit[iai] << endl;
    cout << endl << endl;
    interauditing = false;
    #endif

    btot += ligand_waters_energy;
    btots += ligand_waters_energy;
    if (btot < -100*ligand->get_atom_count()) btot = btots = 0;

    kJmol        = btots;
    kJmol_raw     = btot;
    ikJmol         = 0;
    polsat          = pstot;
    metric           = new char*[nmetrics+4];
    this->mkJmol      = new float[nmetrics];
    this->mstab        = new float[nmetrics];
    this->imkJmol       = new float[nmetrics];
    #if compute_vdw_repulsion
    this->mvdWrepl        = new float[metcount];
    this->imvdWrepl        = new float[metcount];
    #endif
    this->mb_atom1_name      = new const char*[nmetrics];
    this->mb_atom2_name       = new const char*[nmetrics];
    this->mc_atom1_name        = new const char*[nmetrics];
    this->mc_atom2_name         = new const char*[nmetrics];
    this->is_mcr = new bool[nmetrics];
    #if compute_missed_connections
    this->missed_connections = new float[nmetrics];
    #endif
    ligand_self = fmax(0, ligand->get_intermol_binding(ligand).summed() - ligand->get_base_clashes()
    #if include_eclipses
        + ligand->total_eclipses()
    #endif
                      );
    // cout << "Ligand base clashes " << ligand->get_base_clashes() << endl << "Ligand self " << ligand_self << endl << endl;
    A100 = protein->A100();
    if (mbbr) estimated_TDeltaS = mbbr->estimate_DeltaS() * temperature;
    // kJmol += ligand_self;
    #if _dbg_internal_energy
    cout << "Ligand internal = " << ligand_self << endl;
    #endif
    protclash = protein->get_rel_int_clashes();

    int itn;
    for (itn=0; itn<_INTER_TYPES_LIMIT; itn++)
    {
        i = itn;
        bytype[i] = total_binding_by_type[i];
        ibytype[i] = init_total_binding_by_type[i];
        ikJmol += init_total_binding_by_type[i];
        kJmol += fin_total_binding_by_type[i];
        kJmol_raw += fin_total_binding_by_type[i];
        // cout << drcount << "|" << i << " ";
    }

    // Populate the array.
    for (i=0; i<nmetrics; i++)
    {
        metric[i] = new char[max(8,(int)strlen(metrics[i])+4)];
        strcpy(metric[i], metrics[i]);
        mkJmol[i] = lmkJmol[i];
        mstab[i] = lmstab[i];
        imkJmol[i] = limkJmol[i];
        #if compute_vdw_repulsion
        mvdWrepl[i] = lmvdWrepl[i];
        imvdWrepl[i] = limvdWrepl[i];
        #endif
        is_mcr[i] = lis_mcr[i];
        mb_atom1_name[i] = lmba1n[i];
        mb_atom2_name[i] = lmba2n[i];
        mc_atom1_name[i] = lmca1n[i];
        mc_atom2_name[i] = lmca2n[i];
        #if compute_missed_connections
        missed_connections[i] = lmc[i];
        #endif
        // cout << "*" << metric[i] << ": " << mkJmol[i] << endl;
    }

    // Terminate with an empty string and a null pointer.
    metric[i] = new char[2];
    metric[i][0] = 0;
    metric[i+1] = 0;

    std::ostringstream lpdbdat;

    // Prepare a partial PDB of the ligand atoms and all involved residue sidechains.
    n = ligand->get_atom_count();
    int offset = n;
    for (l=0; l<n; l++)
    {
        Atom* a = ligand->get_atom(l);
        if (!a) continue;
        a->residue = pose;
        a->stream_pdb_line(lpdbdat, 9000+l, true);
    }

    #if compute_lsrb
    for (l=0; l<nlsrb_points; l++)
    {
        if (lsrb_points[l].weight > 0.01)
            lpdbdat << "REMARK 503 " << lsrb_points[l].x << " " << lsrb_points[l].y << " " << lsrb_points[l].z << " " << lsrb_points[l].weight << endl;
    }
    #endif

    if (waters)
    {
        for (k=0; waters[k]; k++)
        {
            for (l=0; l<3; l++)
            {
                Atom* a = waters[k]->get_atom(l);
                if (!a) continue;
                a->residue = pose;
                a->stream_pdb_line(lpdbdat, 9000+offset+l+3*k, true);
            }
        }
    }

    int en = protein->get_end_resno();
    int resno;
    for (resno = protein->get_start_resno(); resno <= en; resno++)
    {
        AminoAcid* laa = protein->get_residue(resno);
        if (!laa) continue;
        if (!laa->been_flexed && !laa->coordmtl)
        {
            if (laa->distance_to(ligand) > 5) continue;
            for (k=0; reaches_spheroid[k]; k++)
            {
                if (!protein->aa_ptr_in_range(reaches_spheroid[k])) continue;
                if (reaches_spheroid[k] == laa) goto _afterall;
            }
            continue;
        }
        _afterall:
        n = laa->get_atom_count();
        for (l=0; l<n; l++)
        {
            laa->get_atom(l)->stream_pdb_line(
                lpdbdat,
                laa->atno_offset+l
            );
        }

        if (laa->coordmtl && laa->coordmtl->residue == laa->get_residue_no())
        {
            laa->coordmtl->stream_pdb_line(lpdbdat, laa->coordmtl->residue);
        }
    }

    this->pdbdat = lpdbdat.str();
}

bool DockResult::clashes_with(DockResult *o)
{
    Molecule lig1(*mlig), ligA(*mlig);
    ligpos.restore_state(&lig1);
    o->ligpos.restore_state(&ligA);
    float clash = lig1.get_intermol_clashes(&ligA);
    if (clash > clash_limit_per_aa) return true;
    return false;
}

DockResult DockResult::merge(DockResult *other)
{
    // Have to identify residues by which copy of the ligand they're more interactive with,
    // and then position the residues according to the proper starting object.
    //
    // Drat.
    //
    // This is going to be much more involved than it seemed at first.
    //
    return DockResult();
}

std::ostream& operator<<(std::ostream& output, const DockResult& dr)
{
    int l;

    if (dr.isomer.length())
    {
        output << "Isomer: " << dr.isomer << endl << endl;
    }

    if (dr.mbbr)
    {
        output << "Best-binding targets:" << endl;
        output << *(dr.mbbr) << endl << endl;
    }

    if (dr.out_per_res_e || dr.out_per_btyp_e)
    {
        output << "# Binding energies:" << endl;
    }

    if (dr.out_per_res_e)
    {
        output << "BENERG:" << endl;
        for (	l=0;

                dr.mkJmol
                && dr.metric
                && dr.metric[l]
                && dr.metric[l][0]
                ;

                l++
            )
        {
            if (!dr.is_mcr[l] && fabs(dr.mkJmol[l]) < dr.out_itemized_e_cutoff) continue;

            if (dr.do_output_colors) colorize(dr.mkJmol[l]);
            output << dr.metric[l] << ": " << dr.mkJmol[l]*dr.energy_mult << "{" << dr.mstab[l] << "}";
            if (dr.display_binding_atoms) output << " |";
            if (dr.display_binding_atoms) output << " " << (dr.mb_atom1_name[l] ? dr.mb_atom1_name[l] : "-");
            if (dr.display_binding_atoms) output << " " << (dr.mb_atom2_name[l] ? dr.mb_atom2_name[l] : "-");
            if (dr.display_clash_atoms) output << " |";
            if (dr.display_clash_atoms) output << " " << (dr.mc_atom1_name[l] ? dr.mc_atom1_name[l] : "-");
            if (dr.display_clash_atoms && dr.mc_atom1_name[l] && dr.mc_atom2_name[l]) output << " clash";
            if (dr.display_clash_atoms) output << " " << (dr.mc_atom2_name[l] ? dr.mc_atom2_name[l] : "-");
            output << endl;
            if (dr.do_output_colors) colorless();
        }
        output << endl;
        // output << "Worst energy: " << dr.worst_energy << endl;
    }

    if (dr.out_per_btyp_e) for (l=0; l<_INTER_TYPES_LIMIT; l++)
    {
        if (fabs(dr.bytype[l]) < dr.out_itemized_e_cutoff) continue;

        char lbtyp[64];
        switch (l+covalent)
        {
            case covalent:
                continue; /*strcpy(lbtyp, "Total covalent: ");		break;*/
            case ionic:
                strcpy(lbtyp, "Total ionic: ");
                break;
            case hbond:
                strcpy(lbtyp, "Total H-bond: ");
                break;
            case pi:
                strcpy(lbtyp, "Total pi stack: ");
                break;
            case polarpi:
                strcpy(lbtyp, "Total polar-pi and cation-pi: ");
                break;
            case mcoord:
                strcpy(lbtyp, "Total metal coordination: ");
                break;
            case vdW:
                strcpy(lbtyp, "Total van der Waals: ");
                break;
            default:
                goto _btyp_unassigned;
        }

        if (dr.do_output_colors) colorize(dr.bytype[l]);
        output << lbtyp << dr.bytype[l]*dr.energy_mult << endl;
        if (dr.do_output_colors) colorless();
    }
    output << endl;

_btyp_unassigned:

    output << "Atom-level interactions:" << endl;
    output << dr.atomlvl;
    output << endl;

    // output << "Ligand volume: " << dr.ligvol << endl;
    if (dr.ligand_waters_energy) output << "Pocket ligand-water energy: " << dr.ligand_waters_energy << endl;
    if (dr.do_output_colors) colorize(dr.kJmol);
    output << "Total: " << dr.kJmol*dr.energy_mult << endl;
    if (dr.do_output_colors) colorless();
    if (dr.do_output_colors) colorize(dr.kJmol_raw);
    output << "Total ignoring stability: " << dr.kJmol_raw*dr.energy_mult << endl;
    if (dr.do_output_colors) colorless();
    output << "Ligand pocket occlusion: " << dr.ligand_pocket_occlusion << endl;
    output << "Ligand solvation energy: " << dr.ligand_solvation_energy*dr.energy_mult << endl;
    output << "Ligand pocket solvation energy: " << dr.ligand_pocket_wet_energy*dr.energy_mult << endl;
    output << "Pocket hydration energy: " << dr.pocket_wet_solvation_energy*dr.energy_mult << endl;
    output << "Pocket bound hydration energy: " << dr.pocket_bound_solvation_energy*dr.energy_mult << endl;
    output << "Solvation energy stabilizing internal contacts: " << dr.pocket_ic_DeltaG_solvation*dr.energy_mult << endl;
    #if compute_lsrb
    output << "Ligand surface receptor binding: " << dr.ligand_surface_receptor_binding << endl;
    #endif
    if (dr.mcoord_charge_repulsion) output << "Metal-coordinated charge repulsion: " << dr.mcoord_charge_repulsion << endl;
    if (dr.stay_close_ligand)
    {
        output << "Maintained " << dr.stay_close_ligand->name;
        Atom* heavy = nullptr;
        if (dr.stay_close_ligand->Z == 1 && dr.stay_close_ligand->get_bonded_heavy_atoms_count()) heavy = dr.stay_close_ligand->get_bond_by_idx(0)->atom2;
        if (heavy) output << " (" << heavy->name << ")"; 
        output << " ~ "
            << dr.stay_close_protein->aa3let << dr.stay_close_protein->residue << ":" << dr.stay_close_protein->name
            << " at " << dr.stay_close_protein->distance_to(dr.stay_close_ligand) << " A." << endl;
    }
    output << "Worst atom clash: " << dr.worst_energy*dr.energy_mult << endl;

    if (dr.out_lig_int_e) output << "Ligand internal energy: " << dr.ligand_self*dr.energy_mult << endl << endl;

    if (dr.estimated_TDeltaS)
    {
        output << "Estimated TDeltaS: " << dr.estimated_TDeltaS << endl;
    }

    if (dr.miscdata.size())
    {
        output << dr.miscdata << endl;
    }
    output << "A100 score: " << dr.A100 << endl;
    output << endl;

    if (dr.out_lig_pol_sat)
    {
        output << "Ligand polar satisfaction: " << dr.polsat << endl;
        output << endl;
    }

    if (dr.out_prox) output << "Proximity: " << dr.proximity << endl << endl;

    if (dr.out_pro_clash) output << "Protein clashes: " << dr.protclash << endl << endl;

    #if compute_missed_connections
    if (dr.out_mc)
    {
        output << "# Missed Connections" << endl << "MC:" << endl;
        output << endl;
        
        for (	l=0;

                dr.metric
                && dr.metric[l]
                && dr.metric[l][0];

                l++
            )
        {
            if (fabs(dr.missed_connections[l]) < dr.out_itemized_e_cutoff) continue;
            if (dr.metric[l]) output << dr.metric[l] << ": " << dr.missed_connections[l]*dr.energy_mult << endl;
        }
        output << endl;
    }
    #endif

    #if compute_vdw_repulsion
    if (dr.out_vdw_repuls)
    {
        output << "# van der Waals repulsion" << endl << "vdWRPL:" << endl;
        output << endl;
        
        for (	l=0;

                dr.metric
                && dr.metric[l]
                && dr.metric[l][0];

                l++
            )
        {
            if (fabs(dr.mvdWrepl[l]) < dr.out_itemized_e_cutoff) continue;
            if (dr.metric[l]) output << dr.metric[l] << ": " << dr.mvdWrepl[l]*dr.energy_mult << endl;
        }
        output << endl;
    }
    #endif

    if (dr.disqualified)
    {
        if (dr.disqualify_reason.length()) output << "Disqualified because: " << dr.disqualify_reason << endl;
        else output << "Disqualified for reasons." << endl;
    }

    if (dr.include_pdb_data)
    {
        if (!dr.pdbdat.length())
        {
            output << "WARNING: Missing PDB data." << endl;
        }
        else
        {
            output << "# PDB Data" << endl << "PDBDAT:" << endl;

            output << dr.pdbdat << endl;

            output << "TER" << endl << "END" << endl << endl << endl;
        }
    }

    return output;
}
