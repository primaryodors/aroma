
#include "search.h"
#include "moiety.h"

Point loneliest;
Point size(10,10,10);
std::vector<int> exclusion;

int agqty = 0;
AminoAcid* cs_res[MAX_CS_RES];
intera_type cs_bt[MAX_CS_RES];
int cs_res_qty = 0;
int cs_idx = -1;
bool Search::any_resnos_priority = false;

void Search::do_tumble_spheres(Protein* protein, Molecule* ligand, Point l_pocket_cen)
{
    int i, j, l, n;
    float lig_min_int_clsh = ligand->get_internal_clashes();

    // Begin tumble sphere behavior.
    AminoAcid* tsphres[SPHREACH_MAX+4];
    int tsphsz = protein->fetch_residues_near(l_pocket_cen, size.magnitude()+6, tsphres);
    float outer_sphere[tsphsz+4], inner_sphere[tsphsz+4];

    for (i=0; i<tsphsz+4; i++) outer_sphere[i] = inner_sphere[i] = 0;

    Point pocketsize = protein->estimate_pocket_size(tsphres);
    Point ligbbox = ligand->get_bounding_box();

    for (i=0; !ligbbox.fits_inside(pocketsize) && i<100; i++)
    {
        ligand->crumple(fiftyseventh*30);
    }

    for (i=0; i<tsphsz; i++)
    {
        #if use_exclusions
        if (exclusion.size()
                &&
                std::find(exclusion.begin(), exclusion.end(), tsphres[i]->get_residue_no())!=exclusion.end()
        )
        {
            for (j=i; j<tsphsz; j++) tsphres[j] = tsphres[j+1];
            i--;
            tsphsz--;
            continue;
        }
        #endif

        // TODO: Algorithmically determine more accurate values based on interaction type, etc.
        outer_sphere[i] = tsphres[i]->get_reach() + 2.5;
        inner_sphere[i] = tsphres[i]->get_reach() / 3 + 1;
    }

    const SCoord xaxis = Point(1,0,0), yaxis = Point(0,1,0), zaxis = Point(0,0,1);
    float loneliness=0, blone=0, xrad, yrad, zrad, lrad, step, bestxr, bestyr, bestzr, score, worth, weight, bestscore;
    const int ac = ligand->get_atom_count();
    Pose besp(ligand);
    #if _DBG_TUMBLE_SPHERES
    std::string tsdbg = "", tsdbgb = "";
    #endif

    if (ligbbox.x > ligbbox.y && pocketsize.x < pocketsize.y) ligand->rotate(zaxis, square);
    if (ligbbox.x > ligbbox.z && pocketsize.x < pocketsize.z) ligand->rotate(yaxis, square);
    if (ligbbox.y > ligbbox.x && pocketsize.y < pocketsize.x) ligand->rotate(zaxis, square);
    if (ligbbox.y > ligbbox.z && pocketsize.y < pocketsize.z) ligand->rotate(xaxis, square);
    if (ligbbox.z > ligbbox.x && pocketsize.z < pocketsize.x) ligand->rotate(yaxis, square);
    if (ligbbox.z > ligbbox.y && pocketsize.z < pocketsize.y) ligand->rotate(xaxis, square);

    step = fiftyseventh*30;
    bestscore = -Avogadro;
    float lonely_step = 1.0 / loneliest.get_3d_distance(l_pocket_cen);
    #if _DBG_LONELINESS
    cout << "Loneliest point " << loneliest << " is " << loneliest.get_3d_distance(l_pocket_cen) << "A from pocketcen " << l_pocket_cen << "." << endl;
    cout << "Pocket size is " << pocketsize << " vs. ligand bounding box " << ligbbox << endl;
    #endif
    if (isnan(lonely_step) || lonely_step < 0.1) lonely_step = 0.1;

    #if pocketcen_is_loneliest
    if (1)
    {
        ligand->recenter(l_pocket_cen);
    #else
    for (loneliness=0; loneliness <= 1; loneliness += lonely_step)
    {
        float centeredness = 1.0 - loneliness;
        Point tmpcen(loneliest.x * loneliness + l_pocket_cen.x * centeredness,
                     loneliest.y * loneliness + l_pocket_cen.y * centeredness,
                     loneliest.z * loneliness + l_pocket_cen.z * centeredness
                    );
        ligand->recenter(tmpcen);
    #endif

        #if _DBG_LONELINESS && !pocketcen_is_loneliest
        cout << "Ligand is " << loneliness << " lonely centered at " << tmpcen << "." << endl;
        #endif

        for (xrad=0; xrad <= M_PI*2; xrad += step)
        {
            for (yrad=0; yrad <= M_PI*2; yrad += step)
            {
                for (zrad=0; zrad <= M_PI*2; zrad += step)
                {
                    ligbbox = ligand->get_bounding_box();

                    if (ligbbox.x > ligbbox.y && pocketsize.x < pocketsize.y) continue;
                    if (ligbbox.x > ligbbox.z && pocketsize.x < pocketsize.z) continue;
                    if (ligbbox.y > ligbbox.x && pocketsize.y < pocketsize.x) continue;
                    if (ligbbox.y > ligbbox.z && pocketsize.y < pocketsize.z) continue;
                    if (ligbbox.z > ligbbox.x && pocketsize.z < pocketsize.x) continue;
                    if (ligbbox.z > ligbbox.y && pocketsize.z < pocketsize.y) continue;

                    Bond** rb = ligand->get_rotatable_bonds();

                    if (!rb) n = 0;
                    else for (n=0; rb[n]; n++);		// Get count.

                    l = 0;
                    lrad = 0;
                _xyzl_loop:
                    if (ligand->get_internal_clashes() >= lig_min_int_clsh*5+5) goto _xyzl_skip_loop;

                    score = 0;
                    #if _DBG_TUMBLE_SPHERES
                    tsdbg = "";
                    // cout << ligand->get_internal_clashes() << " vs. " << lig_min_int_clsh << endl;
                    #endif
                    for (i=0; i<ac; i++)
                    {
                        Atom* a = ligand->get_atom(i);
                        intera_type it = vdW;

                        for (j=0; j<tsphsz; j++)
                        {
                            worth = 0.4;
                            if (fabs(a->get_charge()) >= hydrophilicity_cutoff && fabs(tsphres[j]->get_charge()) >= hydrophilicity_cutoff
                                    &&
                                    sgn(a->get_charge()) == -sgn(tsphres[j]->get_charge())
                            )
                            {
                                it = ionic;
                                worth = 100;
                            }
                            else if (fabs(a->get_charge()) >= hydrophilicity_cutoff || fabs(a->is_polar()) >= hydrophilicity_cutoff)
                            {
                                it = hbond;
                                worth = 40;
                            }
                            else if (a->is_pi())
                            {
                                it = pi;
                                worth = 7;
                            }

                            if (tsphres[j]->capable_of_inter(it))
                            {
                                float r = a->loc.get_3d_distance(tsphres[j]->get_atom_location("CA"));
                                if (r <= outer_sphere[j])
                                {
                                    if (r > inner_sphere[j])
                                    {
                                        weight = 1;

                                        if (tsphres[j]->priority)
                                        {
                                            weight = ts_priority_coefficient;		// Extra weight for residues mentioned in a CEN RES or PATH RES parameter.
                                        }

                                        if (tsphres[j]->ring_is_aromatic(0) && fabs(a->is_polar()) >= hydrophilicity_cutoff) weight /= 2;

                                        #if !tumble_spheres_include_vdW
                                        if ((worth*weight) < 1) continue;
                                        #endif

                                        score += worth * weight;
                                        #if _DBG_TUMBLE_SPHERES
                                        tsdbg += std::string("+ ")
                                                +  std::string(a->name) + std::string(" reaches ") + std::string(tsphres[j]->get_3letter())
                                                +  std::to_string(tsphres[j]->get_residue_no()) + std::string(".  ");
                                        #endif
                                    }
                                    else
                                    {
                                        score -= 200;
                                        #if _DBG_TUMBLE_SPHERES
                                        tsdbg += std::string("- ")
                                                +  std::string(a->name) + std::string(" clashes ") + std::string(tsphres[j]->get_3letter())
                                                +  std::to_string(tsphres[j]->get_residue_no()) + std::string(".  ");
                                        #endif
                                    }
                                }
                            }
                        }
                    }

                    #if !pocketcen_is_loneliest
                    if (score > 0) score *= 1.0 + 0.1 * centeredness;
                    #endif

                    if (score > bestscore)
                    {
                        besp.copy_state(ligand);
                        blone = loneliness;
                        bestxr = xrad;
                        bestyr = yrad;
                        bestzr = zrad;
                        bestscore = score;

                        #if _DBG_TUMBLE_SPHERES
                        tsdbgb = tsdbg;

                        cout << "Tumble score " << score << " for ligand box " << ligand->get_bounding_box() << endl;


                        #if output_tumble_debug_docs
                        int u, v, w;
                        char protfttl[1000];
                        strcpy(protfttl, protfname);

                        char** lwords = chop_spaced_words(protfttl, '/');

                        for (u=0; lwords[u]; u++);
                        u--;

                        char fname[1000];
                        sprintf(fname, "output/tumble_%s_%d_%d_%d_%f.dock",
                                lwords[u],
                                (int)(xrad*fiftyseven),
                                (int)(yrad*fiftyseven),
                                (int)(zrad*fiftyseven),
                                score);
                        cout << fname << endl;
                        std::ofstream tspdbdat(fname, std::ofstream::out);

                        tspdbdat << "PDB file: " << protfname << endl;
                        tspdbdat << "Pose: 1\nNode: 0\nPDBDAT:\n";

                        int lac = ligand->get_atom_count();
                        for (u=0; u<lac; u++) ligand->get_atom(u)->stream_pdb_line(tspdbdat, 9000+u);

                        int pseql = protein->get_seq_length();
                        v = 1;
                        for (u = 1; u < pseql; u++)
                        {
                            AminoAcid* dbgaa = protein->get_residue(u);
                            if (dbgaa)
                            {
                                int aaac = dbgaa->get_atom_count();
                                for (w=0; w<aaac; w++)
                                {
                                    Atom* dbga = dbgaa->get_atom(w);
                                    if (!strcmp(dbga->name, "CA") || !strcmp(dbga->name, "CB")) dbga->stream_pdb_line(tspdbdat, v++);
                                }
                            }
                        }
                        tspdbdat << "END" << endl;
                        tspdbdat.close();
                        #endif
                        #endif
                    }

                _xyzl_skip_loop:

                    if (rb && rb[l])
                    {
                        rb[l]->rotate(step);

                        lrad += step;
                        if (lrad >= M_PI*2)
                        {
                            l++;
                            if (l < n) goto _xyzl_loop;
                        }
                        else goto _xyzl_loop;
                    }

                    ligand->rotate(zaxis, step);
                }		// zrad.
                ligand->rotate(yaxis, step);
            }			// yrad.
            ligand->rotate(xaxis, step);
        }				// xrad.

        #if !pocketcen_is_loneliest
        if (bestscore >= (ligand->get_atom_count()*13)) break;
        #endif

        #if _DBG_LONELINESS
        cout << "Best score: " << bestscore << endl;
        #endif
    }					// loneliness.

    #if _DBG_TUMBLE_SPHERES
    cout << "Tumble sphere best score " << bestscore << " for "
        << "x" << bestxr*fiftyseven << "deg, "
        << "y" << bestyr*fiftyseven << "deg, "
        << "z" << bestzr*fiftyseven << "deg."
        << " (" << blone << " lonely)."
        << endl;
    cout << tsdbgb << endl;
    #endif

    // Load the best ligand conformer.
    besp.restore_state(ligand);

    // Minimize ligand clashes.
    #if prerot_sidechains_from_ligand
    for (i=0; i<tsphsz; i++)
    {
        Bond** tsphb = tsphres[i]->get_rotatable_bonds();
        if (tsphb)
        {
            for (j=0; tsphb[j]; j++)
            {
                float rad=0, bestrad=0, clash, bestclash=6.25e24;
                for (rad=0; rad < M_PI*2; rad += step)
                {
                    clash = tsphres[i]->get_intermol_clashes(ligand);

                    if (clash < bestclash)
                    {
                        bestrad = rad;
                        bestclash = clash;
                    }

                    tsphb[j]->rotate(step);
                }

                tsphb[j]->rotate(bestrad);
            }
            // delete[] tsphb;
        }
    }
    #endif
    // End tumble sphere behavior.
}

void Search::copy_ligand_position_from_file(Protein* protein, Molecule* ligand, const char* filename, const char* ligname, int resno)
{
    char buffer[4096];
    FILE* fp = fopen(filename, "rb");
    if (!fp)
    {
        cout << "Failed to open " << filename << " for reading." << endl;
        throw 0xbadf12e;
    }

    int lused = rand();
    bool copying = false;
    while (!feof(fp))
    {
        char* got = fgets(buffer, 4090, fp);
        if (!got) break;
        char** words = chop_spaced_words(buffer);
        if (!words[0] || strcmp(words[0], "HETATM"))
        {
            if (copying && (resno >= 0 || frand(0,1) < 0.4)) break;
            else
            {
                copying = false;
                continue;
            }
        }
        if (ligname && strlen(ligname) && strcmp(ligname, words[3]))
        {
            cout << ligname << " != " << words[3] << endl;
            if (copying && (resno >= 0 || frand(0,1) < 0.4)) break;
            else
            {
                copying = false;
                continue;
            }
        }
        if (resno >= 0 && resno != atoi(words[4]))
        {
            cout << resno << " != " << words[4] << endl;
            if (copying && (resno >= 0 || frand(0,1) < 0.4)) break;
            else
            {
                copying = false;
                continue;
            }
        }

        Atom* a = ligand->get_atom(words[2]);
        if (!a) continue;
        Point pt(atof(words[5]), atof(words[6]), atof(words[7]));
        a->move(pt);
        a->used = lused;
        copying = true;
    }

    fclose(fp);

    int i, n = ligand->get_atom_count();
    bool any_unmoved = false;
    for (i=0; i<n; i++)
    {
        Atom* a = ligand->get_atom(i);
        if (!a) continue;
        if (a->used != lused)
        {
            if (a->Z > 1)
            {
                cout << "ERROR: Source file does not contain all ligand atoms (missing " << a->name << ")." << endl;
                throw 0xbadda7a;
            }
            any_unmoved = true;
            break;
        }
    }

    if (any_unmoved)
    {
        ligand->dehydrogenate();
        ligand->hydrogenate();
    }
}

int Search::identify_ligand_pairing_targets(Molecule *ligand, LigandTarget *results, int max_results)
{
    int found = 0, i, j, n = ligand->get_heavy_atom_count();
    for (i=0; i<n; i++)
    {
        Atom* a = ligand->get_atom(i);
        if (a->Z < 2) continue;
        if (!a->conjugation
            || (a->conjugation->num_heavy_atoms > 6)
            || (fabs(a->is_polar()) > hydrophilicity_cutoff)
            || (a->get_family() == PNICTOGEN)
            || (a->get_family() == CHALCOGEN)
           )
        {
            if (a->get_family() == TETREL)
            {
                if (a->conjugation) continue;
                int nhb = a->get_bonded_heavy_atoms_count();
                if (nhb > 1) continue;
            }

            results[found].single_atom = a;
            results[found].conjgrp = nullptr;
            found++;
        }
        if (a->conjugation)
        {
            bool already = false;
            if (a->conjugation->num_heavy_atoms < 0.666*n)
            {
                for (j=0; j<found; j++) if (results[j].conjgrp == a->conjugation) already = true;
                if (!already)
                {
                    results[found].conjgrp = a->conjugation;
                    results[found].single_atom = nullptr;
                    found++;
                }
            }
            if (a->num_rings())
            {
                Ring** r = a->get_rings();
                if (r[0]->get_atom_count() <= 6)
                {
                    Conjugation* rc = new Conjugation(r[0]);
                    already = false;
                    for (j=0; j<found; j++) if (results[j].conjgrp && rc->contains(results[j].conjgrp)) already = true;
                    if (!already)
                    {
                        results[found].conjgrp = rc;
                        results[found].single_atom = nullptr;
                        found++;
                    }
                    else delete rc;
                }
            }
        }
    }
    return found;
}

void Search::pair_targets(Molecule *ligand, LigandTarget *targets, AminoAcid **pocketres, Point loneliest, 
    BestBindingResult* output, Cavity* container, bool allow_thiolation)
{
    int i, j, k, l, m, n, ii;
    float best = 0;
    bool best_ionic = false, best_neutral_polar = false, best_pi = false;
    Point ligcen = ligand->get_barycenter();
    int ntarg, npr;
    bool override_target_compatibility = false;
    BestBindingResult bbr;
    float score;
    int multichalc = 0;
    Atom* mcatoms[256];

    // TODO:
    // Right now you are considering multiple possible starting poses and then selecting hopably the best (most energetically favorable) one.
    // What you can do is rule out the ones that cannot avoid clashes with the protein, and then store what's left over in order to estimate
    // the probabilities of all the possible states. Then the result will allow estimating the TΔS of the optimal state.
    // https://en.wikipedia.org/wiki/Boltzmann_distribution

    memset(mcatoms, 0, 256*sizeof(Atom*));
    Moiety y;
    if (!ligand->num_monomers && !ligand->get_charge())
    {
        y.pattern = "OCCO";
        if (!y.contained_by(ligand, mcatoms)) y.pattern = "OCCCO";
        multichalc = y.contained_by(ligand, mcatoms);

        if (multichalc)
        {
            Atom *O1 = nullptr, *O2 = nullptr;
            for (ii=0; mcatoms[ii]; ii++)
            {
                // cout << mcatoms[ii]->name << " ";
                if (mcatoms[ii]->Z == 8)
                {
                    if (!O1) O1 = mcatoms[ii];
                    else if (!O2) O2 = mcatoms[ii];
                    if (O2) break;
                }

            }

            if (O1 && O2)
            {
                ligand->conform_atom_to_location(O2, O1, 50, 2);
                // cout << O1->distance_to(O2);
                for (ii=1; mcatoms[ii]; ii++)
                {
                    if (ii)
                    {
                        if (mcatoms[ii-1]->Z != 8 && mcatoms[ii]->Z != 8)
                        {
                            Bond* b = mcatoms[ii-1]->get_bond_between(mcatoms[ii]);
                            if (b) b->can_rotate = b->can_flip = false;
                            b = mcatoms[ii]->get_bond_between(mcatoms[ii-1]);
                            if (b) b->can_rotate = b->can_flip = false;
                        }
                    }
                }
            }
            // cout << endl;
        }
    }

    for (ntarg=0; targets[ntarg].conjgrp || targets[ntarg].single_atom; ntarg++);       // count
    for (npr=0; pocketres[npr]; npr++);                                                 // count

    #if _dbg_bb_scoring
    cout << ntarg << " targets, " << i << " pocket residues." << endl;
    #endif

    _new_attempt_no_tgt_comp:
    for (i=0; i<ntarg; i++)
    {
        float ichg = targets[i].charge();
        float ipol = targets[i].polarity();
        if (ipol < hydrophilicity_cutoff) ipol = 0;
        int ifam = targets[i].single_atom ? targets[i].single_atom->get_family() : -1;
        int iZ = targets[i].single_atom ? targets[i].single_atom->Z : -1;
        if (ifam == PNICTOGEN) ipol += hydrophilicity_cutoff*2;
        int ipi = targets[i].conjgrp ? targets[i].conjgrp->num_heavy_atoms : 0;
        bool ihba = targets[i].has_hb_acceptors();
        bool ihbd = targets[i].has_hb_donors();
        Point icen = targets[i].barycenter();
        bool imc = false;
        if (multichalc) for (ii=0; mcatoms[ii]; ii++) if (targets[i].contains(mcatoms[ii])) imc = true;

        for (j=0; pocketres[j]; j++)
        {
            float jchg = pocketres[j]->get_charge();
            float jpol = pocketres[j]->hydrophilicity();
            if (jpol < hydrophilicity_cutoff) jpol = 0;
            int jpi = pocketres[j]->has_pi_atoms(false);
            float jcba = pocketres[j]->CB_angle(loneliest);
            bool jhba = pocketres[j]->has_hbond_acceptors();
            bool jhbd = pocketres[j]->has_hbond_donors();
            Atom* jmtl = pocketres[j]->coordmtl;

            bool ijmc = jmtl && (ifam == CHALCOGEN || ifam == PNICTOGEN) && ((iZ != 8 && ichg <= 0) || imc);

            #if 0
            int jrno = pocketres[j]->get_residue_no();
            if (jrno == 111) 
                cout << jrno 
                    << " ichg " << ichg << " jchg " << jchg
                    << " ipol " << ipol << " jpol " << jpol
                    << " ipi " << ipi << " jpi " << jpi
                    << " ihba " << ihba << " jhba " << jhba
                    << " ihbd " << ihbd << " jhbd " << jhbd
                    << (target_compatibility(ichg, jchg, ipol, jpol, ipi, jpi, ihba, jhba, ihbd, jhbd) ? " comp" : " inc")
                    << endl;
            #endif

            if (!override_target_compatibility) if (!ijmc && !target_compatibility(ichg, jchg, ipol, jpol, ipi, jpi, ihba, jhba, ihbd, jhbd)) continue;

            if (ntarg < 2 || npr < 2)
            {

                bbr.ligand = ligand;
                bbr.pri_tgt = &targets[i];
                bbr.pri_res = pocketres[j];
                bbr.sec_tgt = nullptr;
                bbr.sec_res = nullptr;
                bbr.tert_tgt = nullptr;
                bbr.tert_res = nullptr;

                if (container) align_targets(ligand, bbr.barycenter(), &bbr);
                score = bbr.score(ligcen, container) * (1.0 + frand(0,bb_stochastic));

                if (score > best)
                {
                    *output = bbr;
                    output->probability = 1;
                    best = score;
                }

                continue;
            }

            for (k=0; k<ntarg; k++)
            {
                if (ntarg>=2 && k<=i) continue;
                if (targets[k].contains(&targets[i])) continue;
                if (targets[i].contains(&targets[k])) continue;
                float kchg = targets[k].charge();
                float kpol = targets[k].polarity();
                if (kpol < hydrophilicity_cutoff) kpol = 0;
                int kfam = targets[k].single_atom ? targets[k].single_atom->get_family() : -1;
                int kZ = targets[k].single_atom ? targets[k].single_atom->Z : -1;
                if (kfam == PNICTOGEN) kpol += hydrophilicity_cutoff*2;
                int kpi = targets[k].conjgrp ? targets[k].conjgrp->num_heavy_atoms : 0;
                Point kcen = targets[k].barycenter();
                bool khba = targets[k].has_hb_acceptors();
                bool khbd = targets[k].has_hb_donors();
                bool kmc = false;
                if (multichalc) for (ii=0; mcatoms[ii]; ii++) if (targets[k].contains(mcatoms[ii])) kmc = true;
                for (l=0; pocketres[l]; l++)
                {
                    if (l==j) continue;
                    float lchg, lpol, lcba;
                    int lpi;
                    bool lhba, lhbd, klmc;
                    Atom* lmtl;

                    lchg = pocketres[l]->get_charge();
                    lpol = pocketres[l]->hydrophilicity();
                    if (lpol < hydrophilicity_cutoff) lpol = 0;
                    lpi = pocketres[l]->has_pi_atoms(false);
                    lcba = pocketres[l]->CB_angle(loneliest);
                    lhba = pocketres[l]->has_hbond_acceptors();
                    lhbd = pocketres[l]->has_hbond_donors();
                    lmtl = pocketres[l]->coordmtl;

                    klmc = lmtl && (kfam == CHALCOGEN || kfam == PNICTOGEN) && ((kZ != 8 && kchg <= 0) || kmc);
                    if (!override_target_compatibility) if (!klmc && !target_compatibility(kchg, lchg, kpol, lpol, kpi, lpi, khba, lhba, khbd, lhbd)) continue;

                    for (m=-1; m<ntarg; m++)
                    {
                        float mchg, mpol;
                        int mfam, mZ, mpi;
                        Point mcen;
                        bool mhba, mhbd;

                        if (m >= 0)
                        {
                            if (/*m==i ||*/ m<=k) continue;
                            if (targets[m].contains(&targets[i])) continue;
                            if (targets[i].contains(&targets[m])) continue;
                            if (targets[m].contains(&targets[k])) continue;
                            if (targets[k].contains(&targets[m])) continue;

                            mchg = targets[m].charge();
                            mpol = targets[m].polarity();
                            if (mpol < hydrophilicity_cutoff) mpol = 0;
                            mfam = targets[m].single_atom ? targets[m].single_atom->get_family() : -1;
                            mZ = targets[m].single_atom ? targets[m].single_atom->Z : -1;
                            if (mfam == PNICTOGEN) mpol += hydrophilicity_cutoff*2;
                            mpi = targets[m].conjgrp ? targets[m].conjgrp->num_heavy_atoms : 0;
                            mcen = targets[m].barycenter();
                            mhba = targets[m].has_hb_acceptors();
                            mhbd = targets[m].has_hb_donors();
                        }
                        for (n=0; pocketres[n]; n++)
                        {
                            if (m >= 0 && ntarg>=3 && (n==j || n==l)) continue;

                            bool has_ionic = false, has_neutral_polar = false, has_pi = false;
                            if (ichg && jchg) has_ionic = true;
                            else if (ipol && jpol) has_neutral_polar = true;
                            else if (ipi && jpi) has_pi = true;
                            if (kchg && lchg) has_ionic = true;
                            else if (lpol && kpol) has_neutral_polar = true;
                            else if (kpi && lpi) has_pi = true;

                            float nchg, npol, ncba;
                            int npi;
                            bool nhba, nhbd, mnmc;
                            Atom* nmtl;
                            if (m >= 0)
                            {
                                nchg = pocketres[n]->get_charge();
                                npol = pocketres[n]->hydrophilicity();
                                npi = pocketres[n]->has_pi_atoms(false);
                                ncba = pocketres[n]->CB_angle(loneliest);
                                nhba = pocketres[n]->has_hbond_acceptors();
                                nhbd = pocketres[n]->has_hbond_donors();
                                nmtl = pocketres[n]->coordmtl;

                                mnmc = pocketres[n]->coordmtl && (mfam == CHALCOGEN || mfam == PNICTOGEN) && mZ != 8 && mchg <= 0;
                                if (!override_target_compatibility) 
                                    if (!mnmc && !target_compatibility(mchg, nchg, mpol, npol, mpi, npi, mhba, nhba, mhbd, nhbd))
                                        continue;

                                if (mchg && nchg) has_ionic = true;
                                else if (mpol && npol) has_neutral_polar = true;
                                else if (mpi && npi) has_pi = true;
                            }

                            bbr.ligand = ligand;
                            bbr.pri_tgt = &targets[i];
                            bbr.pri_res = pocketres[j];
                            bbr.sec_tgt = &targets[k];
                            bbr.sec_res = pocketres[l];
                            bbr.tert_tgt = (m<0) ? nullptr : &targets[m];
                            bbr.tert_res = (m<0) ? nullptr : pocketres[n];

                            if (container) align_targets(ligand, bbr.barycenter(), &bbr);
                            score = bbr.score(ligcen, container) * (1.0 + frand(0,bb_stochastic));

                            #if bb_avoid_eclipsing_contacts
                            // Prefer contacts without intervening residues.
                            int h;
                            float jldist = bbr.pri_res->distance_to(bbr.sec_res), jndist, lndist;
                            Atom *aj, *al, *an;
                            aj = bbr.pri_res->get_nearest_atom(loneliest);
                            al = bbr.sec_res->get_nearest_atom(loneliest);
                            if (m >= 0)
                            {
                                jndist = bbr.pri_res->distance_to(bbr.tert_res);
                                lndist = bbr.sec_res->distance_to(bbr.tert_res);
                                an = bbr.tert_res->get_nearest_atom(loneliest);
                            }
                            if (aj && al)
                            {
                                for (h=0; pocketres[h]; h++)
                                {
                                    if (h == j || h == l) continue;
                                    if (m >= 0 && h == n) continue;

                                    Atom* ha = pocketres[h]->get_nearest_atom(loneliest);
                                    if (!ha) continue;

                                    float hjdist, hldist, hndist;
                                    hjdist = pocketres[h]->distance_to(bbr.pri_res);
                                    hldist = pocketres[h]->distance_to(bbr.sec_res);
                                    if (m >= 0)
                                    {
                                        hndist = pocketres[h]->distance_to(bbr.tert_res);
                                    }

                                    if (hjdist < jldist && hldist < jldist)
                                    {
                                        ha = pocketres[h]->get_nearest_atom_to_line(aj->loc, al->loc);
                                        if (ha)
                                        {
                                            float hr = ha->loc.get_distance_to_line(aj->loc, al->loc);
                                            if (hr < bb_eclipsing_limit)
                                            {
                                                score /= bb_eclipsing_divisor;
                                                #if _dbg_eclipsing_contacts
                                                cout << bbr.pri_res->get_name() << " --- " << bbr.sec_res->get_name()
                                                    << " eclipsed by " << pocketres[h]->get_name() << " r=" << hr << endl;
                                                #endif
                                            }
                                        }
                                    }

                                    if (m >= 0)
                                    {
                                        if (hjdist < jndist && hndist < jndist)
                                        {
                                            ha = pocketres[h]->get_nearest_atom_to_line(aj->loc, an->loc);
                                            if (ha)
                                            {
                                                float hr = ha->loc.get_distance_to_line(aj->loc, an->loc);
                                                if (hr < bb_eclipsing_limit)
                                                {
                                                    score /= bb_eclipsing_divisor;
                                                    #if _dbg_eclipsing_contacts
                                                    cout << bbr.pri_res->get_name() << " --- " << bbr.tert_res->get_name()
                                                        << " eclipsed by " << pocketres[h]->get_name() << " r=" << hr << endl;
                                                    #endif
                                                }
                                            }
                                        }
                                        if (hldist < lndist && hndist < lndist)
                                        {
                                            ha = pocketres[h]->get_nearest_atom_to_line(al->loc, an->loc);
                                            if (ha)
                                            {
                                                float hr = ha->loc.get_distance_to_line(al->loc, an->loc);
                                                if (hr < bb_eclipsing_limit)
                                                {
                                                    score /= bb_eclipsing_divisor;
                                                    #if _dbg_eclipsing_contacts
                                                    cout << bbr.sec_res->get_name() << " --- " << bbr.tert_res->get_name()
                                                        << " eclipsed by " << pocketres[h]->get_name() << " r=" << hr << endl;
                                                    #endif
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            #endif

                            #if _dbg_bb_scoring
                            cout << " ";
                            if (has_ionic) cout << "i";
                            if (has_neutral_polar) cout << "n";
                            if (has_pi) cout << "p";
                            cout << " ";
                            if (best_ionic) cout << "I";
                            if (best_neutral_polar) cout << "N";
                            if (best_pi) cout << "P";
                            #endif

                            if (score > best
                                && (has_ionic || !best_ionic)
                                && (has_neutral_polar || !best_neutral_polar)
                                && (has_pi || !best_pi))
                            {
                                // Assign by importance.
                                int i1, i2, i3, j1, j2, j3;
                                float iimp, kimp, mimp;
                                iimp = targets[i].importance(ligand);
                                kimp = targets[k].importance(ligand);
                                mimp = (m >= 0) ? targets[m].importance(ligand) : 0;

                                if (ijmc) iimp += 1000;
                                if (klmc) kimp += 1000;
                                if ((m >= 0) && mnmc) mimp += 1000;

                                /* cout << targets[i] << " " << iimp << " "
                                     << targets[k] << " " << kimp << " "
                                     << targets[m] << " " << mimp << endl; */

                                if (iimp > kimp && iimp > mimp)
                                {
                                    i1 = i; 
                                    j1 = j;

                                    if (kimp > mimp) { i2 = k; j2 = l; i3 = m; j3 = n; }
                                    else { i2 = m; j2 = n; i3 = k; j3 = l; }
                                }
                                else if (kimp > mimp)
                                {
                                    i1 = k;
                                    j1 = l;

                                    if (mimp > iimp) { i2 = m; j2 = n; i3 = i; j3 = j; }
                                    else { i2 = i; j2 = j; i3 = m; j3 = n; }
                                }
                                else
                                {
                                    i1 = m;
                                    j1 = n;

                                    if (kimp > iimp) { i2 = k; j2 = l; i3 = i; j3 = j; }
                                    else { i2 = i; j2 = j; i3 = k; j3 = l; }
                                }

                                if (m >= 0)
                                {
                                    float AB = targets[i1].barycenter().get_3d_distance(targets[i2].barycenter());
                                    float AC = targets[i1].barycenter().get_3d_distance(targets[i3].barycenter());

                                    float A_B_ = pocketres[j1]->distance_to(pocketres[j2]);
                                    float A_C_ = pocketres[j1]->distance_to(pocketres[j3]);

                                    // Nearby atoms should not pair with far away residues and vice versa.
                                    #if _dbg_bb_scoring
                                    if ((AC > AB && A_C_ < A_B_) || (AC < AB && A_C_ > A_B_))
                                    {
                                        cout << "X4" << endl;
                                        continue;
                                    }
                                    #else
                                    if (AC > AB && A_C_ < A_B_) continue;
                                    if (AC < AB && A_C_ > A_B_) continue;
                                    #endif

                                    if (AC > AB)
                                    {
                                        int swap = i2;
                                        i2 = i3;
                                        i3 = swap;
                                        swap = j2;
                                        j2 = j3;
                                        j3 = swap;
                                        // cout << "swapped" << endl;
                                    }
                                }

                                output->pri_res = pocketres[j1];
                                output->pri_tgt = &targets[i1];
                                output->sec_res = pocketres[j2];
                                output->sec_tgt = &targets[i2];
                                output->tert_res = (m >= 0) ? pocketres[j3] : nullptr;
                                output->tert_tgt = (m >= 0) ? &targets[i3] : nullptr;
                                output->probability = 1;
                                best = score;

                                if (has_ionic) best_ionic = true;
                                if (has_neutral_polar) best_neutral_polar = true;
                                if (has_pi) best_pi = true;

                                #if _dbg_bb_scoring
                                cout << " *** ";
                                cout << score;
                                #endif
                            }
                            
                            #if _dbg_bb_scoring
                            cout << endl;
                            #endif

                            if (m<0) break;
                        }
                    }
                }
            }
        }
    }

    if (!output->pri_res)
    {
        if (!override_target_compatibility)
        {
            override_target_compatibility = true;
            goto _new_attempt_no_tgt_comp;
        }
        cout << "No candidate BB pairings found!" << endl;
        throw 0x90ba1125;
    }

    // Make thiolates when binding to metals.
    if (allow_thiolation && (output->pri_res)->coordmtl && (output->pri_tgt)->single_atom)
    {
        int fam1 = (output->pri_tgt)->single_atom->get_family(), Z1 = (output->pri_tgt)->single_atom->Z;
        float chg1 = (output->pri_tgt)->charge();
        if (fam1 == CHALCOGEN && Z1 != 8 && chg1 <= 0)
        {
            Atom* H = (output->pri_tgt)->single_atom->is_bonded_to("H");
            if (H)
            {
                ligand->delete_atom(H);
                (output->pri_tgt)->single_atom->increment_charge(-1);
            }
        }
    }

    #if _dbg_bb_pairs
    if (output->pri_res && output->pri_tgt) cout << "Primary pairing: " << *output->pri_tgt << " to " << (output->pri_res)->get_name() << endl;
    if (*outbbr->sec_res && *outtgt2) cout << "Secondary pairing: " << **outtgt2 << " to " << (*outbbr->sec_res)->get_name() << endl;
    if (*outbbr->tert_res && *outtgt3) cout << "Tertiary pairing: " << **outtgt3 << " to " << (*outbbr->tert_res)->get_name() << endl;
    cout << endl;
    #endif
}

Point BestBindingResult::barycenter()
{
    Point pt(0,0,0);
    int divisor = 0;

    if (pri_res)
    {
        pt = pt.add(pri_res->get_reach_atom_location());
        divisor++;
    }
    if (sec_res)
    {
        pt = pt.add(sec_res->get_reach_atom_location());
        divisor++;
    }
    if (tert_res)
    {
        pt = pt.add(tert_res->get_reach_atom_location());
        divisor++;
    }

    if (divisor) pt.scale(pt.magnitude() / divisor);
    return pt;
}

#define bb_min_res_dist_minus_tgt_dist -1.0
#define bb_max_res_dist_minus_tgt_dist 4.6
#define bb_target_compatibility_boost 1.5

void Search::scan_protein(Protein *prot, Molecule *ligand, LigandTarget *targets, BestBindingResult *results, int max_results, Box limit)
{
    Pose putitback(ligand);
    int i, j, l, ii, jj, ll, m=0, mm, n = prot->get_end_resno(), nn;

    if (!limit.volume()) limit = Box(-Avogadro, Avogadro, -Avogadro, Avogadro, -Avogadro, Avogadro);

    for (nn=0; targets[nn].single_atom || targets[nn].conjgrp; nn++);

    Interaction result_energy[max_results+4];
    for (i=0; i<max_results; i++)
    {
        results[i].probability = 0;
        result_energy[i].clash = Avogadro;
    }

    AminoAcid* aaarr[SPHREACH_MAX+4];

    bool require_charged = false, require_metal = false, require_2x_charged = false;

    if (ligand->get_charge()) require_charged = true;
    if (prot->get_metals_count()) require_metal = true;

    Point protcen = prot->get_region_center(1, prot->get_end_resno());

    BestBindingResult bbr;
    for (ii=0; ii<nn; ii++)
    {
        bbr.pri_tgt = &targets[ii];
        float pri_imp = bbr.pri_tgt->importance(ligand);
        float tchg1 = bbr.pri_tgt->charge();

        if (require_charged && !bbr.pri_tgt->charge()) continue;

        for (i=1; i<=n; i++)
        {
            bbr.pri_res = prot->get_residue(i);
            if (!bbr.pri_res) continue;

            float rchg1 = bbr.pri_res->get_charge();
            if (require_charged && !rchg1) continue;
            if (require_metal && !bbr.pri_res->coordmtl) continue;
            if (rchg1 > 0 && tchg1 > 0) continue;

            Point CA1 = bbr.pri_res->get_CA_location();
            if (!limit.contains(CA1)) continue;

            #if 0
            if (i == 111) 
                cout << i 
                    << (target_compatibility(bbr.pri_res, bbr.pri_tgt) ? " comp" : " inc")
                    << endl;
            #else
            cout << "Scanning " << ii*n+i << "/" << n*nn << "..." << endl << "\033[A" << flush;
            #endif

            for (jj=0; jj<nn; jj++)
            {
                if (jj == ii) continue;
                bbr.sec_tgt = &targets[jj];

                float sec_imp = bbr.sec_tgt->importance(ligand);
                if (sec_imp > pri_imp) continue;
                float rtgt = bbr.pri_tgt->barycenter().get_3d_distance(bbr.sec_tgt->barycenter());
                float tchg2 = bbr.sec_tgt->charge();

                for (j=i+1; j<=n; j++)
                {
                    if (j == i) continue;
                    bbr.sec_res = prot->get_residue(j);
                    if (!bbr.sec_res) continue;

                    Point CA2 = bbr.sec_res->get_CA_location();
                    if (!limit.contains(CA2)) continue;

                    Atom* ext1 = bbr.pri_res->get_reach_atom();
                    Atom* ext2 = bbr.sec_res->get_reach_atom();
                    float rres = (ext1 && ext2) ? ext1->distance_to(ext2) : bbr.pri_res->distance_to(bbr.sec_res);
                    if (rres > bbr.pri_res->get_CA_location().get_3d_distance(bbr.sec_res->get_CA_location())) continue;
                    if (rres < rtgt + bb_min_res_dist_minus_tgt_dist) continue;
                    if (rres > rtgt + bb_max_res_dist_minus_tgt_dist) continue;
                    float rchg2 = bbr.sec_res->get_charge();

                    if (rchg1 && tchg1 && rchg2 && tchg2)
                        require_2x_charged = true;
                    else if (require_2x_charged) continue;

                    bbr.tert_res = nullptr;
                    bbr.tert_tgt = nullptr;

                    for (ll=-1; ll<nn; ll++)
                    {
                        if (ll == ii) continue;
                        if (ll == jj) continue;
                        bbr.tert_tgt = (ll>=0) ? &targets[ll] : nullptr;
                        float tert_imp = (ll>=0) ? bbr.tert_tgt->importance(ligand) : 0;
                        if (ll>=0 && tert_imp > sec_imp) continue;

                        float rtgt2, rtgt3;
                        if (ll>=0)
                        {
                            rtgt2 = bbr.pri_tgt->barycenter().get_3d_distance(bbr.tert_tgt->barycenter());
                            rtgt3 = bbr.tert_tgt->barycenter().get_3d_distance(bbr.sec_tgt->barycenter());
                        }

                        // int lmin = (ll<0) ? -1 : 0;
                        // int lmax = (ll<0) ? -1 : n-1;
                        for (l=1; l<=n; l++)
                        {
                            // if (l>=0 && l <= j) continue;
                            if (l==i || l==j) continue;

                            bbr.tert_res = (l<0) ? nullptr : prot->get_residue(l);
                            if (l>=0 && !bbr.tert_res) continue;

                            Point CA3 = bbr.tert_res ? bbr.tert_res->get_CA_location() : Point(0,0,0);
                            if (bbr.tert_res && !limit.contains(CA3)) continue;

                            if (ll>=0)
                            {
                                float rres2 = bbr.pri_res->distance_to(bbr.tert_res);
                                float rres3 = bbr.tert_res->distance_to(bbr.sec_res);
                                if (rres2 < rtgt2 + bb_min_res_dist_minus_tgt_dist) continue;
                                if (rres2 > rtgt2 + bb_max_res_dist_minus_tgt_dist) continue;
                                if (rres3 < rtgt3 + bb_min_res_dist_minus_tgt_dist) continue;
                                if (rres3 > rtgt3 + bb_max_res_dist_minus_tgt_dist) continue;
                            }
                            else
                            {
                                if (bbr.tert_res->hydrophilicity() >= hydrophilicity_cutoff) continue;
                                if (bbr.tert_res->pi_stackability() > 0.1) continue;
                            }

                            Point avg = bbr.pri_res->get_CA_location();
                            avg = avg.add(bbr.sec_res->get_CA_location());
                            if (bbr.tert_res) avg = avg.add(bbr.tert_res->get_CA_location());
                            /*else
                            {
                                Point avg2 = avg;
                                avg2.multiply(0.5);
                                Point lonely = prot->find_loneliest_point(avg2, Point(4,4,4));
                                avg = avg.add(lonely);
                            }*/
                            avg.multiply(1.0 / (bbr.tert_res ? 3 : 2));

                            float score = bbr.score(avg);

                            // Prefer compatible targets
                            if (target_compatibility(bbr.pri_res, bbr.pri_tgt)) score *= bb_target_compatibility_boost;
                            if (target_compatibility(bbr.sec_res, bbr.sec_tgt)) score *= bb_target_compatibility_boost;
                            if (ll >= 0 && target_compatibility(bbr.tert_res, bbr.tert_tgt)) score *= bb_target_compatibility_boost;

                            // Prefer pockets near protein center
                            float r = protcen.get_3d_distance(avg);
                            score /= sqrt(r+1);

                            float kchg = bbr.sec_tgt->charge(), lchg = bbr.sec_res->get_charge();
                            if (kchg && lchg && sgn(kchg) == -sgn(lchg)) score *= bb_target_compatibility_boost;

                            for (m=0; m<max_results; m++)
                            {
                                if (score > 0.8*results[m].probability)
                                {
                                    align_targets(ligand, avg, &bbr, 1);
                                    if (ll<0)
                                    {
                                        SCoord nudge = loneliest.subtract(ligand->get_barycenter());
                                        nudge.r = 1;
                                        ligand->move(nudge);
                                    }

                                    int exclres[4];
                                    exclres[0] = bbr.pri_res->get_residue_no();
                                    exclres[1] = bbr.sec_res->get_residue_no();
                                    exclres[2] = bbr.tert_res ? bbr.tert_res->get_residue_no() : 0;
                                    exclres[3] = 0;
                                    prot->get_residues_can_clash_ligand(aaarr, ligand, avg, Point(7,7,7), nullptr, false, exclres);

                                    ligand->movability = MOV_ALL;
                                    ligand->quick_conform((Molecule**)aaarr, 5);

                                    Interaction e = ligand->get_intermol_binding((Molecule**)aaarr);

                                    if (e.clash >= 100)
                                    {
                                        int it;
                                        Interaction f;
                                        float nudgeamt = 0.03;
                                        for (it=0; it<10; it++)
                                        {
                                            if (!ligand->clash1 || !ligand->clash2) break;
                                            SCoord nudge = ligand->clash1->loc.subtract(ligand->clash2->loc);
                                            nudge.r = nudgeamt;
                                            ligand->move(nudge);
                                            f = ligand->get_intermol_binding((Molecule**)aaarr);

                                            // if (ligand->clash2->residue == 115) cout << "                     " << f.clash << endl;

                                            if (f.improved(e)) e = f;
                                            else
                                            {
                                                nudge.r *= -1;
                                                ligand->move(nudge);
                                                nudgeamt *= 0.8;
                                            }
                                            if (!e.clash) break;
                                        }
                                    }

                                    /* if (ligand->clash1 && ligand->clash2 && ligand->clash1->Z == 1 && ligand->clash2->Z == 1)
                                        e.clash = 0; */

                                    #if 0
                                    if ((bbr.pri_res->get_residue_no() == 201 || bbr.pri_res->get_residue_no() == 111)
                                        && (bbr.sec_res->get_residue_no() == 201 || bbr.sec_res->get_residue_no() == 111)) 
                                    {
                                        if (1) // e.clash > 1000)
                                        {
                                            cout << "                                    ";
                                            if (ligand->clash1) cout << ligand->clash1->residue << ":" << ligand->clash1->name;
                                            cout << " X ";
                                            if (ligand->clash2) cout << ligand->clash2->residue << ":" << ligand->clash2->name;
                                            cout << " " << e.clash;
                                            cout << endl;
                                        }

                                        cout << "                                    ";
                                        cout << *(bbr.pri_tgt) << "~" << bbr.pri_res->get_name() << ", "
                                            << *(bbr.sec_tgt) << "~" << bbr.sec_res->get_name();
                                        if (bbr.tert_res && bbr.tert_tgt) cout << ", " << *(bbr.tert_tgt) << "~" << bbr.tert_res->get_name();
                                        else if (bbr.tert_res) cout << ", (null)~" << bbr.tert_res->get_name();
                                        cout << " @ " << rres << "/" << rtgt;
                                        cout << " = " << score << "/" << e.summed() << endl << flush;
                                    }
                                    #endif

                                    if ((e.summed() < 0 && score > results[m].probability) || e.clash < result_energy[m].clash)
                                    {
                                        #if 0
                                        cout << *(bbr.pri_tgt) << "~" << bbr.pri_res->get_name() << ", "
                                            << *(bbr.sec_tgt) << "~" << bbr.sec_res->get_name();
                                        if (bbr.tert_res) cout << ", " << *(bbr.tert_tgt) << "~" << bbr.tert_res->get_name();
                                        cout << " @ " << score << endl << flush;
                                        #endif

                                        if (results[m].pri_res != bbr.pri_res || results[mm].sec_res != bbr.sec_res)
                                        {
                                            int mmstart;

                                            for (mmstart=m+1; mmstart<max_results; mmstart++)
                                                if (results[mmstart].pri_res == bbr.pri_res && results[mmstart].sec_res == bbr.sec_res) break;
                                            if (mmstart >= max_results) mmstart = max_results-1;
    
                                            for (mm=mmstart; mm>m; mm--)
                                            {
                                                results[mm] = results[mm-1];
                                                result_energy[mm] = result_energy[mm-1];

                                                if (require_2x_charged && results[mm].sec_res)
                                                {
                                                    if (!results[mm].sec_res->get_charge()) results[mm].probability = 0;
                                                }
                                            }
                                        }

                                        results[m] = bbr;
                                        results[m].probability = score;
                                        result_energy[m] = e;

                                        break;
                                    }
                                    else if (results[m].pri_res == bbr.pri_res && results[mm].sec_res == bbr.sec_res) break;
                                }
                                else if (results[m].pri_res == bbr.pri_res && results[mm].sec_res == bbr.sec_res) break;
                            }
                        }
                    }
                }
            }
        }
    }

    float prob_total = 0;
    for (m=0; m<max_results; m++)
    {
        // if (!m) results[m].probability = 1;
        // else results[m].probability = equilibrium(-result_energy[0].attractive, -result_energy[m].attractive);
        prob_total += results[m].probability;
    }

    for (m=0; m<max_results; m++) results[m].probability /= prob_total;

    // Now sort the results. Bubble sort will be adequate.
    for (i=0; i<max_results; i++)
    {
        for (m=1; m<max_results; m++)
        {
            if (results[m-1].probability < results[m].probability)
            {
                BestBindingResult swap = results[m-1];
                results[m-1] = results[m];
                results[m] = swap;
            }
        }
    }

    putitback.restore_state(ligand);
}

void Search::align_targets(Molecule *ligand, Point pocketcen, BestBindingResult* bbr, float amt)
{
    ligand->movability = MOV_ALL;

    int i;
    Bond** bb;
    #if bb_stochastic_ligflex
    bb = ligand->get_rotatable_bonds();
    for (i=0; bb[i]; i++)
    {
        if (frand(0,1) < bb_stochastic/3) bb[i]->rotate(frand(-hexagonal/2, hexagonal/2));
    }
    #endif

    if (amt >= 0.9999) ligand->recenter(pocketcen);
    else
    {
        Point cen = ligand->get_barycenter();
        SCoord rel = pocketcen.subtract(cen);
        rel.r *= amt;
        ligand->move(rel);
    }

    // Align primary.
    Atom* aaaa = bbr->pri_res->coordmtl ? bbr->pri_res->coordmtl : bbr->pri_res->get_nearest_atom(pocketcen);       // amino acid alignment atom
    #if bb_find_empty_space
    Point nearby = bbr->protein ? bbr->protein->find_loneliest_point(aaaa->loc, Point(bb_lonely_radius,bb_lonely_radius,bb_lonely_radius)) : aaaa->loc;
    Rotation rot = align_points_3d(bbr->pri_tgt->barycenter(), nearby.randomize(bb_stochastic_A), pocketcen);
    #else
    Rotation rot = align_points_3d(bbr->pri_tgt->barycenter(), aaaa->loc, pocketcen);
    #endif
    rot.a *= amt;
    if (!aaaa->vanished) ligand->rotate(rot, pocketcen);

    // cout << "Rotated " << *bbr->pri_tgt << " " << (rot.a*fiftyseven) << "deg to face " << *bbr->pri_res << endl;

    // cout << endl << endl; return;

    // Scooch.
    SCoord rel;
    #if enable_bb_scooch
    aaaa = bbr->pri_res->get_nearest_atom(pocketcen);
    if (!aaaa->vanished)
    {
        rel = aaaa->loc.subtract(bbr->pri_tgt->barycenter().randomize(bb_stochastic_A));
        rel.r -= 2.5;
        if (rel.r > 0) 
        {
            rel.r *= amt;
            ligand->move(rel);
        }
    }
    #endif

    // return;

    // Align secondary.
    if (!bbr->sec_res || !bbr->sec_tgt) return;

    #if bb_find_empty_space
    // aaaa = bbr->pri_res->get_nearest_atom(pocketcen);
    Point ref = bbr->pri_tgt->barycenter();
    nearby = bbr->sec_res->get_nearest_atom(pocketcen)->loc;
    #if _dbg_bb_contact_lonely
    cout << "Secondary residue location: " << nearby;
    #endif
    if (bbr->protein)
    {
        nearby = bbr->protein->find_loneliest_point(nearby, Point(bb_lonely_radius,bb_lonely_radius,bb_lonely_radius));
        #if _dbg_bb_contact_lonely
        cout << " moved to " << nearby;
        #endif
    }
    #if _dbg_bb_contact_lonely
    cout << endl;
    #endif
    #else
    aaaa = bbr->sec_res->get_nearest_atom(pocketcen);
    Point ref = bbr->pri_tgt->barycenter();
    #endif
    if (!aaaa->vanished)
    {
        #if bb_find_empty_space
        rot = align_points_3d(bbr->sec_tgt->barycenter(), nearby.randomize(bb_stochastic_A), ref);
        #else
        rot = align_points_3d(bbr->sec_tgt->barycenter(), aaaa->loc.randomize(bb_stochastic_A), ref);
        #endif
        rot.a *= amt;
        ligand->rotate(rot, ref);
    }

    // Rotisserie align tertiary.
    if (!bbr->tert_res || !bbr->tert_tgt) return;
    ref = bbr->pri_tgt->barycenter();
    SCoord axis = bbr->sec_tgt->barycenter().subtract(ref);
    float theta = find_angle_along_vector(bbr->tert_tgt->barycenter(), 
        bbr->tert_res->get_nearest_atom(pocketcen)->loc.randomize(bb_stochastic_A), ref, axis);
    LocatedVector lv = axis;
    lv.origin = ref;
    if (!aaaa->vanished) ligand->rotate(lv, theta*amt);

    // Tertiary scooch
    if (bbr->sec_tgt->polarity() < hydrophilicity_cutoff && bbr->tert_tgt->polarity() >= hydrophilicity_cutoff)
    {
        aaaa = bbr->tert_res->get_nearest_atom(pocketcen);
        rel = aaaa->loc.subtract(bbr->tert_tgt->barycenter()).randomize(bb_stochastic_A);
        rel.r -= 2.5;
        if (rel.r > 0 && !aaaa->vanished) 
        {
            rel.r *= 0.5 * amt;
            ligand->move(rel);
        }
    }
}

bool Search::target_compatibility(float chg1, float chg2, float pol1, float pol2, int pi1, int pi2, bool hba1, bool hba2, bool hbd1, bool hbd2)
{
    if (chg1 && sgn(chg1) == sgn(chg2)) return false;
    if (chg1 && sgn(chg1) == -sgn(chg2)) return true;
    if (!chg1 && !chg2 && pi1 > 5 && pi2 < 2) return false;
    if (!chg1 && !chg2 && pi2 > 5 && pi1 < 2) return false;
    if ((pol1 >= hydrophilicity_cutoff) != (pol2 >= hydrophilicity_cutoff) && (pi1 < 5 || pi2 < 5)) return false;
    if (pol1 >= hydrophilicity_cutoff && pol2 >= hydrophilicity_cutoff)
    {
        if (hba1 && hbd2) return true;
        if (hbd1 && hba2) return true;
        return false;
    }
    return true;
}

bool Search::target_compatibility(AminoAcid *aa, LigandTarget *lt)
{
    return target_compatibility(aa->get_charge(), lt->charge(),
        fabs(aa->hydrophilicity()), fabs(lt->polarity()),
        aa->has_pi_atoms()/aa->get_heavy_atom_count(), lt->is_pi() ? 1 : 0,
        aa->has_hbond_acceptors(), lt->has_hb_acceptors(),
        aa->has_hbond_donors(), lt->has_hb_donors()
    );
}

float LigandTarget::charge()
{
    if (single_atom) return single_atom->get_orig_charge();
    else if (conjgrp) return conjgrp->get_net_charge();
    return 0.0f;
}

float LigandTarget::polarity()
{
    if (single_atom) return fabs(single_atom->is_polar());
    else if (conjgrp) return conjgrp->get_sum_polarity();
    return 0.0f;
}

Point LigandTarget::barycenter()
{
    if (single_atom) return single_atom->loc;
    else if (conjgrp) return conjgrp->get_barycenter();
    return Point(0,0,0);
}

bool LigandTarget::is_pi()
{
    if (single_atom) return single_atom->is_pi();
    else if (conjgrp) return true;
    return false;
}

bool LigandTarget::has_hb_acceptors()
{
    if (single_atom) return (single_atom->is_polar() <= -hydrophilicity_cutoff || single_atom->get_num_lone_pairs() > 0);
    else if (conjgrp) return conjgrp->has_hb_acceptors();
    return false;
}

bool LigandTarget::has_hb_donors()
{
    if (single_atom)
    {
        int afam = single_atom->get_family();
        if (afam <= TETREL) return false;
        Atom* H = single_atom->is_bonded_to("H");
        if (H) return true; // (H->is_polar() >= hydrophilicity_cutoff);
        else return false;
    }
    else if (conjgrp) return conjgrp->has_hb_donors();
    return false;
}

float LigandTarget::importance(Molecule *mol)
{
    Point bary = mol->get_barycenter();
    float result = 0;
    if (single_atom)
    {
        result += 60.0 * fabs(single_atom->get_charge());
        result += 35.0 * fabs(single_atom->is_polar());
        if (single_atom->is_pi()) result += 12;
        result += single_atom->num_bonded_to("H");
        result += bary.get_3d_distance(single_atom->loc);
    }
    else if (conjgrp)
    {
        result += 60.0 * fabs(conjgrp->get_net_charge());
        result += 35.0 * fabs(conjgrp->get_sum_polarity());
        result += 12;
        result += bary.get_3d_distance(conjgrp->get_barycenter());
    }
    return result;
}

bool LigandTarget::contains(LigandTarget *lt)
{
    if (conjgrp && lt->single_atom) return conjgrp->contains(lt->single_atom);
    else if (single_atom && lt->single_atom) return single_atom == lt->single_atom;
    return false;
}

bool LigandTarget::contains(Atom *a)
{
    if (single_atom && single_atom == a) return true;
    if (conjgrp && conjgrp->contains(a)) return true;
    return false;
}

std::string LigandTarget::to_std_string()
{
    if (single_atom) return std::string(single_atom->name);
    else if (conjgrp) return conjgrp->to_std_string();
    else return std::string("(null)");
}

std::ostream& operator<<(std::ostream& os, const LigandTarget& lt)
{
    if (lt.single_atom) os << lt.single_atom->name;
    else if (lt.conjgrp) os << *lt.conjgrp;
    else os << "(null)";

    return os;
}

std::ostream &operator<<(std::ostream &os, const BestBindingResult &bbr)
{
    if (!bbr.pri_res || !bbr.pri_tgt)
    {
        os << "(empty)" << endl;
        return os;
    }
    os << "Pri: " << *(bbr.pri_tgt) << "..." << *(bbr.pri_res) << endl;
    if (!bbr.sec_res || !bbr.sec_tgt) return os;
    os << "Sec: " << *(bbr.sec_tgt) << "..." << *(bbr.sec_res) << endl;
    if (!bbr.tert_res || !bbr.tert_tgt) return os;
    os << "Ter: " << *(bbr.tert_tgt) << "..." << *(bbr.tert_res) << endl;
    return os;
}

float BestBindingResult::score(Point ligcen, Cavity* container)
{
    float lchg=0, lpol=0, lcba;
    int lpi;
    bool lhba, lhbd, klmc;
    Atom* lmtl;
    float mchg, mpol;
    int mfam, mZ, mpi;
    Point mcen;
    bool mhba, mhbd;
    float nchg, npol, ncba;
    int npi;
    bool nhba, nhbd, mnmc;
    Atom* nmtl;

    float ichg = pri_tgt->charge();
    float ipol = pri_tgt->polarity();
    if (ipol < hydrophilicity_cutoff) ipol = 0;
    int ifam = pri_tgt->single_atom ? pri_tgt->single_atom->get_family() : -1;
    int iZ = pri_tgt->single_atom ? pri_tgt->single_atom->Z : -1;
    if (ifam == PNICTOGEN) ipol += hydrophilicity_cutoff*2;
    int ipi = pri_tgt->conjgrp ? pri_tgt->conjgrp->num_heavy_atoms : 0;
    bool ihba = pri_tgt->has_hb_acceptors();
    bool ihbd = pri_tgt->has_hb_donors();
    Point icen = pri_tgt->barycenter();

    float jchg = pri_res->get_charge();
    float jpol = pri_res->hydrophilicity();
    if (jpol < hydrophilicity_cutoff) jpol = 0;
    int jpi = pri_res->has_pi_atoms(false);
    float jcba = pri_res->CB_angle(loneliest);
    bool jhba = pri_res->has_hbond_acceptors();
    bool jhbd = pri_res->has_hbond_donors();
    Atom* jmtl = pri_res->coordmtl;

    bool ijmc = jmtl && (ifam == CHALCOGEN || ifam == PNICTOGEN) && iZ != 8 && ichg <= 0;

    float kchg=0, kpol=0;
    int kfam, kZ, kpi;
    Point kcen;
    bool khba, khbd;

    if (sec_tgt)
    {
        kchg = sec_tgt->charge();
        kpol = sec_tgt->polarity();
        if (kpol < hydrophilicity_cutoff) kpol = 0;
        kfam = sec_tgt->single_atom ? sec_tgt->single_atom->get_family() : -1;
        kZ = sec_tgt->single_atom ? sec_tgt->single_atom->Z : -1;
        if (kfam == PNICTOGEN) kpol += hydrophilicity_cutoff*2;
        kpi = sec_tgt->conjgrp ? sec_tgt->conjgrp->num_heavy_atoms : 0;
        kcen = sec_tgt->barycenter();
        khba = sec_tgt->has_hb_acceptors();
        khbd = sec_tgt->has_hb_donors();
    }

    if (sec_res)
    {
        lchg = sec_res->get_charge();
        lpol = sec_res->hydrophilicity();
        if (lpol < hydrophilicity_cutoff) lpol = 0;
        lpi = sec_res->has_pi_atoms(false);
        lcba = sec_res->CB_angle(loneliest);
        lhba = sec_res->has_hbond_acceptors();
        lhbd = sec_res->has_hbond_donors();
        lmtl = sec_res->coordmtl;
    }

    klmc = lmtl && (kfam == CHALCOGEN || kfam == PNICTOGEN) && kZ != 8 && kchg <= 0;

    if (tert_res && tert_tgt)
    {
        mchg = tert_tgt->charge();
        mpol = tert_tgt->polarity();
        if (mpol < hydrophilicity_cutoff) mpol = 0;
        mfam = tert_tgt->single_atom ? tert_tgt->single_atom->get_family() : -1;
        mZ = tert_tgt->single_atom ? tert_tgt->single_atom->Z : -1;
        if (mfam == PNICTOGEN) mpol += hydrophilicity_cutoff*2;
        mpi = tert_tgt->conjgrp ? tert_tgt->conjgrp->num_heavy_atoms : 0;
        mcen = tert_tgt->barycenter();
        mhba = tert_tgt->has_hb_acceptors();
        mhbd = tert_tgt->has_hb_donors();
        nchg = tert_res->get_charge();
        npol = tert_res->hydrophilicity();
        npi = tert_res->has_pi_atoms(false);
        ncba = tert_res->CB_angle(loneliest);
        nhba = tert_res->has_hbond_acceptors();
        nhbd = tert_res->has_hbond_donors();
        nmtl = tert_res->coordmtl;

        mnmc = nmtl && (mfam == CHALCOGEN || mfam == PNICTOGEN) && mZ != 8 && mchg <= 0;
    }
    
    #if _dbg_bb_scoring
    cout << *pri_tgt << ">" << pri_res->get_name() << ","
        << *sec_tgt << ">" << sec_res->get_name();
    if (tert_res && tert_tgt) cout << "," << *tert_tgt << ">" << tert_res->get_name();
    cout << " ";
    #endif

    // Score based on strength of contacts.
    float score = 0;
    if (ijmc) score += 200;
    if (ichg && jchg && sgn(ichg) == -sgn(jchg)) score += 60; 
    if (ipol && jpol) 
    {
        score += 35;
        if (pri_res->is_amide() && ichg <= 0) score += 10;
    }
    if (!ipol && !jpol) score += 0.2 * sqrt(1.0 + icen.get_3d_distance(ligcen));
    if ((!ipol || !jpol) && ipi && jpi) score += 2.0*ipi*jpi; 

    if (sec_res && sec_tgt)
    {
        if (klmc) score += 200;
        if (kchg && lchg && sgn(kchg) == -sgn(lchg)) score += 60; 
        if (lpol && kpol) 
        {
            score += 35;
            if (sec_res->is_amide() && kchg <= 0) score += 10;
        } 
        if (!kpol && !lpol) score += 0.2 * sqrt(1.0 + kcen.get_3d_distance(ligcen));
        if ((!lpol || !kpol) && kpi && lpi) score += 2.0*kpi*lpi;
    }

    if (tert_res && tert_tgt)
    {
        if (mnmc) score += 200;
        if (mchg && nchg && sgn(mchg) == -sgn(nchg)) score += 60; 
        if (mpol && npol) 
        {
            score += 35;
            if (tert_res->is_amide() && mchg <= 0) score += 10;
        }
        if (!mpol && !npol) score += 0.2 * sqrt(1.0 + mcen.get_3d_distance(ligcen));
        if ((!mpol || !npol) && mpi && npi) score += 2.0*mpi*npi;
    }

    #if _dbg_bb_scoring
    cout << "a:" << score << " ";
    #endif

    // Prefer aminos that point their sidechains toward the pocket center.
    score *= 1.0+cos(jcba);
    if (sec_res && sec_tgt) score *= 1.0+cos(lcba);
    if (tert_res && tert_tgt) score *= 0.5+0.5*cos(ncba);

    #if _dbg_bb_scoring
    cout << "b:" << score << " ";
    #endif

    // Match contacts by their relative spacings.
    float ikdist, jldist, imdist, jndist, kmdist, lndist;

    if (sec_res && sec_tgt)
    {
        if (kpol && lpol)
        {
            ikdist = pri_tgt->barycenter().get_3d_distance(sec_tgt->barycenter());
            jldist = fmax(1, pri_res->distance_to(sec_res) - bb_pocket_res_extra_spacing);
            if (jldist > ikdist) jldist -= fmin(jldist-ikdist, bb_pocket_res_spacing_allowance); 

            score *= (fmin(ikdist, jldist) / fmax(ikdist, jldist));
        }

        if (tert_res && tert_tgt && mpol && npol)
        {
            imdist = pri_tgt->barycenter().get_3d_distance(tert_tgt->barycenter());
            jndist = fmax(1, pri_res->distance_to(tert_res) - bb_pocket_res_extra_spacing);
            if (jndist > imdist) jndist -= fmin(jndist-imdist, bb_pocket_res_spacing_allowance);

            score *= (fmin(imdist, jndist) / fmax(imdist, jndist));

            kmdist = sec_tgt->barycenter().get_3d_distance(tert_tgt->barycenter());
            lndist = fmax(1, sec_res->distance_to(tert_res) - bb_pocket_res_extra_spacing);
            if (lndist > kmdist) lndist -= fmin(lndist-kmdist, bb_pocket_res_spacing_allowance);

            score *= (fmin(kmdist, lndist) / fmax(kmdist, lndist));
        }
        else imdist = jndist = kmdist = lndist = 1;
    }

    #if _dbg_bb_scoring
    cout << "c:" << score << " ";
    #endif

    // Avoid contacts all in a line.
    if (sec_res && tert_res)
    {
        Atom* ra1 = pri_res->get_nearest_atom(ligcen);
        Atom* ra2 = sec_res->get_nearest_atom(ligcen);
        Atom* ra3 = tert_res->get_nearest_atom(ligcen);
        float theta = fmin(find_3d_angle(ra1->loc, ra2->loc, ra3->loc),
            fmin(find_3d_angle(ra1->loc, ra3->loc, ra2->loc), find_3d_angle(ra3->loc, ra2->loc, ra1->loc)));
        score *= fabs(sin(theta*2));
        int res1 = pri_res->get_residue_no(), res2 = sec_res->get_residue_no(), res3 = tert_res->get_residue_no();
        /*if (max(res1, max(res2, res3)) == 262 && min(res1, min(res2, res3)) == 158)
            cout << "good " << res1 << ":" << ra1->name << "|" << res2 << ":" << ra2->name 
                << "|" << res3 << ":" << ra3->name << " " << theta*fiftyseven << endl;
        if (max(res1, max(res2, res3)) == 262 && min(res1, min(res2, res3)) == 184)
            cout << "not good " << res1 << ":" << ra1->name << "|" << res2 << ":" << ra2->name 
                << "|" << res3 << ":" << ra3->name << " " << theta*fiftyseven << endl;*/
    }

    // Prefer contacts that span the entire molecule.
    if (sec_tgt) score *= pow(ikdist + imdist + kmdist, 2) / 1000;

    #if _dbg_bb_scoring
    cout << "d:" << score << " ";
    #endif

    // Prefer results that include hydrogen bonds and ionic bonds.
    if (ipol && jpol) score *= 2;
    if (ichg && jchg) score *= 5;
    if (kpol && lpol) score *= 1.5;
    if (kchg && lchg) score *= 3;
    if (mpol && npol) score *= 1.25;

    // If a containing cavity is set, prefer results that are a good fit for the cavity.
    if (container)
    {
        float lip = container->molecule_inside_pocket(ligand, true) * (1.0+bb_stochastic);
        // cout << lip << endl;
        score *= (1.5+lip);
    }

    // Prefer contacts with residues whose positions average near the pocket center.
    /*Point rismediu[3];
    rismediu[0] = pri_res->get_CA_location();
    rismediu[1] = sec_res->get_CA_location();
    if (tert_res && tert_tgt) rismediu[2] = tert_res->get_CA_location();
    Point avg = average_of_points(rismediu, (tert_res && tert_tgt) ? 3 : 2);
    score /= (1.0+loneliest.get_3d_distance(avg));*/

    #if _dbg_bb_scoring
    cout << "e:" << score << " ";
    #endif

    // Effect of priority.
    if (pri_res->priority) score *= priority_weight_group;
    if (sec_res && sec_tgt && sec_res->priority) score *= priority_weight_group;
    if (tert_res && tert_tgt && tert_res->priority) score *= priority_weight_group;

    #if _dbg_bb_scoring
    cout << "f:" << score << " ";
    #endif

    return score;
}
