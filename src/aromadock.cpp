#include <cstring>
#include <iostream>
#include <stdio.h>
#include <string>
#include <regex>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sstream>
#include <algorithm>
#include <unistd.h>
#include "classes/dynamic.h"
#include "classes/search.h"
#include "classes/scoring.h"
#include "classes/cavity.h"
#include "classes/soft.h"
#include "classes/appear.h"
#include "classes/reshape.h"
#include "classes/progress.h"

using namespace std;

std::string get_fttl_no_path_or_ext(const char* filename)
{
    char buffer[256];
    strcpy(buffer, filename);
    int i, n = strlen(buffer);
    for (i=n; i>0; i--)
    {
        if (buffer[i] == '.')
        {
            buffer[i] = 0;
            break;
        }
    }
    for (; i>0; i--)
    {
        if (buffer[i] == '/')
        {
            std::string result = &buffer[i+1];
            return result;
        }
    }

    std::string result = filename;
    return result;
}

char* get_file_ext(char* filename)
{
    int i = strlen(filename)-1;
    for (; i>=0; i--)
    {
        if (filename[i] == '.') return &filename[i+1];
    }
    return 0;
}

bool output_each_iter = false;
bool output_something_even_if_it_is_wrong = false;
bool progressbar = false;
Progressbar progb;
std::string itersfname;
int movie_offset = 0;
char configfname[256];
char protfname[256];
char protstrand = '\0';
char protafname[256];
char tplfname[256];
char tplrfnam[256];
char ligfname[256];
char smiles[256];
char outfname[256];
char cvtyfname[256];
Point pocketcen;
std::ofstream *output = NULL;

std::vector<std::string> CEN_buf;
int cenbuf_idx = 0;
std::vector<std::string> pathstrs;
std::vector<std::string> states;

Cavity cvtys[256];
int ncvtys = 0;

bool configset=false, protset=false, tplset=false, tprfset=false, ligset=false, ligcmd=false, smset = false, smcmd = false, pktset=false;
bool success_sofar = false;

bool out_per_res_e = true;
bool out_per_btyp_e = true;
float out_itemized_e_cutoff = 0.01;
bool out_lig_int_e = true;
bool out_bb_pairs = false;
bool out_lig_pol_sat = false;
bool out_prox = false;
bool out_pro_clash = false;
bool display_clash_atoms = false;
bool out_mc = false;
bool out_vdw_repuls = false;
bool out_pdbdat_lig = true;
bool out_pdbdat_res = true;

Protein* protein;
Protein* ptemplt;
Protein* ptplref;
Protein* metald_prot = nullptr;
int seql = 0;
int addl_resno[256];
const Region* regions;
Vector region_clashes[85][3];
Molecule* ligand;
std::vector<std::string> isomers;
int copylig = 0;
bool estimate_copylig = false;
int isono = 0;
Molecule** waters = nullptr;
Molecule** owaters = nullptr;
Pose* cfmol_known_good[SPHREACH_MAX+100];
Point ligcen_target;
Vector path[256];
int pathnodes = 0;				// The pocketcen is the initial node.
int poses = 10;
int iters = 50;
int maxh2o = 0;
int omaxh2o = 0;
bool wet_contact = false;       // If true, place a water molecule between the strongest affinity residue and ligand atom pair.
bool flex = true;
float kJmol_cutoff = -0.01;
int pose=1, nodeno=0, iter=0;
bool kcal = false;
float drift = initial_drift;
Molecule** gcfmols = NULL;
int activation_node = -1;		// Default is never.
int found_poses = 0;
int triesleft = 0;				// Default is no retry.
bool echo_progress = false;
bool hydrogenate_pdb = false;
std::string temp_pdb_file = "";
int pid = getpid();
bool append_pdb = false;
bool do_output_colors = false;

LigandTarget* g_ligtargs;
BestBindingResult g_bbr[16];

MCoord mtlcoords[16];
int nmtlcoords;

bool softdock = false;
float softness = 0;

SoftRegion softrgns[32];
int nsoftrgn = 0;
int softrgn_allowed[32];
int nsoftrgna = 0;
ResiduePlaceholder soft_contact_a[256], soft_contact_b[256];
int n_soft_contact = 0;
ResiduePlaceholder soft_nodel_start[32], soft_nodel_end[32];
int nsoftnodel = 0;
char splash[16384];

AminoAcid*** reaches_spheroid = nullptr;
int sphres = 0;

std::string origbuff = "";
std::string optsecho = "";

// Switch to enable "best-binding" algorithm rather than "tumble spheres" algorithm.
PoseSearchType pdpst = default_search_algorithm;
std::string copyfrom_filename;
Pose preconforms[50];
int npreconf;
char copyfrom_ligname[5] = {0,0,0,0,0};
int copyfrom_resno = -1;
Bond retain_bindings[4];
std::vector<int> center_resnos;
std::vector<int> priority_resnos;

Atom* pivotal_hbond_aaa = nullptr;
Atom* pivotal_hbond_la = nullptr;
float pivotal_hbond_r = 0;

Pose pullaway_undo;
float last_ttl_bb_dist = 0;

std::vector<ResiduePlaceholder> required_contacts;
std::vector<std::string> bridges;
std::vector<std::string> atomto;
std::string nfwarned = "";
std::string outpdb;
std::string auditfn = "";
int outpdb_poses = 0;
bool outdisq = false;
bool warn_absent_residue = true;
std::string outdisqrs = "";

std::vector<int>flexible_resnos;
std::vector<ResiduePlaceholder>forced_flexible_resnos;
std::vector<ResiduePlaceholder>forced_static_resnos;

Appear appears[32];
int nappears = 0;
int last_appear_iter = 0;

#if reuse_best_pose
Pose ligand_best_pose[64];
Pose* aa_best_pose = nullptr;
float best_pose_energy;

Atom *best_ligand_stays, *best_ligand_stays_other;
#endif

ICHelixGroup hg;

#if _dummy_atoms_for_debug
std::vector<Atom> dummies;
#endif


void append_dummy(Point pt)
{
    #if _dummy_atoms_for_debug
    Atom a("He");
    a.move(pt);

    a.name = new char[8];
    int i=dummies.size()+1;
    sprintf(a.name, "HE%i", i);

    strcpy(a.aa3let, "DMY");

    dummies.push_back(a);
    #endif
}

void delete_water(Molecule* mol)
{
    if (!waters) return;
    int i, j;
    for (i=0; waters[i]; i++)
    {
        if (waters[i] == mol)
        {
            for (j=i+1; waters[j]; j++)
            {
                waters[j-1] = waters[j];
            }
            waters[j-1] = nullptr;
            maxh2o--;

            #if _DBG_H2O_TELEPORT
            cout << "Deleted water molecule " << mol << endl;
            #endif
        }
        break;
    }
}

float teleport_water(Molecule* mol)
{
    if (!waters) return -1000;

    int i, j;
    float e;
    for (j=0; waters[j]; j++)
        if (waters[j] == mol) break;
    
    if (!waters[j]) return -1000;

    for (i=0; i<_water_teleport_tries; i++)
    {
        Point teleport(
            ligcen_target.x + frand(-search_size.x, search_size.x),
            ligcen_target.y + frand(-search_size.y, search_size.y),
            ligcen_target.z + frand(-search_size.z, search_size.z)
                      );
        waters[j]->recenter(teleport);
        e = waters[j]->get_intermol_binding(gcfmols).summed();
        if (e < _water_satisfaction_threshold) break;
    }
    if (e > _water_satisfaction_threshold) delete_water(mol);
    #if _DBG_H2O_TELEPORT
    else cout << "Teleported water molecule " << mol << endl;
    #endif

    return e;
}

int interpret_resno(const char* field)
{
    char buffer[strlen(field)+4];
    strcpy(buffer, field);

    char* offset = buffer;
    while (*offset >= 'A' && *offset <= 'Z') offset++;

    int retval = 0;
    char* dot = strchr(buffer, '.');
    char* bang = strchr(buffer, '!');
    if (dot)
    {
        *(dot++) = 0;
        int b = atoi(offset);
        int w = atoi(dot);
        int _50 = protein->get_bw50(b);
        if (_50 < 1)
        {
            cerr << "Error: unknown BW number " << b << "." << w << ", please ensure PDB file has REMARK 800 SITE BW words." << endl;
            throw 0xbad12e5;
        }
        if (_50 < 1) return 0;
        else retval = _50 + w - 50;
    }
    else retval = atoi(offset);

    AminoAcid* aa = protein->get_residue(retval);
    if (!aa) return 0;

    if (offset == buffer)
    {
        if (bang)
        {
            aa->priority = true;
            priority_resnos.push_back(aa->get_residue_no());
        }
        return retval;
    }

    int i;
    for (i=0; buffer[i] >= 'A' && buffer[i] <= 'Z'; i++)
    {
        if (buffer[i] == aa->get_letter())
        {
            if (bang)
            {
                aa->priority = true;
                priority_resnos.push_back(aa->get_residue_no());
            }
            return retval;
        }
    }

    return bang ? retval : 0;
}

void freeze_bridged_residues()
{
    int i, l;

    if (bridges.size())
    {
        for (i=0; i<bridges.size(); i++)
        {
            int resno1 = interpret_resno(bridges[i].c_str());
            if (!resno1) continue;
            const char* r2 = strchr(bridges[i].c_str(), '|');
            if (!r2) throw 0xbadc0de;
            r2++;
            int resno2 = interpret_resno(r2);
            if (!resno2) continue;
            
            AminoAcid *aa1 = protein->get_residue(resno1), *aa2 = protein->get_residue(resno2);
            if (aa1)
            {
                aa1->movability = MOV_PINNED;
                aa1->been_flexed = true;
            }
            if (aa2)
            {
                aa2->movability = MOV_PINNED;
                aa2->been_flexed = true;
            }
        }
    }

    if (forced_static_resnos.size())
    {
        for (i=0; i<forced_static_resnos.size(); i++)
        {
            forced_static_resnos[i].resolve_resno(protein);
            int resno = forced_static_resnos[i].resno;
            if (!resno) continue;

            AminoAcid *aa = protein->get_residue(resno);
            if (!aa) continue;

            aa->movability = MOV_PINNED;
            aa->been_flexed = true;
        }
    }

    protein->set_clashables();
}

void reconnect_bridges()
{
    int i;
    for (i=0; i<bridges.size(); i++)
    {
        int resno1 = interpret_resno(bridges[i].c_str());
        if (!resno1) continue;
        const char* r2 = strchr(bridges[i].c_str(), '|');
        if (!r2) throw 0xbadc0de;
        r2++;
        int resno2 = interpret_resno(r2);
        if (!resno2) continue;

        #if _dbg_bridges
        cout << "Bridging " << resno1 << " and " << resno2 << "..." << endl;
        #endif

        protein->bridge(resno1, resno2);

        AminoAcid *aa1 = protein->get_residue(resno1), *aa2 = protein->get_residue(resno2);
        if (aa1) aa1->movability = MOV_PINNED;
        if (aa2) aa2->movability = MOV_PINNED; 

        #if _dbg_bridges
        if (!aa1) cout << endl << resno1 << " not found." << endl;
        if (!aa2) cout << endl << resno2 << " not found." << endl;
        if (aa1 && aa2)
        {
            float tb = aa1->get_intermol_binding(aa2).summed();
            cout << "Bridge energy " << tb << " kJ/mol." << endl;
        }
        #endif
    }

    freeze_bridged_residues();
}

void do_pivotal_hbond_rot_and_scoot()
{
    // Rotate ligand so atom faces side chain atom.
    ligand->movability = MOV_ALL;
    float r = pivotal_hbond_aaa->distance_to(pivotal_hbond_la);
    Point cen = ligand->get_barycenter();
    Rotation rot = align_points_3d(pivotal_hbond_la->loc, pivotal_hbond_aaa->loc, cen);
    LocatedVector lv = rot.v;
    lv.origin = cen;
    ligand->rotate(lv, rot.a);

    #if _dbg_priority_hbond
    cout << "Rotated ligand " << (rot.a*fiftyseven) << "deg." << endl;
    cout << "Atoms were " << r << "Å apart, now " << pivotal_hbond_aaa->distance_to(pivotal_hbond_la) << "Å" << endl;
    #endif

    // Measure the distance between atoms and move the ligand to optimize that distance.
    Vector scooch = pivotal_hbond_aaa->loc.subtract(pivotal_hbond_la->loc);
    if (scooch.r > 2.0)
    {
        scooch.r -= 2.0;
        ligand->move(scooch);

        #if _dbg_priority_hbond
        cout << "Scooched ligand " << scooch.r << "Å." << endl;
        cout << "Ligand was centered at " << cen << "; now " << ligand->get_barycenter() << endl;
        cout << "Atoms are now " << pivotal_hbond_aaa->distance_to(pivotal_hbond_la) << "Å apart" << endl;
        #endif
    }

    // Rotate the ligand about the hbond in order to minimize the intermolecular energy.
    Vector v = pivotal_hbond_la->loc.subtract(pivotal_hbond_aaa->loc);
    lv = v;
    lv.origin = pivotal_hbond_aaa->loc;
    float theta = 0, th, step = M_PI/50, clash;
    AminoAcid* reaches[SPHREACH_MAX+4];
    protein->get_residues_can_clash_ligand(reaches, ligand, ligand->get_barycenter(), Point(5,5,5), nullptr);
    for (th=0; th<circle; th += step)
    {
        float c = ligand->get_intermol_clashes((Molecule**)reaches);

        if (!th || c < clash)
        {
            clash = c;
            theta = th;
        }

        ligand->rotate(lv, step);
    }

    ligand->rotate(lv, theta);
}

void output_iter(int iter, Molecule** mols, std::string mesg)
{
    itersfname = (std::string)"tmp/" + (std::string)"_iters.dock";
    int i, liter = iter + movie_offset;
    FILE* fp = fopen(itersfname.c_str(), ((liter == 0 && pose == 1) ? "wb" : "ab") );
    if (fp)
    {
        if (!liter && (pose == 1))
        {
            fprintf(fp, "PDB file: %s\n", protfname);
        }
        fprintf(fp, "Pose: %d\nNode: %d\nMesg: %s\n\n", pose, liter, mesg.c_str());
        int foff = 0;

        #if _dock_result_in_iter
        DockResult ldr(protein, ligand, search_size, nullptr, pose, waters, true);
        ldr.include_pdb_data = false;
        ldr.display_clash_atoms = true;
        std::stringstream stst;
        stst << ldr;
        fprintf(fp, "%s\n", stst.str().c_str());
        #endif

        fprintf(fp, "\nPDBDAT:\n");

        for (i=0; reaches_spheroid[nodeno][i]; i++)
        {
            reaches_spheroid[nodeno][i]->save_pdb(fp, foff);
            foff += reaches_spheroid[nodeno][i]->get_atom_count();
        }

        for (i=0; mols[i]; i++)
        {
            if (mols[i]->is_residue()) continue;
            mols[i]->save_pdb(fp, foff, false);
            foff += mols[i]->get_atom_count();
        }

        protein->end_pdb(fp);

        fclose(fp);
    }
}

void output_iter(int iter, Molecule** mols)
{
    output_iter(iter, mols, "");
}

void abhor_vacuum(int iter, Molecule** mols)
{
    #if allow_abhor_vacuum
    int i, n = ligand->get_atom_count();
    Point rel(0,0,0);
    LocatedVector lv;

    for (i=0; i<n; i++)
    {
        Atom* a1 = ligand->get_atom(i);
        Point vail = a1->loc;
        Atom* a2 = protein->get_nearest_atom(vail);
        Vector d = a2->loc.subtract(vail);
        d.r -= a1->vdW_radius;
        d.r -= a2->vdW_radius;
        if (d.r > 0)
        {
            rel = rel.add(d);
        }
        if (ligand->stay_close_mine && (a1 != ligand->stay_close_mine))
        {
            lv = lv.add(d);
            lv.origin = lv.origin.add(vail);
        }
    }

    if (ligand->stay_close_mine)
    {
        Point stayloc = ligand->stay_close_mine->loc;
        Rotation rot = align_points_3d(lv.origin, lv.origin.add(lv), stayloc);
        if (rot.a > hexagonal/3) rot.a = hexagonal/3;
        Pose putitback(ligand);
        Interaction before = ligand->get_intermol_binding(mols);
        lv = rot.v;
        lv.origin = stayloc;
        rot.a *= frand(0,1);
        ligand->rotate(lv, rot.a);
        Interaction after = ligand->get_intermol_binding(mols);
        if (!after.accept_change(before)) putitback.restore_state(ligand);
        else if (audit) fprintf(audit, "Iter %d accepted vacuum abhorrence rotation of %f deg from %g (%g/%g) to %g (%g/%g).\n",
            iter, rot.a*fiftyseven,
            before.summed(), before.attractive, before.clash,
            after.summed(), after.attractive, after.clash
            );
    }

    Vector motion = rel;
    motion.r *= frand(0,1);
    if (motion.r > speed_limit/5) motion.r = speed_limit/5;
    Pose putitback(ligand);
    Interaction before = ligand->get_intermol_binding(mols);
    ligand->move(motion);
    Interaction after = ligand->get_intermol_binding(mols);
    if (!after.accept_change(before)) putitback.restore_state(ligand);
    else if (audit) fprintf(audit, "Iter %d accepted vacuum abhorrence translation of %f A from %g (%g/%g) to %g (%g/%g).\n",
        iter, motion.r,
        before.summed(), before.attractive, before.clash,
        after.summed(), after.attractive, after.clash
        );
    #endif
}

#if _dbg_atom_mov_to_clash
void check_moved_atom_for_clashes(Atom* caller, void* prot)
{
    int resno = caller->residue;
    if (resno)
    {
        Protein* p = reinterpret_cast<Protein*>(prot);
        AminoAcid* aa = p->get_residue(resno);
        if (aa)
        {
            if (!aa->mclashables) p->set_clashables(resno);
            float c = aa->get_intermol_clashes(aa->mclashables);
            if (c > clash_limit_per_aa)
            {
                if (!movclash_justtesting)
                {
                    throw 0xbadc0de;
                }
                movclash_justtesting = false;
            }
        }
    }
}
#endif

void iteration_callback(int iter, Molecule** mols)
{
    int i, j, l, n, ac;
    float progress, bbest;
    Point bary;

    float occl = ligand->surface_occlusion(mols);
    if (occl < frand(0.4, 0.7))
    {
        ligand->recenter(loneliest);
        Search::align_targets(ligand, loneliest, &g_bbr[0], 1);
    }

    Interaction e = ligand->get_intermol_binding(mols);
    if (e.summed() >= 50*max(1, (int)fabs(ligand->get_charge()))
        && ligand->get_worst_clash() <= clash_limit_per_atom)
    {
        end_iterations = true;
        goto _oei;
    }

    freeze_bridged_residues();

    abhor_vacuum(iter, mols);

    // Stochastically force flexion on some side chains that get clashes.
    #if stochastic_flexion_of_clashing_residues
    for (l=0; mols[l]; l++)
    {
        if (mols[l]->movability & MOV_PINNED) continue;
        Interaction li = mols[l]->get_intermol_binding(mols);
        float lf = li.summed(), lc = li.clash, ptnl = mols[l]->get_intermol_potential(mols);

        int lres = mols[l]->is_residue();
        float prob = 1;
        AminoAcid* aa = nullptr;
        if (lres)
        {
            aa = protein->get_residue(lres);
            if (aa)
            {
                prob = aa->get_aa_definition()->flexion_probability;
                AminoAcid** cr = protein->clashable_residues(aa->get_residue_no());
                if (cr)
                    prob *= 1.0 - aa->estimated_stability((Molecule**)cr);
                if (aa->coordmtl) continue;
                if (aa->is_ic_res) prob /= 2;
            }
        }

        float mc = mols[l]->lastmc;

        if (flex
            // && iter < 5
            && (lc > clash_limit_per_aa || lf > -0.3 * ptnl || (lc > 0 && lf > 5) || mc >= 1)
            && mols[l]->movability == MOV_FLXDESEL
            && lres
            && frand(0,1) < flexion_probability_multiplier * prob
            )
        {
            mols[l]->movability = MOV_FORCEFLEX;
        }
    }
    #endif

    // Attempt to connect hydrogen bonds to ligand.
    #if attempt_to_connect_hydrogen_bonds_to_ligand
    if (flex)
    {
        n = ligand->get_atom_count();
        for (l=0; mols[l]; l++)
        {
            if (mols[l]->movability & MOV_PINNED) continue;
            if (!mols[l]->is_residue()) continue;
            if (fabs(mols[l]->hydrophilicity()) < hydrophilicity_cutoff) continue;

            AminoAcid* hbaa = reinterpret_cast<AminoAcid*>(mols[l]);
            Atom* reach = hbaa->get_reach_atom();
            if (!reach) continue;
            if (fabs(reach->is_polar()) < hydrophilicity_cutoff) continue;
            Bond* rapb = reach->get_bond_by_idx(0);
            if (!rapb) continue;
            Atom* prev = rapb->atom2;
            if (!prev) continue;
            Atom* CB = hbaa->get_atom("CB");
            if (!CB) continue;
            int reachz = reach->Z;
            Atom* target = nullptr;
            float nearest = Avogadro;

            for (i=0; i<n; i++)
            {
                Atom* la = ligand->get_atom(i);
                if (!la) continue;
                if (fabs(la->is_polar()) < hydrophilicity_cutoff) continue;
                if (reachz > 1 && la->Z > 1) continue;
                float rCB = la->distance_to(CB);
                if (rCB > _DEFAULT_INTERA_R_CUTOFF + hbaa->get_reach()) continue;
                float rCA = la->loc.get_3d_distance(hbaa->get_CA_location());
                if (rCB > rCA) continue;
                float rp = la->distance_to(prev);
                if (rp > _DEFAULT_INTERA_R_CUTOFF) continue;
                float r = fmin(rp, rCB/1.25);
                if (r < nearest)
                {
                    nearest = r;
                    target = la;
                }
            }

            if (target)
            {
                if (reachz == 1 && target->Z == 1)
                {
                    Atom* oshibka = nullptr;
                    if (rand() & 1) reach = reach->get_bond_by_idx(0)->atom2;
                    else oshibka = target->get_bond_by_idx(0)->atom2;
                    if (oshibka) target = oshibka;
                }

                if (target->distance_to(reach) < 3) continue;

                Point pttgt = target->loc;
                Pose hbwas(mols[l]);
                hbwas.copy_state(mols[l]);
                Interaction before = ligand->get_intermol_binding(mols, true, true);
                hbaa->conform_atom_to_location(reach->name, pttgt, 10, 2);
                Interaction after = ligand->get_intermol_binding(mols, true, true);
                if (!after.improved(before)) hbwas.restore_state(mols[l]);
                else if (audit) fprintf(audit, "Iter %d accepted direction of %s:%s toward LIG:%s from %g (%g/%g) to %g (%g/%g).\n",
                    iter, mols[l]->get_name(), reach->name, target->name,
                    before.summed(), before.attractive, before.clash, after.summed(), after.attractive, after.clash);
            }
        }
    }
    #endif

    if (pivotal_hbond_aaa && pivotal_hbond_la) do_pivotal_hbond_rot_and_scoot();

    bary = ligand->get_barycenter();

    ac = ligand->get_atom_count();
    bbest = 0;
    Atom *atom1, *atom2;

    progress = (float)iter / iters;

    // Soft docking.
    if (n && iter>soft_iter_min) soft_docking_iteration(protein, ligand, nsoftrgn, softrgns, softness);

    #if bb_realign_iters
    for (l=0; g_bbr[l].pri_res && g_bbr[l].pri_tgt; l++)
    {
        Molecule* llig = ligand->get_monomer(l);
        Interaction before = llig->get_intermol_binding(mols);
        Pose was(llig);
        #if bb_realign_only_hydro
        Search::align_targets(llig, loneliest, &g_bbr[l], bb_realign_amount);
        #else
        Search::align_targets(llig, loneliest, &g_bbr[l], bb_realign_amount);
        #endif
        Interaction after = llig->get_intermol_binding(mols);
        if (!after.improved(before)) was.restore_state(llig);
    }
    #endif

    if (!iter) goto _oei;
    if (iter == (iters-1)) goto _oei;

    // BB enforcement.
    #if enforce_no_bb_pullaway
    if (pdpst == pst_best_binding)
    {
        // TODO: This functionality was never written. Has it been obviated by the more recent "stays" feature?
    }
    #endif

    bary = ligand->get_barycenter();

    for (i=0; i < ac; i++)
    {
        if (ligand->get_atom(i)->strongest_bind_energy > bbest)
        {
            atom1 = ligand->get_atom(i);
            bbest = atom1->strongest_bind_energy;
            atom2 = atom1->strongest_bind_atom;
        }
    }

    if (bbest >= 15)
    {
        bary = ligand->get_barycenter();
    }
    else
    {
        #if allow_drift

        Pose predrift(ligand);
        Interaction predrift_binding = ligand->get_intermol_binding(mols, true, true);

        if (predrift_binding.summed() > drift_energy_threshold)
        {
            #if !pocketcen_is_loneliest
            if (ligand->lastbind >= 100)
            {
                ligcen_target.x += (loneliest.x - ligcen_target.x) * drift;
                ligcen_target.y += (loneliest.y - ligcen_target.y) * drift;
                ligcen_target.z += (loneliest.z - ligcen_target.z) * drift;
            }
            #endif

            Point was = ligand->get_barycenter();
            if (bary.get_3d_distance(ligcen_target) > size.magnitude())
            {
                //cout << "Wrangle! " << bary << ": " << bary.get_3d_distance(ligcen_target) << " vs. " << size.magnitude() << endl;
                bary = ligcen_target;
                // ligand->reset_conformer_momenta();
            }
            else
            {
                if (ligand->lastbind > 0)
                {
                    bary.x += (ligcen_target.x - bary.x) * drift;
                    bary.y += (ligcen_target.y - bary.y) * drift;
                    bary.z += (ligcen_target.z - bary.z) * drift;
                }
                else drift *= (1.0 - drift_decay_rate/iters);
            }

            Interaction before = ligand->get_intermol_binding(mols);
            ligand->recenter(bary);
            Interaction after = ligand->get_intermol_binding(mols);
            if (!after.accept_change(before)) ligand->recenter(was);
            else if (audit) fprintf(audit, "Iter %d applied ligand drift of %f from %f to %f.\n", iter, drift, before.summed(), after.summed());
        }

        #endif
    }

    #if _teleport_dissatisfied_waters
    if (waters && (iter % 5) == 4)
    {
        for (i=0; waters[i]; i++)
        {
            float r = waters[i]->get_barycenter().get_3d_distance(ligcen_target);
            if (r > size.magnitude()) teleport_water(waters[i]);
            else
            {
                float e = 0;
                int j;
                for (j=0; j<10; j++)
                {
                    if (!j || waters[i]->lastbind_history[j] < e) e = waters[i]->lastbind_history[j];
                }
                if (e > _water_satisfaction_threshold) teleport_water(waters[i]);
            }
        }
    }
    #endif

    if (waters)
    {
        for (i=0; waters[i]; i++)
        {
            if (waters[i]->movability & MOV_PINNED) continue;
            float rpl = ligand->get_intermol_binding(waters[i]).clash;
            if (rpl > 50)
            {
                for (j=i+1; waters[j]; j++) waters[j-1] = waters[j];
                waters[j-1] = nullptr;
                continue;
            }

            Atom *a, *b;
            ligand->mutual_closest_atoms(waters[i], &a, &b);
            if (!a || !b) continue;

            float optimal = InteratomicForce::optimal_distance(a, b);
            if (optimal < 1) optimal = a->vdW_radius + b->vdW_radius + 1;
            float r = a->distance_to(b);
            if (r <= optimal) continue;

            Vector movamt = a->loc.subtract(b->loc);
            movamt.r = fmin(speed_limit, (r-optimal)/5);
            waters[i]->move(movamt);

            rpl = ligand->get_intermol_binding(waters[i]).clash;
            if (rpl > clash_limit_per_atom)
            {
                movamt.r *= -1;
                waters[i]->move(movamt);
            }
        }
    }

    _oei:
    ;

    #if recapture_ejected_ligand
    Point lig_center = ligand->get_barycenter();
    float r = lig_center.get_3d_distance(ligcen_target);
    float recapture_distance = size.magnitude() / 2;
    if (r >= recapture_distance) ligand->recenter(ligcen_target);
    #endif

    for (i=0; i<nappears; i++)
    {
        if (iter == appears[i].disappear_iter || iter == appears[i].reappear_iter)
        {
            appears[i].dispcen = protein->get_region_center(1, protein->get_end_resno()); // loneliest;
            appears[i].update(protein, iter);
        }
    }

    if (output_each_iter) output_iter(iter, mols, "dock iteration");
}

void do_pose_output(DockResult* drjk, int lnodeno, float energy_mult, Pose* tmp_pdb_water, Point* tmp_pdb_metal_loc)
{

    cout << "Pose: " << pose << endl << "Node: " << lnodeno << endl;
    if (output) *output << "Pose: " << pose << endl << "Node: " << lnodeno << endl;

    drjk->energy_mult = energy_mult;
    drjk->do_output_colors = do_output_colors;
    drjk->include_pdb_data = (output == nullptr);
    cout << *drjk;
    drjk->do_output_colors = false;
    drjk->include_pdb_data = true;
    if (output) *output << *drjk;

    if (!lnodeno && outpdb.length() && pose <= outpdb_poses)
    {
        char protn[64];
        strcpy(protn, strrchr(protfname, '/')+1);
        char* dot = strchr(protn, '.');
        if (dot) *dot = 0;

        char lign[64];
        strcpy(lign, strrchr(ligfname, '/')+1);
        dot = strchr(lign, '.');
        if (dot) *dot = 0;

        // std::string temp_pdb_fn = (std::string)"tmp/pose" + std::to_string(j+1) + (std::string)".pdb";
        std::string out_pdb_fn = std::regex_replace(outpdb, std::regex("[%][p]"), protn);
        out_pdb_fn = std::regex_replace(out_pdb_fn, std::regex("[%][l]"), lign);
        out_pdb_fn = std::regex_replace(out_pdb_fn, std::regex("[%][o]"), to_string(pose));

        // FILE* pftmp = fopen(temp_pdb_fn.c_str(), "rb");
        FILE* pfout = fopen(out_pdb_fn.c_str(), "wb");
        if (/*!pftmp ||*/ !pfout)
        {
            cerr << "Failed to open output PDB file " << out_pdb_fn << " for writing." << endl;
            throw 0xbadf12e;
        }

        int j1;
        int n1 = protein->get_end_resno();
        for (j1=1; j1 <= n1; j1++)
        {
            AminoAcid* aa = protein->get_residue(j1);
            if (aa) aa->set_pdb_chain('A');
            if (aa /* && aa->been_flexed */)
            {
                // tmp_pdb_residue[j+1][j1].restore_state_relative(aa, "CA");
                #if _dbg_residue_poses
                cout << "tmp_pdb_residue[" << (j+1) << "][" << j1 << "].restore_state_relative(" 
                    << aa->get_name() << ")" << endl;
                #endif
            }
        }
        // tmp_pdb_ligand[j+1].restore_state(ligand);

        #if recapture_ejected_ligand
        float r = ligand->get_barycenter().get_3d_distance(nodecens[k]);
        if (r > size.magnitude()/2) goto _next_pose;
        #endif

        protein->save_pdb(pfout, ligand);

        int atno_offset = protein->last_saved_atom_number;
        if (waters)
        {
            for (j1=0; waters[j1]; j1++)
            {
                tmp_pdb_water[j1].restore_state(waters[j1]);
                waters[j1]->save_pdb(pfout, atno_offset);
                atno_offset += waters[j1]->get_atom_count();
            }
        }

        n1 = nmtlcoords;
        for (j1=0; j1 < n1; j1++)
        {
            mtlcoords[j1].mtl->move(tmp_pdb_metal_loc[j1]);
            mtlcoords[j1].mtl->save_pdb_line(pfout, ++atno_offset);
        }

        protein->end_pdb(pfout);

        fclose(pfout);
    }
}

void search_callback(std::string mesg)
{
    if (output_each_iter)
    {
        output_iter(0, gcfmols, mesg);
        movie_offset++;
    }
}

void update_progressbar(float percentage)
{
    int n = max(1, pathnodes);
    percentage = percentage/poses/n + (float)(pose-1+(float)nodeno/n)*100.0/poses;
    if (success_sofar)
    {
        if (strstr(protfname, ".inactive.pdb"))
            progb.set_color(98, 176, 224);
        else
            progb.set_color(96, 216, 128);
    }
    else progb.set_color(160, 168, 176);
    progb.update(percentage);
}

Point pocketcen_from_config_words(char** words, Point* old_pocketcen)
{
    int i=1;
    Point local_pocketcen;
    center_resnos.clear();
    priority_resnos.clear();
    if (!strcmp(words[i], "RES"))
    {
        i++;
        for (; words[i]; i++)
        {
            int j = interpret_resno(words[i]);
            if (!j) continue;

            center_resnos.push_back(j);
        }

        int sz = center_resnos.size(), div=0;
        Point foravg[sz + 2];
        for (i=0; i<sz; i++)
        {
            #if pocketcen_from_reach_atoms
            AminoAcid* aa = protein->get_residue(center_resnos[i]);
            if (aa)
            {
                foravg[i] = aa->get_reach_atom_location();
                div++;
            }
            #else
            foravg[div++] = protein->get_atom_location(center_resnos[i], "CA");
            #endif
        }

        return average_of_points(foravg, div?:1);
    }
    else if (!strcmp(words[i], "REL"))
    {
        if (!old_pocketcen)
        {
            cerr << "Error: relative coordinates not supported for CEN." << endl << flush;
            throw 0xbadb19d;
        }
        else
        {
            i++;
            local_pocketcen.x = old_pocketcen->x + atof(words[i++]);
            local_pocketcen.y = old_pocketcen->y + atof(words[i++]);
            local_pocketcen.z = old_pocketcen->z + atof(words[i++]);
            return local_pocketcen;
        }
    }
    else
    {
        if (!strcmp(words[i], "ABS")) i++;
        local_pocketcen.x = atof(words[i++]);
        local_pocketcen.y = atof(words[i++]);
        local_pocketcen.z = atof(words[i++]);
        return local_pocketcen;
    }
}

int interpret_config_line(char** words)
{
    int i;

    optsecho = "";

    if (0) { ; }
    else if (!strcmp(words[0], "APPENDPROT") || !strcmp(words[0], "OPEND"))
    {
        append_pdb = true;
    }
    else if (!strcmp(words[0], "ATOMTO"))
    {
        atomto.push_back(origbuff);
    }
    else if (!strcmp(words[0], "NORESWARN"))
    {
        warn_absent_residue = false;
    }
    else if (!strcmp(words[0], "AUDIT"))
    {
        auditfn = "tmp/dock_audit";
    }
    else if (!strcmp(words[0], "BRIDGE"))
    {
        std::string str = words[1];
        str += (std::string)"|" + (std::string)words[2];
        bridges.push_back(str);
        return 2;
    }
    else if (!strcmp(words[0], "CEN"))
    {
        CEN_buf.push_back(origbuff);
        optsecho = (std::string)"Center " + (std::string)origbuff;
        return 0;
    }
    else if (!strcmp(words[0], "CNTCT"))
    {
        if (file_exists(words[1]))
        {
            hg.load_ic_file(words[1]);
        }
        else
        {
            soft_contact_a[n_soft_contact].set(words[1]);
            soft_contact_b[n_soft_contact].set(words[2]);
            n_soft_contact++;
        }
        return 2;
    }
    else if (!strcmp(words[0], "COLORS"))
    {
        do_output_colors = true;
        return 0;
    }
    else if (!strcmp(words[0], "COLORLESS"))
    {
        do_output_colors = false;
        return 0;
    }
    else if (!strcmp(words[0], "DEACVNODE"))
    {
        cout << "Notice: DEACVNODE has been deprecated in favor of the HXR option. Please update your config files." << endl;
        return 0;
    }
    else if (!strcmp(words[0], "DEBUG"))
    {
        if (!words[1])
        {
            cout << "Missing debug file name; check config file." << endl << flush;
            throw 0xbadf12e;
        }
        #if _DBG_STEPBYSTEP
        cout << "Starting a debug outstream." << endl;
        #endif
        debug = new std::ofstream(words[1], std::ofstream::out);
        optsecho = "Debug file: " + (std::string)words[1];
        return 1;
    }
    else if (!strcmp(words[0], "DISAPPEAR"))
    {
        int diter = atoi(words[1]);
        appears[nappears].disappear_iter = diter;
        if (diter > last_appear_iter) last_appear_iter = diter;
        appears[nappears].start.set(words[2]);
        appears[nappears].end.set(words[3]);
        nappears++;
    }
    else if (!strcmp(words[0], "REAPPEAR"))
    {
        int riter = atoi(words[1]);
        appears[nappears].reappear_iter = riter;
        if (riter > last_appear_iter) last_appear_iter = riter;
        appears[nappears].start.set(words[2]);
        appears[nappears].end.set(words[3]);
        if (words[4]) appears[nappears].displacement = atof(words[4]);
        nappears++;
    }
    else if (!strcmp(words[0], "ECHO"))
    {
        echo_progress = true;
        optsecho = "Echo on.";
    }
    else if (!strcmp(words[0], "ELIM"))
    {
        kJmol_cutoff = atof(words[1]);
        optsecho = "Energy limit: " + to_string(kJmol_cutoff);
        return 1;
    }
    else if (!strcmp(words[0], "EXCL"))
    {
        i=1;
        int excls = atoi(words[i++]);
        int excle = atoi(words[i++]);

        for (i=excls; i<=excle; i++) exclusion.push_back(i);
        optsecho = "Exclude range " + to_string(excls) + (std::string)"-" + to_string(excle);
        return i-1;
    }
    else if (!strcmp(words[0], "FLEX"))
    {
        flex = (atoi(words[1]) != 0);
        optsecho = "Flex: " + (std::string)(flex ? "ON" : "OFF");
        return 1;
    }
    else if (!strcmp(words[0], "FLXR"))
    {
        i = 1;
        while (words[i])
        {
            ResiduePlaceholder rph;
            rph.set(words[i]);
            forced_flexible_resnos.push_back(rph);
            i++;
        }
        return i-1;
    }
    else if (!strcmp(words[0], "STCR"))
    {
        i = 1;
        while (words[i])
        {
            ResiduePlaceholder rph;
            rph.set(words[i]);
            forced_static_resnos.push_back(rph);
            i++;
        }
        return i-1;
    }
    else if (!strcmp(words[0], "H2O"))
    {
        maxh2o = omaxh2o = atoi(words[1]);
        optsecho = "Water molecules: " + to_string(maxh2o);
        return 1;
    }
    else if (!strcmp(words[0], "WETBB"))
    {
        maxh2o++;
        omaxh2o = maxh2o;
        wet_contact = true;
        return 0;
    }
    else if (!strcmp(words[0], "DRYBB"))
    {
        wet_contact = false;
        return 0;
    }
    else if (!strcmp(words[0], "HYDRO"))
    {
        hydrogenate_pdb = true;
    }
    else if (!strcmp(words[0], "ITER") || !strcmp(words[0], "ITERS"))
    {
        iters = atoi(words[1]);
        optsecho = "Iterations: " + to_string(iters);
        return 1;
    }
    else if (!strcmp(words[0], "KCAL"))
    {
        kcal = true;
        optsecho = "Output units switched to kcal/mol.";
        return 0;
    }
    else if (!strcmp(words[0], "LIG"))
    {
        strcpy(ligfname, words[1]);
        ligset = true;
        if (ligcmd) smset = smcmd;
        return 1;
    }
    else if (!strcmp(words[0], "COPYLIG"))
    {
        if (!strcmp(words[1], "E")) estimate_copylig = true;
        else copylig = atoi(words[1]);
        return 1;
    }
    else if (!strcmp(words[0], "ISO"))
    {
        isomers.push_back(words[1]);
        return 1;
    }
    else if (!strcmp(words[0], "FORM"))
    {
        // Even though no constraint requires a form to be an isomer,
        // e.g. hydrolyzed forms of lactones are allowed, in the code we
        // treat them identically to isomers since there's no difference in processing.
        isomers.push_back(words[1]);
        return 1;
    }
    else if (!strcmp(words[0], "SMILES"))
    {
        strcpy(smiles, words[1]);
        smset = true;
        return 1;
    }
    else if (!strcmp(words[0], "MCOORD"))
    {
        int j=0;
        // optsecho = "Metal coordination on residues ";
        MCoord mcr;
        i=1; if (!words[i]) throw 0xbad372;
        mcr.Z = Atom::Z_from_esym(words[1]);
        if (!mcr.Z) throw 0xbad372;

        i++; if (!words[i]) throw 0xbad372;
        mcr.charge = atoi(words[i]);

        i++; if (!words[i]) throw 0xbad372;
        for (; words[i]; i++)
        {
            if (words[i][0] == '-' && words[i][1] == '-') break;
            if (words[i][0] == '#') break;
            ResiduePlaceholder rp;
            rp.set(words[i]);
            mcr.coordres[mcr.ncoordres++] = rp;
        }

        mtlcoords[nmtlcoords++] = mcr;
        return i-1;
    }
    else if (!strcmp(words[0], "MOVIE"))
    {
        output_each_iter = true;
    }
    else if (!strcmp(words[0], "NODEPDB"))
    {
        activation_node = atoi(words[1]);
        strcpy(protafname, words[2]);
        optsecho = "Active PDB " + (std::string)protafname + " for node " + to_string(activation_node);
    }
    else if (!strcmp(words[0], "NOFAIL"))
    {
        output_something_even_if_it_is_wrong = true;
    }
    else if (!strcmp(words[0], "OUT"))
    {
        if (!words[1])
        {
            cout << "Missing output file name; check config file." << endl << flush;
            throw 0xbadf12e;
        }
        strcpy(outfname, words[1]);
        optsecho = "Output file is " + (std::string)outfname;
        return 1;
    }
    else if (!strcmp(words[0], "OUTPDB"))
    {
        outpdb_poses = atoi(words[1]);
        outpdb = words[2];
        return 2;
    }
    else if (!strcmp(words[0], "PATH"))
    {
        i=1;
        int nodeno = atoi(words[i]);

        if (!strcmp(words[i], "ABS")) i++;
        if (!strcmp(words[i], "REL")) i++;
        if (!strcmp(words[i], "RES")) i++;

        if (nodeno > 255)
        {
            cout << "Binding path is limited to 255 nodes." << endl << flush;
            throw 0xbad90de;
        }
        if (nodeno)
        {
            if ((nodeno) > pathnodes) pathnodes = nodeno;

            pathstrs.resize(nodeno+1);
            pathstrs[nodeno] = origbuff;
        }
        optsecho = "Path node set #" + to_string(nodeno);
        return i-1;
    }
    else if (!strcmp(words[0], "PERRES"))
    {
        out_per_res_e = atoi(words[1]);
        return 1;
    }
    else if (!strcmp(words[0], "PERBTYP"))
    {
        out_per_btyp_e = atoi(words[1]);
        return 1;
    }
    else if (!strcmp(words[0], "ELIMITEM"))
    {
        out_itemized_e_cutoff = atoi(words[1]);
        return 1;
    }
    else if (!strcmp(words[0], "LIGINTE"))
    {
        out_lig_int_e = words[1] ? atoi(words[1]) : true;
        return 1;
    }
    else if (!strcmp(words[0], "OUTBBP"))
    {
        out_bb_pairs = words[1] ? atoi(words[1]) : true;
        return 1;
    }
    else if (!strcmp(words[0], "OUTDISQ"))
    {
        outdisq = true;
        return 0;
    }
    else if (!strcmp(words[0], "OUTLPS"))
    {
        out_lig_pol_sat = atoi(words[1]);
        return 1;
    }
    else if (!strcmp(words[0], "OUTPROX"))
    {
        out_prox = atoi(words[1]);
        return 1;
    }
    else if (!strcmp(words[0], "OUTCLSHA"))
    {
        display_clash_atoms = atoi(words[1]);
        return 1;
    }
    else if (!strcmp(words[0], "OUTPCLSH"))
    {
        out_pro_clash = atoi(words[1]);
        return 1;
    }
    else if (!strcmp(words[0], "OUTMC"))
    {
        out_mc = atoi(words[1]);
        return 1;
    }
    else if (!strcmp(words[0], "OUTVDWR"))
    {
        out_vdw_repuls = atoi(words[1]);
        return 1;
    }
    else if (!strcmp(words[0], "OUTPDBL"))
    {
        out_pdbdat_lig = atoi(words[1]);
        return 1;
    }
    else if (!strcmp(words[0], "OUTPDBR"))
    {
        out_pdbdat_res = atoi(words[1]);
        return 1;
    }
    else if (!strcmp(words[0], "OUTRSHP"))
    {
        rshp_verbose = true;
        return 0;
    }
    else if (!strcmp(words[0], "POSE"))
    {
        poses = atoi(words[1]);
        optsecho = "Number of poses: " + to_string(poses);
        return 1;
    }
    else if (!strcmp(words[0], "PROGRESS"))
    {
        progressbar = true;
        return 0;
    }
    else if (!strcmp(words[0], "CONGRESS"))
    {
        progressbar = false;
        return 0;
    }
    else if (!strcmp(words[0], "PROT"))
    {
        strcpy(protfname, words[1]);
        char* c = strchr(protfname, ':');
        if (c)
        {
            *c = 0;
            protstrand = *(++c);
        }
        else protstrand = 0;
        protset = true;
        // optsecho = "Protein file is " + (std::string)protfname;
        return 1;
    }
    else if (!strcmp(words[0], "REQSR"))
    {
        int reqrnode = atoi(words[1]);
        if (!strcmp(words[1], "all") || !strcmp(words[1], "ALL")) reqrnode = -1;
        for (i=2; words[i]; i++)
        {
            if (words[i][0] == '-' && words[i][1] == '-') break;
            ResiduePlaceholder rp;
            rp.set(words[i]);
            rp.node = reqrnode;
            required_contacts.push_back(rp);
        }
        return i-1;
    }
    else if (!strcmp(words[0], "RETRY"))
    {
        // triesleft = atoi(words[1]);
        return 1;
    }
    else if (!strcmp(words[0], "RLIM"))
    {
        _INTERA_R_CUTOFF = atof(words[1]);
        optsecho = "Interatomic radius limit: " + to_string(_INTERA_R_CUTOFF);
        return 1;
    }
    else if (!strcmp(words[0], "SEARCH"))
    {
        if (!words[1]) return 0;       // Stay with default.
        if (!strcmp(words[1], "BB"))
        {
            pdpst = pst_best_binding;
            return 1;
        }
        else if (!strcmp(words[1], "TS"))
        {
            pdpst = pst_tumble_spheres;
            return 1;
        }
        else if (!strcmp(words[1], "CP"))
        {
            int lf = 1;
            pdpst = pst_copyfrom;
            if (!words[2])
            {
                cerr << "ERROR: Search mode CP without source file." << endl;
                throw 0xbad19b07;
            }
            copyfrom_filename = words[2];
            if (words[3])
            {
                if (strlen(words[3]) > 3) words[3][3] = 0;
                strcpy(copyfrom_ligname, words[3]);
                lf++;
                if (words[4])
                {
                    copyfrom_resno = atoi(words[4]);
                    lf++;
                }
            }
            return lf;
        }
        else if (!strcmp(words[1], "EX"))
        {
            pdpst = pst_external;
            // Get filename
            if (!words[2])
            {
                cerr << "ERROR: Search mode EX without source file." << endl;
                throw 0xbad19b07;
            }
            copyfrom_filename = words[2];
            return 2;
        }
        else
        {
            cout << "Unknown search method " << words[1] << endl;
            throw 0xbad5eec;
        }
    }
    else if (!strcmp(words[0], "SIZE"))
    {
        search_size.x = atof(words[1]);
        if (words[2])
        {
            search_size.y = atof(words[2]);
            search_size.z = atof(words[3]);
        }
        else search_size.z = search_size.y = search_size.x;
        if (!search_size.x || !search_size.y || !search_size.z)
        {
            cout << "Pocket size cannot be zero in any dimension." << endl << flush;
            throw 0xbad512e;
        }
        optsecho = "Interatomic radius limit: " + to_string(search_size.x) + (std::string)"," + to_string(search_size.y) + (std::string)"," + to_string(search_size.z);
        return 3;
    }
    else if (!strcmp(words[0], "SOFT"))
    {
        softdock = true;
        softness = 1;

        for (i=1; words[i]; i++)
        {
            if (strchr(words[i], '.'))
            {
                softness = atof(words[i]);
                #if _dbg_soft
                cout << "Softness: " << softness << endl;
                #endif
            }
            else
            {
                int j = atoi(words[i]);
                if (j) softrgn_allowed[nsoftrgna++] = j;
                #if _dbg_soft
                cout << "Allowed soft region: " << j << endl;
                #endif
            }
        }
    }
    else if (!strcmp(words[0], "NODEL"))
    {
        ResiduePlaceholder rpsr, rper;
        rpsr.set(words[1]);
        rper.set(words[2]);
        soft_nodel_start[nsoftnodel] = rpsr;
        soft_nodel_end[nsoftnodel] = rper;
        nsoftnodel++;
        return 2;
    }
    else if (!strcmp(words[0], "STATE"))
    {
        states.push_back(origbuff);
        optsecho = "Added state " + (std::string)origbuff;
        return 0;
    }
    else if (!strcmp(words[0], "TEMPLATE"))
    {
        tplset = (strcmp(words[1], "off") != 0) && (strcmp(words[1], "OFF") != 0);
        if (tplset)
        {
            strcpy(tplfname, words[1]);
            if (words[2])
            {
                tprfset = true;
                strcpy(tplrfnam, words[2]);
            }
        }
        return 1;
    }
    else if (!strcmp(words[0], "CAVS"))
    {
        cavity_stuffing = atof(words[1]);
        return 1;
    }
    else if (!strcmp(words[0], "CFLEE"))
    {
        clash_fleeing = atof(words[1]);
        return 1;
    }
    else if (!strcmp(words[0], "VCVTY"))
    {
        strcpy(cvtyfname, words[1]);
        return 1;
    }

    return 0;
}

void read_config_file(FILE* pf)
{
    char buffer[65536];
    int i;

    while (!feof(pf))
    {
        buffer[0] = 0;
        char* wgas = fgets(buffer, 1015, pf);
        char* hash = strchr(buffer, '#');
        if (hash) *hash = 0;
        origbuff = buffer;
        if (buffer[0] >= ' ' && buffer[0] != '#')
        {
            char** words = chop_spaced_words(buffer);
            if (!words) continue;
            if (!words[0]) continue;

            interpret_config_line(words);

            delete[] words;
        }
        buffer[0] = 0;
    }
}

void attempt_priority_hbond()
{
    int i, j, m, n;

    n = priority_resnos.size();
    for (i=0; i<n; i++)
    {
        AminoAcid* aa = protein->get_residue(priority_resnos[i]);
        if (fabs(aa->hydrophilicity()) >= hydrophilicity_cutoff)
        {
            // Find atom of side chain capable of hydrogen bond.
            Atom* aaa = aa->capable_of_inter(hbond);
            if (!aaa) continue;

            #if _dbg_priority_hbond
            cout << aa->get_name() << ":" << aaa->name << " found." << endl;
            #endif

            // Find nearest atom of ligand capable of hbond.
            m = ligand->get_atom_count();
            float r = Avogadro;
            Atom* la = nullptr;
            for (j=0; j<m; j++)
            {
                if (frand(0, 1) < 0.29) continue;
                Atom* a1 = ligand->get_atom(j);
                if (fabs(a1->is_polar()) >= hydrophilicity_cutoff)
                {
                    float r1 = a1->distance_to(aaa);
                    if (r1 < r)
                    {
                        r = r1;
                        la = a1;
                    }
                }
            }
            if (!la) continue;

            #if _dbg_priority_hbond
            cout << "ligand:" << la->name << " found." << endl;
            #endif

            // Ensure one atom is a donor and the other an acceptor. If one donor and one acceptor cannot be found, skip.
            if (sgn(aaa->is_polar()) == sgn(la->is_polar()))
            {
                if (aaa->Z == 1)
                {
                    Atom* heavy = aaa->get_bond_by_idx(0)->atom1;
                    if (heavy->get_family() == CHALCOGEN) aaa = heavy;
                }
                else
                {
                    Atom* hyd = aaa->is_bonded_to("H");
                    if (hyd) aaa = hyd;
                }
            }
            if (sgn(aaa->is_polar()) == sgn(la->is_polar()))
            {
                if (la->Z == 1)
                {
                    Atom* heavy = la->get_bond_by_idx(0)->atom1;
                    if (heavy->get_family() == CHALCOGEN) la = heavy;
                    else if (heavy->get_family() == PNICTOGEN && !heavy->get_charge()) la = heavy;
                }
                else
                {
                    Atom* hyd = la->is_bonded_to("H");
                    if (hyd) la = hyd;
                }
            }
            if (sgn(aaa->is_polar()) == sgn(la->is_polar())) continue;

            #if _dbg_priority_hbond
            cout << aa->get_name() << ":" << aaa->name << " and ligand:" << la->name << endl;
            #endif

            // Pin the ligand so it can only rotate about the hbond atom.
            pivotal_hbond_aaa = aaa;
            pivotal_hbond_la = la;
            pivotal_hbond_r = 2.0;

            do_pivotal_hbond_rot_and_scoot();

            #if _dbg_priority_hbond
            cout << endl;
            #endif
            return;
        }
    }
}

void choose_cen_buf()
{
    int i, n;

    n = protein->get_end_resno();
    for (i=1; i<=n; i++)
    {
        AminoAcid* aa = protein->get_residue(i);
        if (aa) aa->priority = false;
    }

    priority_resnos.clear();

    n = CEN_buf.size();
    if (pose <= n) cenbuf_idx = pose-1;
    else
    {
        for (i=0; i<10; i++)
        {
            cenbuf_idx = rand() % n;
            if (strchr(CEN_buf[cenbuf_idx].c_str(), '!')) return;
            else if (frand(0,1) < 0.1) return;
        }
    }
}

Point loneliest_from_cavities(Point inpt)
{
    Point outpt = inpt;
    if (ncvtys)
    {
        int i, j, npartial, nresno, nbsr, mostbsr=0;
        AminoAcid* cavaa[protein->get_end_resno()+16];
        for (i=0; i<ncvtys; i++)
        {
            nbsr = 0;
            nresno = cvtys[i].resnos(protein, cavaa);
            for (j=0; j<nresno; j++)
            {
                auto it = std::find(center_resnos.begin(), center_resnos.end(), cavaa[j]->get_residue_no());
                if (it != center_resnos.end()) nbsr += cavaa[j]->priority ? 5 : 1;
            }
            if (nbsr > mostbsr)
            {
                mostbsr = nbsr;
                outpt = cvtys[i].get_center();
            }
        }
    }

    return outpt;
}

void apply_protein_specific_settings(Protein* p)
{
    int i, j, n;
    choose_cen_buf();

    char buffer[1024];
    strcpy(buffer, CEN_buf[cenbuf_idx].c_str());
    char** words = chop_spaced_words(buffer);
    pocketcen = pocketcen_from_config_words(words, nullptr);
    loneliest = protein->find_loneliest_point(pocketcen, search_size);
    loneliest = loneliest_from_cavities(loneliest);
    delete[] words;

    if (n = priority_resnos.size()) for (i=0; i<n; i++)
    {
        AminoAcid* aa = protein->get_residue(priority_resnos[i]);
        if (aa) aa->priority = true;
    }

    if (!p->get_metals_count() && metald_prot) p->copy_mcoords(metald_prot);

    j=0;

    n = atomto.size();
    for (i=0; i<n; i++)
    {
        char buffer[1024];
        strcpy(buffer, atomto[i].c_str());
        char** words = chop_spaced_words(buffer);

        if (!words[1]) throw -1;
        AminoAcid* aa = protein->get_residue_bw(words[1]);
        if (!words[2]) throw -1;
        char* aname = words[2];
        if (!words[3]) throw -1;
        AminoAcid* target = protein->get_residue_bw(words[3]);
        if (words[4]) throw -1;

        if (!aa)
        {
            if (warn_absent_residue && !strstr(nfwarned.c_str(), words[1]))
            {
                cout << "Residue " << words[1] << " not present in protein; this ATOMTO command will have no effect." << endl;
                nfwarned += (std::string)" " + (std::string)words[1];
            }
            continue;
        }

        if (!target)
        {
            if (warn_absent_residue && !strstr(nfwarned.c_str(), words[3]))
            {
                cout << "Residue " << words[3] << " not present in protein; this ATOMTO command will have no effect." << endl;
                nfwarned += (std::string)" " + (std::string)words[3];
            }
            continue;
        }

        Atom* a = aa->get_atom(aname);
        if (!strcmp("EXTENT", aname)) a = aa->get_reach_atom();
        if (!a)
        {
            if (!strstr(nfwarned.c_str(), aname))
            {
                cout << "Warning: atom not found " << *aa << ":" << aname << endl;
                nfwarned += (std::string)" " + (std::string)aname;
                if (progressbar) progb.erase();
            }

            continue;
        }

        MovabilityType aamov = aa->movability;
        aa->movability = MOV_FLEXONLY;
        if (!aa->mclashables) protein->set_clashables(aa->get_residue_no());
        aa->conform_atom_to_location(a->name, target->get_CA_location());
        aa->movability = aamov;

        delete[] words;
    }

    if (softdock)
    {
        for (i=0; i<nsoftnodel; i++)
        {
            if (!soft_nodel_start[i].resno) soft_nodel_start[i].resolve_resno(protein);
            if (!soft_nodel_end[i].resno) soft_nodel_end[i].resolve_resno(protein);
        }

        n = protein->get_end_resno();
        for (i=1; i<=n; i++)
        {
            bool helixed = false;
            AminoAcid* aa = protein->get_residue(i);
            if (!aa) continue;
            if (aa->priority) continue;

            // TODO: This should not be hard coded to 6.48-6.59, and it should not depend on the protein having BW numbers.
            BallesterosWeinstein bw = protein->get_bw_from_resno(aa->is_residue());
            if (bw.helix_no == 6 && bw.member_no >= 48 && bw.member_no <= 59) continue;

            std::string rgname = (bw.helix_no < 8 ? "TMR" : (bw.helix_no == 8 ? "HXR" : 
                ((bw.helix_no == 23 || bw.helix_no == 45 || bw.helix_no == 67) ? "EXR" : "CYT")));
            if (bw.helix_no < 10) rgname += std::to_string(bw.helix_no);
            else rgname += std::to_string((int)(bw.helix_no+10)/20);
            if (i >= protein->get_region_start(rgname.c_str()) && i <= protein->get_region_end(rgname.c_str())) continue;

            bool snfound = false;
            for (j=0; j<nsoftnodel; j++)
            {
                if (i >= soft_nodel_start[j].resno && i <= soft_nodel_end[j].resno)
                {
                    snfound = true;
                    break;
                }
            }
            if (snfound) continue;

            if (aa->is_alpha_helix()) helixed = true;
            else for (j=i-2; j<=i+2; j++)
            {
                if (j < 1) continue;
                if (j == i) continue;
                if (j > n) break;
                aa = protein->get_residue(j);
                if (!aa) continue;
                if (aa->is_alpha_helix() || aa->priority) helixed = true;
                if (helixed) break;
            }

            if (!helixed) protein->delete_residue(i);
        }

        if (!nsoftrgn) for (i=1; i<=n; i++)
        {
            if (protein->get_residue(i))
            {
                Region r;
                r.start = i;
                BallesterosWeinstein bwi = protein->get_bw_from_resno(i);
                // cout << "Residue " << i << " has bw " << bwi.helix_no << "." << bwi.member_no << endl;
                for (j=i+1; j<=n; j++)
                {
                    BallesterosWeinstein bwj = protein->get_bw_from_resno(j);
                    if (j==n || !protein->get_residue(j) || bwi.helix_no != bwj.helix_no)
                    {
                        r.end = j-1;
                        if (r.end > r.start+2)
                        {
                            bool allowed = true;
                            int l = (i+j-1)/2;
                            BallesterosWeinstein bw = protein->get_bw_from_resno(l);

                            if (nsoftrgna)
                            {
                                allowed = false;
                                for (l=0; l<nsoftrgna; l++) if (bw.helix_no == softrgn_allowed[l]) allowed = true;
                            }

                            if (allowed)
                            {
                                softrgns[nsoftrgn].rgn = r;
                                #if _dbg_soft
                                cout << "Region " << r.start << "-" << r.end << " bw " << bw.helix_no;
                                cout << endl;
                                #endif
                                nsoftrgn++;
                            }
                            #if _dbg_soft
                            else
                            {
                                cout << "Region " << r.start << "-" << r.end << " bw " << bw.helix_no << " is not allowed to soft pivot." << endl;
                            }
                            #endif
                        }
                        i=j;
                        break;
                    }
                }
            }
        }

        if (hg.n_helix)
        {
            for (i=0; i<hg.n_helix; i++)
            {
                if (hg.helices[i].n_ic)
                {
                    for (j=0; j<hg.helices[i].n_ic; j++)
                    {
                        InternalContact* lic = &(hg.helices[i].ic[j]);
                        lic->res1.resolve_resno(p);
                        lic->res2.resolve_resno(p);
                        if (lic->res1.resno && lic->res2.resno)
                        {
                            AminoAcid* aa1 = p->get_residue(lic->res1.resno);
                            AminoAcid* aa2 = p->get_residue(lic->res2.resno);

                            if (aa1 && aa2
                                && aa1->get_reach_atom(hbond)
                                && aa2->get_reach_atom(hbond)
                                )
                            {
                                aa1->movability = aa2->movability = MOV_PINNED;
                            }
                        }
                    }
                }
            }
        }

        for (i=0; i<n_soft_contact; i++)
        {
            soft_contact_a[i].resolve_resno(p);
            soft_contact_b[i].resolve_resno(p);

            AminoAcid* aa = p->get_residue(soft_contact_a[i].resno);
            if (aa) aa->movability = MOV_PINNED;
            aa = p->get_residue(soft_contact_b[i].resno);
            if (aa) aa->movability = MOV_PINNED;

            bool matched = false;

            for (j=0; j<nsoftrgn; j++)
            {
                if (soft_contact_a[i].resno >= softrgns[j].rgn.start && soft_contact_a[i].resno <= softrgns[j].rgn.end)
                {
                    softrgns[j].add_contact(soft_contact_a[i].resno, soft_contact_b[i].resno, p, matched);
                    if (j) softrgns[j].link_region(&softrgns[j-1]);
                    matched = true;
                }
                else if (soft_contact_b[i].resno >= softrgns[j].rgn.start && soft_contact_b[i].resno <= softrgns[j].rgn.end)
                {
                    softrgns[j].add_contact(soft_contact_a[i].resno, soft_contact_b[i].resno, p, matched);
                    if (j) softrgns[j].link_region(&softrgns[j-1]);
                    matched = true;
                }
            }

            softrgns[i].check_chain_constraints(p);                 // Refreshes AA pointers to current protein.
        }
    }

    protein->optimize_hydrogens();
    reconnect_bridges();
}

int main(int argc, char** argv)
{
    strcpy(splash, "\n        ***     ******      ****    ***    ***     ***     *****       ****      ****   **    **     \n       ** **    **   **    **  **   ****  ****    ** **    **  **     **  **    **  **  **    **     \n      **   **   **    **  **    **  ** **** **   **   **   **   **   **    **  **       **    **     \n     **     **  **   **   **    **  **  **  **  **     **  **    **  **    **  **       **   **      \n     *********  ******    **    **  **  **  **  *********  **    **  **    **  **       ******       \n     **     **  **   **   **    **  **      **  **     **  **    **  **    **  **       **   **      \n     **     **  **    **  **    **  **      **  **     **  **   **   **    **  **       **    **     \n     **     **  **    **   **  **   **      **  **     **  **  **     **  **    **  **  **    **     \n     **     **  **    **    ****    **      **  **     **  *****       ****      ****   **    **     \n");
    char buffer[65536];
    int i, j;

    _momentum_rad_ceiling = fiftyseventh * 5;

    for (i=0; i<65536; i++) buffer[i] = 0;

    for (i=0; i<256; i++)
        configfname[i] = protfname[i] = protafname[i] = ligfname[i] = cvtyfname[i] = 0;

    for (i=0; i<SPHREACH_MAX+100; i++)
        cfmol_known_good[i] = nullptr;

    time_t began = time(NULL);
    struct tm *lbegan = localtime(&began);
    if (lbegan->tm_mon == 3 && lbegan->tm_mday == lbegan->tm_mon-2) for (i=0; i<203; i+=102) for (j=69; j<78; j++) splash[i+j] = splash[i+j+408];

    strcpy(configfname, "example.config");

    smcmd = false;
    for (i=1; i<argc; i++)
    {
        if (argv[i][0] == '-' && argv[i][1] == '-')
        {
            argv[i] += 2;
            for (j=0; argv[i][j]; j++) if (argv[i][j] >= 'a' && argv[i][j] <= 'z') argv[i][j] &= 0x5f;
            j = interpret_config_line(&argv[i]);
            if (ligset) ligcmd = true;
            if (smset) smcmd = true;
            // if (optsecho.size()) cout << optsecho << endl;
            argv[i] -= 2;
            i += j;
        }
        else if (!configset && file_exists(argv[i]))
        {
            char* dot = strrchr(argv[i], '.');
            if (dot)
            {
                char* ext = &dot[1];
                if (!strcmp(ext, "config"))
                {
                    strcpy(configfname, argv[i]);
                    configset = true;
                }
            }
        }
    }

    FILE* pf = fopen(configfname, "r");
    if (!pf)
    {
        cout << "Config file not found: " << configfname << ", exiting." << endl;
        return 0xbadf12e;
    }

    read_config_file(pf);
    fclose(pf);

    char protid[255];
    char* slash = strrchr(protfname, '/');
    if (!slash) slash = strrchr(protfname, '\\');
    strcpy(protid, slash ? slash+1 : protfname );
    char* dot = strchr(protid, '.');
    if (dot) *dot = 0;

    Protein pose_proteins[poses];
    Molecule pose_ligands[poses+1];
    protein = &pose_proteins[0]; // new Protein(protid);
    pf = fopen(protfname, "r");
    if (!pf)
    {
        cout << splash << endl;
        cerr << "Error trying to read " << protfname << endl;
        return 0xbadf12e;
    }
    protein->load_pdb(pf, 0, protstrand ?: 'A');
    fclose(pf);

    j=1;
    std::string seq = protein->get_sequence();
    for (i=0; splash[i]; i++)
    {
        if (splash[i] == '*')
        {
            splash[i] = seq.c_str()[j++];
            if (j >= seq.length()) j = 0;
        }
        else if (splash[i] == ' ') splash[i] = '-';
    }

    #if _dbg_A100
    float init_A100 = protein->A100();
    #endif

    cout << splash << endl;

    for (i=1; i<argc; i++)
    {
        if (argv[i][0] == '-' && argv[i][1] == '-')
        {
            argv[i] += 2;
            for (j=0; argv[i][j]; j++) if (argv[i][j] >= 'a' && argv[i][j] <= 'z') argv[i][j] &= 0x5f;
            j = interpret_config_line(&argv[i]);
            if (optsecho.size()) cout << optsecho << endl;
            argv[i] -= 2;
            i += j;
        }
    }

    if (out_bb_pairs && pdpst != pst_best_binding)
    {
        cout << "Warning: OUTBBP is enabled but the search mode is not the best-binding algorithm." << endl;
    }

    char* pcntp = strstr(outfname, "%p");
    if (pcntp)
    {
        char tmp[512], protn[64];
        *(pcntp++) = 0;
        *(pcntp++) = 0;
        strcpy(protn, strrchr(protfname, '/')+1);
        char* dot = strchr(protn, '.');
        if (dot) *dot = 0;
        sprintf(tmp, "%s%s%s", outfname, protn, pcntp);
        strcpy(outfname, tmp);
    }

    char* pcntl = strstr(outfname, "%l");
    if (pcntl)
    {
        char tmp[512], lign[64];
        *(pcntl++) = 0;
        *(pcntl++) = 0;
        strcpy(lign, strrchr(ligfname, '/')+1);
        char* dot = strchr(lign, '.');
        if (dot) *dot = 0;
        sprintf(tmp, "%s%s%s", outfname, lign, pcntl);
        strcpy(outfname, tmp);
    }

    #if _DBG_SPACEDOUT
    cout << "Starting a file outstream: " << outfname << endl;
    #endif
    output = new std::ofstream(outfname, std::ofstream::out);
    if (!output) return -1;

    pre_ligand_flex_radius = search_size.magnitude();
    pre_ligand_multimol_radius = pre_ligand_flex_radius + (default_pre_ligand_multimol_radius - default_pre_ligand_flex_radius);

    #if _DBG_STEPBYSTEP
    if (debug) *debug << "Loaded config file." << endl;
    #endif

    if (kcal) kJmol_cutoff /= _kcal_per_kJ;
    drift = 1.0 / (iters/25+1);

    if (strlen(cvtyfname))
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

    apply_protein_specific_settings(protein);
    #if _DBG_STEPBYSTEP
    if (debug) *debug << "Loaded protein." << endl;
    #endif

    Pose tmp_pdb_residue[poses+1][protein->get_end_resno()+1];
    Pose tmp_pdb_waters[poses+1][omaxh2o+1];
    Pose tmp_pdb_ligand[poses+1];
    Point tmp_pdb_metal_locs[poses+1][nmtlcoords+1];

    if (tplset)
    {
        ptemplt = new Protein("template");
        pf = fopen(tplfname, "r");
        if (!pf)
        {
            cerr << "Error trying to read " << tplfname << endl;
            return 0xbadf12e;
        }
        ptemplt->load_pdb(pf);
        fclose(pf);

        #if _dbg_homology
        cout << "Homology template is " << tplfname << endl;
        #endif

        if (tprfset)
        {
            ptplref = new Protein("reference");
            pf = fopen(tplrfnam, "r");
            if (!pf)
            {
                cerr << "Error trying to read " << tplrfnam << endl;
                return 0xbadf12e;
            }
            ptplref->load_pdb(pf);
            fclose(pf);
            protein->homology_conform(ptemplt, ptplref);
        }
        else protein->homology_conform(ptemplt, protein);

        if (temp_pdb_file.length()) std::remove(temp_pdb_file.c_str());
        temp_pdb_file = (std::string)"tmp/" + std::to_string(pid) + (std::string)"_homolog.pdb";

        pf = fopen(temp_pdb_file.c_str(), "wb");
        protein->save_pdb(pf);
        fclose(pf);
    }

    if (hydrogenate_pdb)
    {
        int resno, endres = protein->get_end_resno();
        for (resno=1; resno<=endres; resno++)
        {
            AminoAcid* res = protein->get_residue(resno);
            if (res) res->hydrogenate();
        }

        if (temp_pdb_file.length()) std::remove(temp_pdb_file.c_str());
        temp_pdb_file = (std::string)"tmp/" + std::to_string(pid) + (std::string)"_hydro.pdb";

        pf = fopen(temp_pdb_file.c_str(), "wb");
        protein->save_pdb(pf);
        fclose(pf);
    }

    int l, j1, i2, miter;

    if (nmtlcoords)
    {
        protein->pocketcen = pocketcen;
        protein->coordinate_metal(mtlcoords, nmtlcoords);
        metald_prot = protein;

        if (temp_pdb_file.length()) std::remove(temp_pdb_file.c_str());
        temp_pdb_file = (std::string)"tmp/" + std::to_string(pid) + (std::string)"_metal.pdb";

        pf = fopen(temp_pdb_file.c_str(), "wb");
        protein->save_pdb(pf, protein->metals_as_molecule());
        fclose(pf);
    }

    if (bridges.size())
    {
        reconnect_bridges();

        if (temp_pdb_file.length()) std::remove(temp_pdb_file.c_str());
        temp_pdb_file = (std::string)"tmp/" + std::to_string(pid) + (std::string)"_bridged.pdb";

        pf = fopen(temp_pdb_file.c_str(), "wb");
        protein->save_pdb(pf);
        fclose(pf);
    }

    for (i=0; i<required_contacts.size(); i++)
    {
        required_contacts[i].resolve_resno(protein);
    }

    if (!CEN_buf.size())
    {
        cerr << "Error: no binding pocket centers defined." << endl;
        return 0xbadb19d;
    }

    #if pocketcen_is_loneliest
    pocketcen = loneliest;
    #endif

    if (wet_contact && !maxh2o) maxh2o = omaxh2o = 1;
    if (maxh2o > 0)
    {
        waters = new Molecule*[maxh2o+2];
        owaters = new Molecule*[maxh2o+2];
        for (i=0; i<maxh2o; i++)
        {
            waters[i] = new Molecule("H2O");
            waters[i]->from_smiles("O{O1}");
            owaters[i] = waters[i];

            int j;
            for (j=0; j<3; j++) strcpy(waters[i]->get_atom(j)->aa3let, "HOH");
        }
        waters[i] = nullptr;
        owaters[i] = nullptr;
    }

    if (waters)
    {
        float szscale = 0.5;
        for (i=0; i<maxh2o; i++)
        {
            waters[i]->recenter(Point(
                pocketcen.x + frand(-search_size.x * szscale, search_size.x * szscale),
                pocketcen.y + frand(-search_size.y * szscale, search_size.y * szscale),
                pocketcen.z + frand(-search_size.z * szscale, search_size.z * szscale)
            ));
        }
    }

    pktset = true;

    l=0;
    addl_resno[l] = 0;

    // Load the ligand or return an error.
    // Molecule m(ligfname);
    // ligand = &m;
    for (l=0; l<=poses; l++)
    {
        ligand = &pose_ligands[l];
        std::string ligname = get_fttl_no_path_or_ext(ligfname);
        std::string lligfname = ligfname;
        ligand->set_name(ligname.c_str());
        char* ext = get_file_ext(ligfname);
        if (!ext)
        {
            cout << "Ligand file is missing its extension! " << ligfname << endl;
            return 0xbadf12e;
        }

        for (i=0; i<65536; i++) buffer[i] = 0;

        size_t wgaf;
        if (smset) ligand->from_smiles(smiles);
        else switch (ext[0])
        {
        case 's':
        case 'S':
            // SDF
            if (isomers.size())
            {
                isono = rand() % isomers.size();
                ligname = get_fttl_no_path_or_ext(isomers[isono].c_str());
                ligand->set_name(ligname.c_str());
                lligfname = isomers[isono].c_str();
            }
            pf = fopen(lligfname.c_str(), "rb");
            if (!pf)
            {
                cerr << "Error trying to read " << lligfname << endl;
                return 0xbadf12e;
            }
            wgaf = fread(buffer, 1, 65535, pf);
            fclose(pf);
            ligand->from_sdf(buffer);
            if (ligand->get_atom_count() > 1 && ligand->get_bounding_box().size().magnitude() < 0.1)
            {
                cout << "Error in input file " << lligfname << endl;
                throw 0xbadda7a;
            }
            break;

        case 'p':
        case 'P':
            if (isomers.size())
            {
                i = rand() % isomers.size();
                ligname = get_fttl_no_path_or_ext(isomers[i].c_str());
                ligand->set_name(ligname.c_str());
                lligfname = isomers[i].c_str();
            }
            pf = fopen(lligfname.c_str(), "rb");
            if (!pf)
            {
                cerr << "Error trying to read " << lligfname << endl;
                return 0xbadf12e;
            }
            ligand->from_pdb(pf);
            fclose(pf);
            break;

        default:
            cout << "Unrecognized ligand file extension: " << ext << endl;
            return 0xbadf12e;
        }

        ligand->minimize_internal_clashes();
        if (estimate_copylig)
        {
            if (!ncvtys)
            {
                cerr << "Unable to estimate ligand multiplicity without a cavity file." << endl;
                throw 0xbadda7a;
            }
            float rcav = Avogadro;
            int cvi = 0;
            for (i=0; i<ncvtys; i++)
            {
                float lr = cvtys[i].get_center().get_3d_distance(pocketcen);
                if (!i || lr < rcav)
                {
                    rcav = lr;
                    cvi = i;
                }
            }

            float cv = cvtys[cvi].get_volume(), lv = ligand->get_volume();
            copylig = cv / lv / 6;
            // cout << "Ligand volume " << lv << " fits " << (cv/lv) << "x into cavity " << cv << endl << flush;
        }

        if (!l)
        {
            cout << "Ligand multiplicity: " << copylig << endl;
            if (output) *output << "Ligand multiplicity: " << copylig << endl;
        }
        if (copylig>1)
        {
            ligand->make_multimer(copylig);
        }
    }
    ligand = &pose_ligands[0];
    cout << endl;

    #if _DBG_STEPBYSTEP
    if (debug) *debug << "Loaded ligand." << endl;
    #endif

    Point box = ligand->get_bounding_box().size();

    if (debug) *debug << "Ligand bounding box size: " << box.printable() << endl;
    #if _DBG_STEPBYSTEP
    if (debug) *debug << "Ligand bounding box." << endl;
    #endif

    // Identify the ligand atom with the greatest potential binding.
    int k, n;

    // Best-Binding Code

    #if _DBG_STEPBYSTEP
    if (debug) *debug << "Identified best binding ligand atoms." << endl;
    #endif

    Point nodecen = pocketcen;
    seql = protein->get_seq_length();
    int rstart = protein->get_start_resno();

    // Filter residues according to which ones are close enough to the spheroid to "reach" it.
    nodecen = pocketcen;

    // When docking with a metalloprotein, use this temporary Molecule for interactions the same as
    // we use AminoAcid objects, except don't attempt to flex the metals object.
    Molecule* met = protein->metals_as_molecule();
    #if _DBG_STEPBYSTEP
    if (debug) *debug << "Created metals molecule." << endl;
    #endif

    float bclash = 0;

    for (l=0; l<seql; l++)
    {
        AminoAcid* laa = protein->get_residue(l+rstart);
        if (!laa) continue;
        for (n=l+1; n<seql; n++)
        {
            AminoAcid* naa = protein->get_residue(n+rstart);
            if (!naa) continue;
            bclash += laa->get_intermol_clashes(naa);
        }

        bclash += ligand->get_intermol_clashes(laa);

        if (met) bclash += laa->get_intermol_clashes(met);
    }
    if (met) bclash += ligand->get_intermol_clashes(met);
    if (debug) *debug << "Initial clashes: " << bclash << endl;

    // TODO: Output some basic stats: receptor, ligand, etc.
    cout << "PDB file: " << protfname << endl;
    if (output) *output << "PDB file: " << protfname << endl;
    cout << "Ligand: " << ligfname << endl;
    if (output) *output << "Ligand: " << ligfname << endl;
    cout << endl;
    if (output) *output << endl;

    i = poses*(triesleft+1)+8;
    j = pathnodes+2;
    DockResult dr[i][j];

    float rgnxform_r[i][pathnodes+2][PROT_MAX_RGN], rgnxform_theta[i][pathnodes+2][PROT_MAX_RGN], rgnxform_y[i][pathnodes+2][PROT_MAX_RGN];
    float rgnrot_alpha[i][pathnodes+2][PROT_MAX_RGN], rgnrot_w[i][pathnodes+2][PROT_MAX_RGN], rgnrot_u[i][pathnodes+2][PROT_MAX_RGN];

    int drcount = 0, qpr;

    if (pathnodes)
    {
        cout << pathnodes << " path node" << (pathnodes == 1 ? "" : "s") << "." << endl;
        if (output) *output << pathnodes << " path node" << (pathnodes == 1 ? "" : "s") << "." << endl;
    }
    else
    {
        cout << "Static dock - no path nodes." << endl;
        if (output) *output << "Static dock - no path nodes." << endl;
    }
    if (progressbar) cout << endl;

    reaches_spheroid = new AminoAcid**[pathnodes+2];
    for (i=0; i<=pathnodes; i++) reaches_spheroid[i] = new AminoAcid*[SPHREACH_MAX+4];

    found_poses = 0;
    int wrote_acvmx = -1, wrote_acvmr = -1;
    float l_atom_clash_limit = clash_limit_per_atom; // + kJmol_cutoff;

    if (pdpst == pst_external)
    {
        ligand = &pose_ligands[1];

        // Open file.
        FILE* fp = fopen(copyfrom_filename.c_str(), "rb");
        if (!fp)
        {
            cerr << "Failed to open pre-conformation file " << copyfrom_filename << " for reading." << endl;
            throw 0xbadf12e;
        }

        // Populate preconforms with data from file.
        npreconf = 0;
        while (!feof(fp))
        {
            char* prevent_useless_warning = fgets(buffer, 1024, fp);
            if (buffer[0] == 'A' && buffer[1] == 'T' && buffer[2] == 'O' && buffer[3] == 'M' && buffer[4] == ' ')
            {
                char aname[5];
                for (i=0; i<5; i++) aname[i] = 0;
                l = (buffer[12] == ' ') ? 13 : 12;
                for (i=0; i<4; i++)
                {
                    aname[i] = (buffer[l+i] == ' ') ? 0 : buffer[l+i];
                    if (!aname[i]) break;
                }

                if (!aname[0])
                {
                    cerr << "Bad atom name in preconform file." << endl;
                    throw 0xbadda7a;
                }

                Atom* lla = ligand->get_atom(aname);
                if (!lla)
                {
                    cerr << "Atom name " << aname << " in preconform file doesn't match ligand." << endl;
                    throw 0xbadda7a;
                }

                float x, y, z;
                buffer[53] = 0;
                z = atof(buffer+45);
                buffer[45] = 0;
                y = atof(buffer+37);
                buffer[37] = 0;
                x = atof(buffer+29);

                lla->move(Point(x, y, z));
            }
            else if (buffer[0] == 'T' && buffer[1] == 'E' && buffer[2] == 'R')
            {
                preconforms[npreconf].copy_state(ligand);
                npreconf++;
            }
        }

        fclose(fp);

        if (!npreconf)
        {
            cerr << "External preconf file has no poses." << endl;
            throw 0xbadda7a;
        }
    }


_try_again:
    srand(time(NULL));
    Point nodecens[pathnodes+1];

    /////////////////////////////////////////////////////////////////////////////////
    // Main loop.
    /////////////////////////////////////////////////////////////////////////////////

    regions = protein->get_regions();
    for (i=0; regions[i].start; i++)
    {
        region_clashes[i][0] = region_clashes[i][1] = region_clashes[i][2] = Vector(0,0,0);
    }

    n = nmtlcoords;
    Point metal_initlocs[n+4];
    for (i=0; i<n; i++)
    {
        metal_initlocs[i] = mtlcoords[i].mtl->loc;
    }

    float best_energy = 0, best_acc_energy = 0, best_worst_clash = 0;
    for (pose = 1; pose <= poses; pose++)
    {
        soft_contact_elasticity = initial_soft_contact_elasticity;
        ligand = &pose_ligands[pose];
        ligand->movability = MOV_ALL;
        ligand->minimize_internal_clashes();
        if (!ligand->is_chiral() && frand(0,1) < 0.5) ligand->mirror();

        memset(g_bbr, 0, sizeof(BestBindingResult)*16);

        // TODO: Update atom & res pointers for global bb pairs.

        last_ttl_bb_dist = 0;
        ligand->minimize_internal_clashes();
        float lig_min_int_clsh = ligand->get_internal_clashes();

        int rcn = required_contacts.size();
        if (rcn)
        {
            ligand->delete_mandatory_connections();
            ligand->allocate_mandatory_connections(rcn);

            for (i=0; i<rcn; i++)
            {
                if (required_contacts[i].node >= 0 && required_contacts[i].node != nodeno) continue;
                Star s;
                s.paa = protein->get_residue(required_contacts[i].resno);
                if (s.n) ligand->add_mandatory_connection(s.pmol);
            }
        }

        // delete protein;
        // protein = new Protein(protfname);
        protein = &pose_proteins[pose-1];

        if (temp_pdb_file.length())
        {
            pf = fopen(temp_pdb_file.c_str(), "r");
            protein->load_pdb(pf);
            fclose(pf);
            apply_protein_specific_settings(protein);

            if (nmtlcoords)
            {
                for (i=0; i<nmtlcoords; i++)
                {
                    for (j=0; j<mtlcoords[i].ncoordres; j++)
                    {
                        AminoAcid* aa = protein->get_residue(mtlcoords[i].coordres[j].resno);
                        if (aa)
                        {
                            aa->coordmtl = mtlcoords[i].mtl;
                        }
                    }
                }
            }
        }
        else
        {
            pf = fopen(protfname, "r");
            protein->load_pdb(pf, 0, protstrand ?: 'A');
            fclose(pf);
            apply_protein_specific_settings(protein);
        }

        n = nmtlcoords;
        for (i=0; i<n; i++)
        {
            mtlcoords[i].mtl->move(metal_initlocs[i]);
        }

        freeze_bridged_residues();

        ligand->recenter(pocketcen);
        // cout << "Centered ligand at " << pocketcen << endl << endl << flush;

        if (pdpst == pst_tumble_spheres)
        {
            Search::do_tumble_spheres(protein, ligand, pocketcen);
            attempt_priority_hbond();
            pocketcen = ligand->get_barycenter();

            #if debug_stop_after_tumble_sphere
            return 0;
            #endif
        }

        #if _DBG_STEPBYSTEP
        if (debug) *debug << "Pose " << pose << endl;
        #endif
        nodecen = pocketcen;
        nodecen.weight = 1;

        #if _dummy_atoms_for_debug
        dummies.clear();
        #endif

        for (nodeno=0; nodeno<=pathnodes; nodeno++)
        {
            soft_contact_elasticity = initial_soft_contact_elasticity;
            movie_offset = iters * (nodeno /*+ (pose-1)*(pathnodes+1)*/);

            if (waters)
            {
                for (i = 0; i <= omaxh2o; i++)
                {
                    waters[i] = owaters[i];
                }
                maxh2o = omaxh2o;
            }

            // if (pathstrs.size() < nodeno) break;
            drift = initial_drift;

            if (echo_progress) cout << (time(NULL) - began) << " seconds: starting pose " << pose << " node " << nodeno << "..." << endl;

            #if internode_momentum_only_on_activation 
            conformer_momenta_multiplier = 1;
            #else
            conformer_momenta_multiplier = nodeno ? internode_momentum_mult : 1;
            #endif
            conformer_tumble_multiplier = 1;

            allow_ligand_360_tumble = nodes_no_ligand_360_tumble && pdpst != pst_best_binding;
            allow_ligand_360_flex   = nodes_no_ligand_360_flex;

            if (strlen(protafname) && nodeno == activation_node)
            {
                #if internode_momentum_only_on_activation 
                conformer_momenta_multiplier = nodeno ? internode_momentum_mult : 1;
                #endif

                #if prevent_ligand_360_on_activate
                allow_ligand_360_tumble = allow_ligand_360_flex = false;
                #endif

                // Persist the flexions of the side chains. 
                // TODO: Do not persist those residues whose positions are important to activation.
                float* sidechain_bondrots[seql+4];
                int sidechain_bondrotq[seql+4];
                for (i=0; i<seql+4; i++)
                {
                    sidechain_bondrots[i] = nullptr;
                    sidechain_bondrotq[i] = 0;
                }
                for (i=1; i<=seql; i++)
                {
                    Bond** b = protein->get_residue(i)->get_rotatable_bonds();
                    if (b)
                    {
                        int bq;
                        for (bq=0; b[bq]; bq++);                // Get count.
                        sidechain_bondrots[i] = new float[bq+2];
                        for (j=0; j<bq; j++)
                        {
                            sidechain_bondrots[i][j] = b[j]->total_rotations;
                        }
                        sidechain_bondrotq[i] = j;

                        // delete[] b;
                    }
                }

                // delete protein;
                // protein = new Protein(protafname);
                protein = &pose_proteins[pose-1];

                pf = fopen(protafname, "r");
                protein->load_pdb(pf);
                fclose(pf);
                apply_protein_specific_settings(protein);

                freeze_bridged_residues();

                for (i=1; i<=seql; i++)
                {
                    Bond** b = protein->get_residue(i)->get_rotatable_bonds();
                    if (b)
                    {
                        int bq;

                        for (j=0; j<sidechain_bondrotq[i]; j++)
                        {
                            if (!b[j]) break;
                            b[j]->clear_moves_with_cache();
                            b[j]->rotate(sidechain_bondrots[i][j]);
                        }

                        // delete[] b;
                    }

                    delete sidechain_bondrots[i];
                }
            }

            #if _DBG_STEPBYSTEP
            if (debug) *debug << "Pose " << pose << endl << "Node " << nodeno << endl;
            #endif
            if (nodeno)
            {
                for (i=0; i<states.size(); i++)
                {
                    strcpy(buffer, states[i].c_str());
                    char** words = chop_spaced_words(buffer);
                    if (atoi(words[1]) == nodeno)
                    {
                        int sr = atoi(words[2]), er = atoi(words[3]);
                        float theta = atof(words[4]) * fiftyseventh;

                        Point sloc = protein->get_atom_location(sr, "CA"),
                              eloc = protein->get_atom_location(er, "CA");

                        LocatedVector lv = (Vector)(sloc.subtract(eloc));
                        lv.origin = sloc;

                        int resno;
                        for (resno = sr; resno <= er; resno++)
                        {
                            AminoAcid* aa = protein->get_residue(resno);
                            if (aa)
                            {
                                MovabilityType mt = aa->movability;
                                aa->movability = MOV_ALL;

                                aa->rotate(lv, theta);

                                aa->movability = mt;
                            }
                        }
                    }
                }

                // nodecen = nodecen.add(&path[nodeno]);
                strcpy(buffer, pathstrs[nodeno].c_str());
                if (!strlen(buffer))
                {
                    cerr << "Error in config file: path node " << nodeno << " is missing." << endl;
                    return 0xbadc09f;
                }
                char** words = chop_spaced_words(buffer);
                nodecen = pocketcen_from_config_words(&words[1], &nodecen);

                #if _DBG_STEPBYSTEP
                if (debug) *debug << "Added whatever points together." << endl;
                #endif
            }

            loneliest = protein->find_loneliest_point(nodecen, search_size);
            loneliest = loneliest_from_cavities(loneliest);
            // cout << "Loneliest is " << loneliest << endl;

            #if pocketcen_is_loneliest
            if (!pathstrs.size()) nodecen = loneliest;
            #endif

            Point lastnodecen = nodecen;
            ligcen_target = nodecen;
            nodecens[nodeno] = ligcen_target;

            #if redo_tumble_spheres_every_node

            if (pdpst == pst_tumble_spheres && (!prevent_ligand_360_on_activate))
            {
                Search::do_tumble_spheres(protein, ligand, ligcen_target);
                attempt_priority_hbond();
            }
            #endif

            #if _DBG_STEPBYSTEP
            if (debug) *debug << "Saved last nodecen." << endl;
            #endif

            #if recenter_ligand_each_node
            // Move the ligand to the new node center.
            ligand->recenter(nodecen);
            #if _DBG_STEPBYSTEP
            if (debug) *debug << "Molecule recenter (or not)." << endl;
            #endif
            ligand->reset_conformer_momenta();
            #if _DBG_STEPBYSTEP
            if (debug) *debug << "Conformer momenta reset." << endl;
            #endif
            #endif

            #if _dbg_groupsel
            cout << "Priority resnos: ";
            #endif

            if (n = priority_resnos.size()) for (i=0; i<n; i++)
            {
                AminoAcid* aa = protein->get_residue(priority_resnos[i]);
                if (aa) aa->priority = true;
                #if _dbg_groupsel
                cout << *aa << " ";
                #endif
            }

            #if _dbg_groupsel
            cout << endl << endl;
            #endif

            Point lsz = search_size;
            lsz.multiply(0.75);
            sphres = protein->get_residues_can_clash_ligand(reaches_spheroid[nodeno], ligand, nodecen, lsz, addl_resno);
            if (!sphres) cerr << "Warning: No binding pocket residues found." << endl;
            for (i=sphres; i<SPHREACH_MAX; i++) reaches_spheroid[nodeno][i] = NULL;

            #if point_hydro_side_chains_inward
            for (i=0; i<sphres; i++)
            {
                AminoAcid* aa = reaches_spheroid[nodeno][i];
                if (!aa) continue;
                if (aa->movability == MOV_PINNED) continue;
                bool aa_cannot_point_inward = false;
                for (j=0; j<n_soft_contact; j++)
                {
                    if (soft_contact_a[j].resno == aa->get_residue_no()) aa_cannot_point_inward = true;
                    if (soft_contact_b[j].resno == aa->get_residue_no()) aa_cannot_point_inward = true;
                    if (aa_cannot_point_inward) break;
                }
                if (aa_cannot_point_inward) continue;
                if (fabs(aa->hydrophilicity()) < hydrophilicity_cutoff) continue;
                // if (aa->get_charge() < -0.1) continue;                           // better of with this or without it? ¯\_(ツ)_/¯

                Atom* a = aa->get_reach_atom(hbond);
                if (!a) continue;
                aa->movability = MOV_FLEXONLY;
                Pose was(aa);
                AminoAcid** aaclashables = protein->clashable_residues(aa->get_residue_no());
                Interaction e = aa->get_intermol_binding(aaclashables);
                aa->conform_atom_to_location(a->name, loneliest);
                aa->conform_molecules((Molecule**)aaclashables, nullptr, 20);
                Interaction f = aa->get_intermol_binding(protein->clashable_residues(aa->get_residue_no()));
                if (f.clash > clash_limit_per_aa*2 || f.summed() > e.summed()+clash_limit_per_aa*2) was.restore_state(aa);
            }
            #endif

            #if _dbg_bb_scoring
            if (pose==1)
            {
                cout << endl << "BSR:";
                for (i=0; i<sphres; i++) cout << " " << reaches_spheroid[nodeno][i]->get_name();
                cout << endl << endl;
            }
            #endif

            if (nsoftrgn)
            {
                for (i=0; i<nsoftrgn; i++)
                {
                    int nc = softrgns[i].num_contacts();
                    for (j=0; j<nc; j++)
                    {
                        AminoAcid *aac1 = softrgns[i].get_local_contact(j, protein), *aac2 = softrgns[i].get_distant_contact(j, protein);
                        if (!aac1 || !aac2) continue;
                        aac1->is_ic_res = true;
                        aac2->is_ic_res = true;
                    }
                }
            }

            // Flexion Selection
            if (flex && !nodeno)
            {
                if (forced_static_resnos.size())
                {
                    for (i=0; i<forced_static_resnos.size(); i++)
                    {
                        forced_static_resnos[i].resolve_resno(protein);
                        AminoAcid* mvaa = protein->get_residue(forced_static_resnos[i].resno);
                        if (mvaa)
                        {
                            mvaa->movability = MOV_PINNED;
                            #if _dbg_flexion_selection
                            cout << mvaa->get_name() << " forced static." << endl;
                            #endif
                        }
                    }
                }
                #if flexion_selection

                flexible_resnos.clear();
                j = protein->get_end_resno();
                for (i=protein->get_start_resno(); i<=j; i++)
                {
                    AminoAcid* mvaa = protein->get_residue(i);
                    if (mvaa && !(mvaa->movability & MOV_PINNED)) mvaa->movability = min(MOV_FLXDESEL, mvaa->movability);
                }

                #if _dbg_null_flexions
                bool another_flex = false;
                #elif no_zero_flexions
                bool another_flex = true;
                #else
                bool another_flex = (frand(0,1) < 0.6);
                #endif

                while (another_flex)
                {
                    float bestwt = 0;
                    int besti = -1;
                    for (j=0; j<100; j++)
                        for (i=0; i<sphres; i++)
                        {
                            if (reaches_spheroid[nodeno][i]->movability != MOV_FLXDESEL) continue;
                            float weight = reaches_spheroid[nodeno][i]->get_aa_definition()->flexion_probability;
                            if (!weight) continue;

                            if (reaches_spheroid[nodeno][i]->is_ic_res) weight /= 2;

                            // Multiply weight by unrealized ligand binding potential.
                            float potential = reaches_spheroid[nodeno][i]->get_intermol_potential(ligand, true);
                            float adjusted_potential = fmin(1, potential / 1000);
                            Atom *nearest1, *nearest2;
                            nearest1 = reaches_spheroid[nodeno][i]->get_nearest_atom(ligand->get_barycenter());
                            if (!nearest1) throw 0xbadd157;
                            nearest2 = ligand->get_nearest_atom(nearest1->loc);
                            if (!nearest2) throw 0xbadd157;
                            float nearr = fmax(1, nearest1->distance_to(nearest2) / 2);
                            adjusted_potential *= nearr;

                            // weight = (1.0 - ((1.0 - weight) / adjusted_potential)) / 2;
                            weight *= sqrt(adjusted_potential);

                            #if _dbg_flexion_selection
                            if (reaches_spheroid[nodeno][i]->get_residue_no() == 9262)
                                cout << reaches_spheroid[nodeno][i]->get_name() << " has weight " << weight << endl;
                            #endif

                            if ( /*weight >= bestwt &&*/ frand(0,100) < weight )
                            {
                                besti = i;
                                bestwt = weight;
                            }
                        }
                    if (besti >= 0)
                    {
                        if (!(reaches_spheroid[nodeno][besti]->movability & MOV_PINNED))
                            reaches_spheroid[nodeno][besti]->movability = MOV_FLEXONLY;
                        flexible_resnos.push_back(reaches_spheroid[nodeno][besti]->get_residue_no());
                        #if _dbg_flexion_selection
                        cout << "Selected " << reaches_spheroid[nodeno][besti]->get_name()
                                << " for flexion with a weight of " << bestwt << endl;
                        #endif
                    }
                    another_flex = (frand(0,1) < 0.6);
                }
                
                #if _dbg_flexion_selection
                cout << flexible_resnos.size() << " residues selected for flexion." << endl;
                #endif

                freeze_bridged_residues();

                if (forced_flexible_resnos.size())
                {
                    for (i=0; i<forced_flexible_resnos.size(); i++)
                    {
                        forced_flexible_resnos[i].resolve_resno(protein);
                        AminoAcid* mvaa = protein->get_residue(forced_flexible_resnos[i].resno);
                        if (mvaa)
                        {
                            mvaa->movability = MOV_FORCEFLEX;
                            flexible_resnos.push_back(mvaa->get_residue_no());
                            #if _dbg_flexion_selection
                            cout << mvaa->get_name() << " forced flexible." << endl;
                            #endif
                        }
                    }
                }
                #endif
            }

            #if _dbg_flexion_selection
            cout << endl;
            #endif

            #if _dbg_groupsel
            cout << "Candidate binding residues: ";
            for (i=0; i<sphres; i++)
            {
                cout << *reaches_spheroid[nodeno][i];
                if (reaches_spheroid[nodeno][i]->priority) cout << "!";
                cout << " ";
            }
            cout << endl;
            #endif

            for (i=0; i<_INTER_TYPES_LIMIT; i++) total_binding_by_type[i] = 0;
            #if use_exclusions
            for (i=0; i<sphres; i++)
            {
                if (exclusion.size()
                        &&
                        std::find(exclusion.begin(), exclusion.end(), reaches_spheroid[nodeno][i]->get_residue_no())!=exclusion.end()
                   )
                {
                    for (j=i; j<sphres; j++) reaches_spheroid[nodeno][j] = reaches_spheroid[nodeno][j+1];
                    sphres--;
                    reaches_spheroid[nodeno][sphres] = nullptr;
                }
            }
            if (!sphres) cerr << "Warning: No binding pocket residues left after applying exclusions." << endl;
            #endif

            // float driftamt = 1.0 / (iters/25+1);
            // cout << pose << ":" << nodeno << " drift " << driftamt << endl;
            int iters_div = iters*0.259;

            Molecule* cfmols[SPHREACH_MAX+4];
            for (i=0; i<=SPHREACH_MAX; i++) cfmols[i] = nullptr;
            gcfmols = cfmols;
            i=0;
            ligand->movability = MOV_ALL;
            cfmols[i++] = ligand;
            if (met)
            {
                met->movability = MOV_NONE;
                cfmols[i++] = met;
            }
            if (waters)
            {
                for (j=0; j<maxh2o; j++)
                {
                    waters[j]->movability = MOV_ALL;
                    waters[j]->reset_conformer_momenta();
                    cfmols[i++] = waters[j];
                }
            }

            if (n = priority_resnos.size())
            {
                for (j=0; j<n; j++)
                {
                    AminoAcid* aa = protein->get_residue(priority_resnos[j]);
                    if (aa) cfmols[i++] = (Molecule*)aa;
                }
            }


            for (j=0; j<sphres; j++)
            {
                #if ! flexion_selection
                if (reaches_spheroid[nodeno][j]->movability >= MOV_FLEXONLY && !(reaches_spheroid[nodeno][j]->movability & MOV_PINNED) ) 
                    reaches_spheroid[nodeno][j]->movability = MOV_FLEXONLY;
                #endif
                if (!flex) reaches_spheroid[nodeno][j]->movability = MOV_PINNED;
                cfmols[i++] = reaches_spheroid[nodeno][j];
                protein->get_residues_can_clash(reaches_spheroid[nodeno][j]->get_residue_no());      // This sets clashables internally to the protein.
            }

            int cfmolqty = i;
            for (; i<=SPHREACH_MAX; i++) cfmols[i] = NULL;


            if (!nodeno)
            {
                float alignment_distance[5];
                for (l=0; l<3; l++)
                {
                    alignment_distance[l]=0;
                }

                #if _DBG_STEPBYSTEP
                if (debug) *debug << "Initialize null AA pointer." << endl;
                #endif

                protein->set_conditional_basicities();

                #if reuse_best_pose
                if (aa_best_pose && ligand_best_pose[isono].has_data() && frand(0,1) <= reuse_pose_probability)
                {
                    n = protein->get_end_resno();
                    for (i=1; i<n; i++)
                    {
                        AminoAcid* aa = protein->get_residue(i);
                        if (!aa) continue;
                        aa_best_pose[i+n*isono].restore_state_relative(aa, "CA");
                    }

                    ligand_best_pose[isono].restore_state(ligand);
                    ligand->stay_close_mine = best_ligand_stays;
                    ligand->stay_close_other = best_ligand_stays_other;
                }
                else if (pdpst == pst_best_binding)
                #else
                if (pdpst == pst_best_binding)
                #endif
                {
                    int maxlt = ligand->get_heavy_atom_count()+8;
                    g_ligtargs = new LigandTarget[maxlt];
                    int ltargs = Search::identify_ligand_pairing_targets(ligand, g_ligtargs, maxlt);

                    if (ltargs > 0)
                    {
                        if (pose==1) cout << "Performing Best-Binding search..." << endl << flush;

                        Cavity* cv = nullptr;
                        float cvr = Avogadro;
                        if (ncvtys) for (i=0; i<ncvtys; i++)
                        {
                            float r = cvtys[i].get_center().get_3d_distance(ligand->get_barycenter());
                            if (r < cvr)
                            {
                                cv = &cvtys[i];
                                cvr = r;
                            }
                        }

                        Molecule* llig;
                        j = 0;
                        AminoAcid** lrs = new AminoAcid*[SPHREACH_MAX];
                        memcpy(lrs, reaches_spheroid[nodeno], sizeof(AminoAcid**)*SPHREACH_MAX);
                        #if _dbg_bb_scoring
                        cout << "Node " << nodeno << "; sphres = " << sphres << endl;
                        #endif
                        Point searchcen = pathnodes ? nodecen : loneliest;
                        for (l = 0; !l || l < copylig; l++)
                        {
                            llig = ligand->get_monomer(l);
                            LigandTarget llt[maxlt];
                            Search::identify_ligand_pairing_targets(llig, llt, maxlt);
                            upivsmfpoy:
                            try
                            {
                                Search::pair_targets(protein, llig, llt, lrs, searchcen, &g_bbr[l], cv, !l);
                            }
                            catch (int sodivld)
                            {
                                if (sodivld == 0x90ba1125)
                                {
                                    if (pose == 1)
                                    {
                                        cout << "FATAL ERROR." << endl;
                                        throw sodivld;
                                    }
                                    else goto upivsmfpoy;           // try again, we already know it can be done.
                                }
                                else throw sodivld;
                            }
                            g_bbr[l].protein = protein;

                            if (g_bbr[l].pri_res)
                            {
                                for (i=0; lrs[i]; i++)
                                {
                                    if (lrs[i] == g_bbr[l].pri_res)
                                    {
                                        for (; lrs[i]; i++) lrs[i] = lrs[i+1];
                                    }
                                }
                            }

                            if (g_bbr[l].pri_tgt)
                            {
                                g_ligtargs[j] = *g_bbr[l].pri_tgt;      // memberwise copy
                                g_bbr[l].pri_tgt = &g_ligtargs[j];      // repoint
                                j++;
                            }

                            if (g_bbr[l].sec_res)
                            {
                                for (i=0; lrs[i]; i++)
                                {
                                    if (lrs[i] == g_bbr[l].sec_res)
                                    {
                                        for (; lrs[i]; i++) lrs[i] = lrs[i+1];
                                    }
                                }
                            }

                            if (g_bbr[l].sec_tgt)
                            {
                                g_ligtargs[j] = *g_bbr[l].sec_tgt;      // memberwise copy
                                g_bbr[l].sec_tgt = &g_ligtargs[j];      // repoint
                                j++;
                            }

                            if (g_bbr[l].tert_res)
                            {
                                for (i=0; lrs[i]; i++)
                                {
                                    if (lrs[i] == g_bbr[l].tert_res)
                                    {
                                        for (; lrs[i]; i++) lrs[i] = lrs[i+1];
                                    }
                                }
                            }

                            if (g_bbr[l].tert_tgt)
                            {
                                g_ligtargs[j] = *g_bbr[l].tert_tgt;     // memberwise copy
                                g_bbr[l].tert_tgt = &g_ligtargs[j];     // repoint
                                j++;
                            }

                            if ((g_bbr[l].pri_res->has_hbond_acceptors() && g_bbr[l].pri_tgt->has_hb_donors())
                                    || (g_bbr[l].pri_res->has_hbond_donors() && g_bbr[l].pri_tgt->has_hb_acceptors())
                                    || (g_bbr[l].pri_res->coordmtl)
                                    )
                            {
                                if (g_bbr[l].pri_tgt->single_atom) llig->stay_close_mine = g_bbr[l].pri_tgt->single_atom;
                                else if (g_bbr[l].pri_tgt->conjgrp) llig->stay_close_mine = g_bbr[l].pri_tgt->conjgrp->get_nearest_atom(g_bbr[l].pri_res->get_CA_location());
                                llig->stay_close_mol = g_bbr[l].pri_res;
                                if (g_bbr[l].pri_res->coordmtl) llig->stay_close_other = g_bbr[l].pri_res->coordmtl;
                                else llig->stay_close_other = g_bbr[l].pri_res->get_reach_atom(g_bbr[l].pri_tgt->best_interaction());
                                if (!llig->stay_close_other) llig->stay_close_other = g_bbr[l].pri_res->get_nearest_atom(g_bbr[l].pri_tgt->barycenter());
                                if (llig->stay_close_other) llig->stay_close_optimal = InteratomicForce::optimal_distance(llig->stay_close_mine, llig->stay_close_other);
                                llig->stay_close_tolerance = 0.5;

                                if (pose <= 1 && llig->stay_close_other) cout << "Staying " << llig->stay_close_mine->name << " near " 
                                    << g_bbr[l].pri_res->get_residue_no() << ":" << llig->stay_close_other->name << endl;
                            }

                            #if bb_secondary_must_be_farthest_from_primary
                            if (g_bbr[l].tert_res && g_bbr[l].tert_tgt
                                && ((g_bbr[l].tert_res->has_hbond_acceptors() && g_bbr[l].tert_tgt->has_hb_donors())
                                    || (g_bbr[l].tert_res->has_hbond_donors() && g_bbr[l].tert_tgt->has_hb_acceptors())
                                    || (g_bbr[l].tert_res->coordmtl)
                                    )
                                )
                            {
                                if (g_bbr[l].tert_tgt->single_atom) llig->stay_close2_mine = g_bbr[l].tert_tgt->single_atom;
                                else if (g_bbr[l].tert_tgt->conjgrp) llig->stay_close2_mine = g_bbr[l].tert_tgt->conjgrp->get_nearest_atom(g_bbr[l].tert_res->get_CA_location());
                                llig->stay_close2_mol = g_bbr[l].tert_res;
                                if (g_bbr[l].tert_res->coordmtl) llig->stay_close2_other = g_bbr[l].tert_res->coordmtl;
                                else llig->stay_close2_other = g_bbr[l].tert_res->get_reach_atom(g_bbr[l].tert_tgt->best_interaction());
                                if (!llig->stay_close2_other) llig->stay_close2_other = g_bbr[l].tert_res->get_nearest_atom(g_bbr[l].tert_tgt->barycenter());
                                if (llig->stay_close2_other) llig->stay_close2_optimal = InteratomicForce::optimal_distance(llig->stay_close2_mine, llig->stay_close2_other);

                                if (pose <= 1 && llig->stay_close2_other) cout << "Staying " << llig->stay_close2_mine->name << " near " 
                                    << g_bbr[l].tert_res->get_residue_no() << ":" << llig->stay_close2_other->name << endl;
                            }
                            else
                            #endif
                            if (g_bbr[l].sec_res && g_bbr[l].sec_tgt
                                && ((g_bbr[l].sec_res->has_hbond_acceptors() && g_bbr[l].sec_tgt->has_hb_donors())
                                    || (g_bbr[l].sec_res->has_hbond_donors() && g_bbr[l].sec_tgt->has_hb_acceptors())
                                    || (g_bbr[l].sec_res->coordmtl)
                                    )
                                )
                            {
                                if (g_bbr[l].sec_tgt->single_atom) llig->stay_close2_mine = g_bbr[l].sec_tgt->single_atom;
                                else if (g_bbr[l].sec_tgt->conjgrp) llig->stay_close2_mine = g_bbr[l].sec_tgt->conjgrp->get_nearest_atom(g_bbr[l].sec_res->get_CA_location());
                                llig->stay_close2_mol = g_bbr[l].sec_res;
                                if (g_bbr[l].sec_res->coordmtl) llig->stay_close2_other = g_bbr[l].sec_res->coordmtl;
                                else llig->stay_close2_other = g_bbr[l].sec_res->get_reach_atom(g_bbr[l].sec_tgt->best_interaction());
                                if (!llig->stay_close2_other) llig->stay_close2_other = g_bbr[l].sec_res->get_nearest_atom(g_bbr[l].sec_tgt->barycenter());
                                if (llig->stay_close2_other) llig->stay_close2_optimal = InteratomicForce::optimal_distance(llig->stay_close2_mine, llig->stay_close2_other);

                                if (pose <= 1 && llig->stay_close2_other) cout << "Staying " << llig->stay_close2_mine->name << " near " 
                                    << g_bbr[l].sec_res->get_residue_no() << ":" << llig->stay_close2_other->name << endl;
                            }

                            if (llig->stay_close2_mine && !llig->stay_close_mine)
                            {
                                llig->stay_close_mine = llig->stay_close2_mine;
                                llig->stay_close_mol = llig->stay_close2_mol;
                                llig->stay_close_optimal = llig->stay_close2_optimal;
                                llig->stay_close_other = llig->stay_close2_other;
                                llig->stay_close2_mine = nullptr;
                                llig->stay_close2_other = nullptr;
                            }

                            // ligand->propagate_stays();
                            g_bbr[l].protein = protein;
                            Search::align_targets(llig, searchcen, &g_bbr[l]);
                            searchcen = searchcen.add(searchcen.subtract(llig->get_barycenter()));
                        }
                    }
                    else ligand->recenter(loneliest);

                    abhor_vacuum(0, (Molecule**)reaches_spheroid[nodeno]);

                    if (out_bb_pairs && ltargs > 1 && pose<=1)
                    {
                        for (l=0; g_bbr[l].pri_tgt; l++)
                        {
                            if (g_bbr[l].pri_tgt && g_bbr[l].pri_res)
                            {
                                cout << "Primary target: " << *g_bbr[l].pri_tgt << "..." << *g_bbr[l].pri_res << endl;
                                if (output) *output << "Primary target: " << *g_bbr[l].pri_tgt << "..." << *g_bbr[l].pri_res << endl;
                            }
                            if (g_bbr[l].sec_tgt && g_bbr[l].sec_res)
                            {
                                cout << "Secondary target: " << *g_bbr[l].sec_tgt << "..." << *g_bbr[l].sec_res << endl;
                                if (output) *output << "Secondary target: " << *g_bbr[l].sec_tgt << "..." << *g_bbr[l].sec_res << endl;
                            }
                            if (g_bbr[l].tert_tgt && g_bbr[l].tert_res)
                            {
                                cout << "Tertiary target: " << *g_bbr[l].tert_tgt << "..." << *g_bbr[l].tert_res << endl;
                                if (output) *output << "Tertiary target: " << *g_bbr[l].tert_tgt << "..." << *g_bbr[l].tert_res << endl;
                            }
                            cout << endl << endl << flush;
                            if (output) *output << endl;
                        }
                    }
                }
                else if (pdpst == pst_copyfrom)
                {
                    Search::copy_ligand_position_from_file(protein, ligand, copyfrom_filename.c_str(), copyfrom_ligname, copyfrom_resno);
                }
                else if (pdpst == pst_external)
                {
                    i = (pose-1) % npreconf;
                    preconforms[i].restore_state(ligand);
                }

                // else ligand->recenter(ligcen_target);

                // Best-Binding Algorithm
                // Find a binding pocket feature with a strong potential binding to the ligand.
                std::string alignment_name = "";
                
                #if _DBG_STEPBYSTEP
                if (debug) *debug << "Selected an alignment AA." << endl;
                #endif

                #if _DBG_STEPBYSTEP
                if (debug) *debug << "Alignment AA." << endl;
                #endif

                #if _DBG_STEPBYSTEP
                if (debug) *debug << "Aligned ligand to AA." << endl;
                cout << endl;
                #endif

                freeze_bridged_residues();
            }

            ligand->reset_conformer_momenta();

            if (rcn)
            {
                for (i=0; i<rcn; i++)
                {
                    Star s;
                    s.paa = protein->get_residue(required_contacts[i].resno);
                    if (required_contacts[i].node >= 0 && required_contacts[i].node != nodeno) ligand->remove_mandatory_connection(s.pmol);
                    else ligand->add_mandatory_connection(s.pmol);
                }
            }

            protein->find_residue_initial_bindings();
            freeze_bridged_residues();

            n = nmtlcoords;
            if (n)
            {
                int j1;
                for (j=0; j<n; j++)
                {
                    int n1 = mtlcoords[j].ncoordres;
                    for (j1=0; j1<n1; j1++)
                    {
                        AminoAcid* aa = protein->get_residue(mtlcoords[j].coordres[j1].resno);
                        if (aa) aa->movability = MOV_PINNED;
                    }
                }
            }

            #if _dbg_improvements_only_rule
            check_ligand = ligand;
            #endif

            ligand->movability = MOV_ALL;
            if (!flex) for (j=0; j<sphres; j++)
            {
                if (!(reaches_spheroid[nodeno][j]->movability & MOV_PINNED)) reaches_spheroid[nodeno][j]->movability = MOV_FLXDESEL;
            }
            freeze_bridged_residues();
            if (output_each_iter) output_iter(0, cfmols, "initial placement");
            if (pdpst == pst_best_binding) ligand->movability = (MovabilityType)(MOV_CAN_AXIAL | MOV_CAN_RECEN | MOV_CAN_FLEX);

            freeze_bridged_residues();
            if (auditfn.length())
            {
                audit = fopen(auditfn.c_str(), (pose==1 && !nodeno) ? "wb" : "ab");
                if (!audit) cerr << "Warning: failed to open " << auditfn << " for writing." << endl << flush;
            }

            for (i=0; i<sphres; i++)
            {
                if (reaches_spheroid[nodeno][i]->coordmtl)
                {
                    Bond** aabb = reaches_spheroid[nodeno][i]->get_rotatable_bonds();
                        if (aabb) for (l=0; aabb[l]; l++) aabb[l]->can_rotate = aabb[l]->can_flip = false;          // I SAID NO FLEXION OF ANY MCOORD RESIDUE AT ANY TIME.
                }
            }

            if (audit && !nodeno) fprintf(audit, "\nPose Candidate %d\nMovie offset: %d\n", pose, movie_offset);
            float predock_mclashes = protein->total_mclashes();

            for (i=0; i<nappears; i++)
            {
                if (!appears[i].disappear_iter || !appears[i].reappear_iter)
                {
                    appears[i].dispcen = protein->get_region_center(1, protein->get_end_resno()); // loneliest;
                    appears[i].update(protein, 0);
                }
            }

            // Stays rotation
            #if stays_rotation
            if (ligand->stay_close_mine && ligand->stay_close_other
                && fabs(ligand->stay_close_mine->is_polar()) >= hydrophilicity_cutoff
                && fabs(ligand->stay_close_other->is_polar()) >= hydrophilicity_cutoff
                )
            {
                bool do_stays_rot = true;
                if (ligand->stay_close2_mine && ligand->stay_close2_other)
                {
                    if (fabs(ligand->stay_close2_mine->is_polar()) >= hydrophilicity_cutoff
                        && fabs(ligand->stay_close2_other->is_polar()) >= hydrophilicity_cutoff)
                    {
                        do_stays_rot = false;
                    }
                }
                else if (pdpst == pst_best_binding)
                {
                    if (fabs(g_bbr->sec_res->hydrophilicity()) >= hydrophilicity_cutoff
                        && fabs(g_bbr->sec_tgt->polarity()) >= hydrophilicity_cutoff
                        )
                    {
                        do_stays_rot = false;
                    }
                }
                if (do_stays_rot)
                {
                    #if stays_rotation_verbose
                    cout << "Performing stays rotation step..." << endl << flush;
                    #endif
                    float frot = Search::stays_rotate_byinter(protein, ligand, ligand->stay_close_other, ligand->stay_close_mine);
                    frot += Search::stays_rotate_headtotail(protein, ligand, ligand->stay_close_other, ligand->stay_close_mine);
                    #if stays_rotation_verbose
                    cout << "Rotated " << (frot*fiftyseven) << "deg." << endl << endl << flush;
                    #endif
                }
            }
            #endif

            #if _dbg_atom_mov_to_clash
            movclash_prot = protein;
            movclash_cb = &check_moved_atom_for_clashes;
            #endif

            /////////////////////////////////////////////////////////////////////////////////
            // Main call to conformational search function.
            /////////////////////////////////////////////////////////////////////////////////
            Molecule::conform_molecules(cfmols, iters, &iteration_callback, progressbar ? &update_progressbar : nullptr, last_appear_iter?iters:0);
            if (end_program) poses = pose;
            float postdock_mclashes = protein->total_mclashes();
            float mclash_delta = postdock_mclashes - predock_mclashes;
            if (audit) fclose(audit);
            audit = nullptr;

            for (i=0; cfmol_known_good[i]; i++)
            {
                delete cfmol_known_good[i];
                cfmol_known_good[i] = 0;
            }

            #if accommodate_ligand_in_post
            // Adjustments to accommodate ligand
            for (j=0; j<1; j++)
            {
                for (i=0; cfmols[i]; i++)
                {
                    if (!cfmols[i]->is_residue()) continue;
                    float clash = cfmols[i]->get_intermol_clashes(ligand);
                    if (clash < clash_limit_per_aa) continue;
                    BallesterosWeinstein bw = protein->get_bw_from_resno(cfmols[i]->is_residue());
                    Point lc = ligand->get_barycenter();
                    if (bw.helix_no == 5 || bw.helix_no == 6)
                    {
                        // Find the terminus closest to the ligand
                        ResiduePlaceholder rpstart, rpend;
                        rpstart.bw = std::to_string(bw.helix_no) + (std::string)".s";
                        rpend.bw = std::to_string(bw.helix_no) + (std::string)".e";
                        rpstart.resolve_resno(protein);
                        rpend.resolve_resno(protein);
                        if (!rpstart.resno || !rpend.resno) continue;
                        AminoAcid *aastart = protein->get_residue(rpstart.resno), *aaend = protein->get_residue(rpend.resno);
                        if (!aastart || !aaend) continue;
                        bool start_is_closer = (aastart->get_CA_location().get_3d_distance(lc) < aaend->get_CA_location().get_3d_distance(lc));
                        AminoAcid *aaterminus = start_is_closer ? aastart : aaend;

                        // Find the side chain closest to the ligand
                        Atom* rna = protein->get_nearest_atom(lc, rpstart.resno, rpend.resno);
                        if (!rna) continue;
                        AminoAcid* aaneddamon = protein->get_residue(rna->residue);
                        if (!aaneddamon) continue;

                        // Find the internal contact closest to the ligand
                        int ici, rnsc = 0;
                        float rnscr;
                        for (ici=0; ici<n_soft_contact; ici++)
                        {
                            soft_contact_a[ici].resolve_resno(protein);
                            if (soft_contact_a[ici].resno >= rpstart.resno && soft_contact_a[ici].resno <= rpend.resno)
                            {
                                AminoAcid* aasc = protein->get_residue(soft_contact_a[ici].resno);
                                float r = aasc->get_CA_location().get_3d_distance(lc);
                                if (!rnsc || r<rnscr)
                                {
                                    rnsc = soft_contact_a[ici].resno;
                                    rnscr = r;
                                }
                            }
                            soft_contact_b[ici].resolve_resno(protein);
                            if (soft_contact_b[ici].resno >= rpstart.resno && soft_contact_b[ici].resno <= rpend.resno)
                            {
                                AminoAcid* aasc = protein->get_residue(soft_contact_b[ici].resno);
                                float r = aasc->get_CA_location().get_3d_distance(lc);
                                if (!rnsc || r<rnscr)
                                {
                                    rnsc = soft_contact_b[ici].resno;
                                    rnscr = r;
                                }
                            }
                        }
                        if (!rnsc)      // If there is no ic in this region, do a translation instead.
                        {
                            // TODO:
                            continue;
                        }

                        // If the sc is between the ic and the terminus, bend the region from ic to terminus to avoid clash
                        int AC = abs(rnsc - aaterminus->is_residue()),
                            AB = abs(rna->residue - aaterminus->is_residue()),
                            BC = abs(rna->residue - rnsc);
                        if (AC > AB && AC > BC)
                        {
                            ReshapeMotion rshpm;
                            rshpm.rap_start.bw = protein->get_bw_from_resno(rnsc).to_string();
                            rshpm.rap_end.bw = protein->get_bw_from_resno(aaterminus->get_residue_no()).to_string();
                            rshpm.rshpmt = rshpm_bend;
                            rshpm.fixclash = true;
                            rshpm.tgtligand = true;
                            rshpm.ligand = ligand;
                            rshpm.rap_index.bw = protein->get_bw_from_resno(aaneddamon->get_residue_no()).to_string();
                            rshpm.apply(protein);
                        }
                        // Otherwise, use the opposite terminus, but bend it back to its neighbor.
                        else
                        {
                            aaterminus = start_is_closer ? aaend : aastart;
                            if (aastart->get_residue_no() > aaend->get_residue_no()) throw 0xbadc0de;       // assertion
                            int sagitis, ardis = start_is_closer ? 1 : -1, nres = protein->get_end_resno();
                            AminoAcid* neighbor = nullptr;
                            if (!ardis) throw 0xbadc0de;
                            for (sagitis = aaterminus->get_residue_no()+ardis; sagitis > 0 && sagitis <= nres; sagitis += ardis)
                            {
                                neighbor = protein->get_residue(sagitis);
                                if (neighbor) break;
                            }
                            float rneighbs = neighbor->get_CA_location().get_3d_distance(aaterminus->get_CA_location());

                            ReshapeMotion rshpm;
                            rshpm.rap_start.bw = protein->get_bw_from_resno(rnsc).to_string();
                            rshpm.rap_end.bw = protein->get_bw_from_resno(aaterminus->get_residue_no()).to_string();
                            rshpm.rshpmt = rshpm_bend;
                            rshpm.fixclash = true;
                            rshpm.tgtligand = true;
                            rshpm.entire = true;
                            rshpm.morethan = false;
                            rshpm.ligand = ligand;
                            // rshpm.rap_index.bw = protein->get_bw_from_resno(aaneddamon->get_residue_no()).to_string();
                            rshpm.apply(protein);

                            int ciallon = 0;
                            for (sagitis = aaneddamon->get_residue_no()+ardis*3; sagitis > 0 && sagitis <= nres; sagitis += ardis)
                            {
                                AminoAcid* tmpaa = protein->get_residue(sagitis);
                                if (tmpaa)
                                {
                                    ciallon = sagitis;
                                    break;
                                }
                            }

                            if (false && 
                                ciallon >= min(rshpm.rap_start.resno, rshpm.rap_end.resno) 
                                && ciallon <= max(rshpm.rap_start.resno, rshpm.rap_end.resno))
                            {
                                rshpm.rap_start.bw = protein->get_bw_from_resno(ciallon).to_string();
                                rshpm.rap_end.bw = protein->get_bw_from_resno(aaterminus->get_residue_no()).to_string();
                                rshpm.rshpmt = rshpm_bend;
                                rshpm.fixclash = false;
                                rshpm.tgtligand = false;
                                rshpm.morethan = false;
                                rshpm.entire = false;
                                rshpm.rap_index.bw = protein->get_bw_from_resno(aaterminus->get_residue_no()).to_string();
                                rshpm.rap_target.bw = protein->get_bw_from_resno(neighbor->get_residue_no()).to_string();
                                rshpm.morethan = true;
                                rshpm.tgtdist = rneighbs;
                                rshpm.apply(protein);
                                rshpm.morethan = false;
                                rshpm.tgtdist = rneighbs + 2;
                                rshpm.apply(protein);
                            }
                        }
                    }
                    else
                    {
                        Atom* ra = cfmols[i]->get_nearest_atom(lc);
                        if (!ra) continue;
                        Vector mv = lc.subtract(ra->loc);
                        mv.r = 0.25;
                        ligand->move(mv);
                    }
                }

                if (pdpst == pst_best_binding && g_bbr[0].pri_res && g_bbr[0].pri_tgt)
                {
                    Atom* ra = g_bbr[0].pri_res->get_reach_atom(hbond);
                    if (ra)
                    {
                        g_bbr[0].pri_res->conform_atom_to_location(ra->name, g_bbr[0].pri_tgt->barycenter(), 15, 2.5);
                    }
                }
            }
            #endif

            #if optimize_internal_contacts_post_iterations
            if (nsoftrgn)
            {
                for (i=0; i<nsoftrgn; i++)
                {
                    int nc = softrgns[i].num_contacts();
                    for (j=0; j<nc; j++)
                    {
                        AminoAcid *aac1 = softrgns[i].get_local_contact(j, protein), *aac2 = softrgns[i].get_distant_contact(j, protein);
                        if (!aac1 || !aac2) continue;
                        Molecule* mols[3];
                        mols[0] = (Molecule*)aac1;
                        mols[1] = (Molecule*)aac2;
                        mols[2] = nullptr;
                        Molecule::conform_molecules(mols, contact_energy_optimization_iterations);
                    }
                }
            }
            #endif

            if (ligand->stay_close_water && ligand->check_stays_dry()) delete_water(ligand->stay_close_water);

            if (waters) for (i=0; waters[i]; i++)
            {
                float rpl = ligand->get_intermol_binding(waters[i]).clash;
                if (rpl > 50)
                {
                    for (j=i+1; waters[j]; j++) waters[j-1] = waters[j];
                    waters[j-1] = nullptr;
                    i--;
                }
            }

            if (!nodeno) // && outpdb.length())
            {
                protein->get_internal_clashes(1, protein->get_end_resno(), false);

                n = protein->get_end_resno();
                for (j=1; j <= n; j++)
                {
                    AminoAcid* aa = protein->get_residue(j);
                    if (aa)
                    {
                        tmp_pdb_residue[pose][j].copy_state(aa);
                        #if _dbg_residue_poses
                        cout << "tmp_pdb_residue[" << pose << "][" << j << "].copy_state(" << aa->get_name() << ")" << endl;
                        #endif
                    }
                }
                if (waters)
                {
                    for (j=0; waters[j]; j++)
                    {
                        tmp_pdb_waters[pose][j].copy_state(waters[j]);
                    }
                }
                n = nmtlcoords;
                for (j=0; j < n; j++)
                {
                    tmp_pdb_metal_locs[pose][j] = mtlcoords[j].mtl->loc;
                }
                tmp_pdb_ligand[pose].copy_state(ligand);
            }

            if (ligand->stay_close_mine && ligand->stay_close_mine->Z < 2)
            {
                Atom* heavy = ligand->stay_close_mine->get_heaviest_bonded_atom_that_isnt(nullptr);
                Bond* b = heavy->get_bond_by_idx(0);
                if (b) b = b->get_reversed();
                if (b && b->can_rotate && b->atom1)
                {
                    LocatedVector lv = b->get_axis();
                    lv.origin = b->atom1->loc;
                    float theta = find_angle_along_vector(ligand->stay_close_mine->loc, ligand->stay_close_other->loc,
                        lv.origin, lv);
                    float was = ligand->stay_close_mine->distance_to(ligand->stay_close_other);
                    b->rotate(theta);
                    float isnow = ligand->stay_close_mine->distance_to(ligand->stay_close_other);
                    if (isnow > was) b->rotate(theta*-2);
                }
            }

            // Add the current pose/path sequentially to the dr[][] array.
            // If the path node # is zero:
            // If it is the first (zeroth) entry, set the pose number to 1.
            // Otherwise, go through all of the preceding entries and:
            // Any entry with a smaller kJ/mol, increment its pose# but remember the smallest pre-increment pose #
            // from the lot of them;
            // Claim that new smallest pose# (which might be 1) as your own.
            #if _DBG_STEPBYSTEP
            if (debug) *debug << "Preparing output." << endl;
            #endif

            // Set the dock result properties and allocate the arrays.
            protein->set_conditional_basicities();
            #if _dbg_cond_basic
            AminoAcid* aadbg = protein->get_residue(155);
            cout << aadbg->get_name() << " charge = " << aadbg->get_charge() << endl;
            #endif

            protein->optimize_hydrogens();

            #if _dbg_atom_mov_to_clash
            movclash_prot = nullptr;
            movclash_cb = nullptr;
            #endif

            dr[drcount][nodeno] = DockResult(protein, ligand, search_size, addl_resno, drcount, waters);
            dr[drcount][nodeno].out_per_res_e = out_per_res_e;
            dr[drcount][nodeno].out_per_btyp_e = out_per_btyp_e;
            dr[drcount][nodeno].out_itemized_e_cutoff = out_itemized_e_cutoff;
            dr[drcount][nodeno].out_lig_int_e = out_lig_int_e;
            dr[drcount][nodeno].out_lig_pol_sat = out_lig_pol_sat;
            dr[drcount][nodeno].out_prox = out_prox;
            dr[drcount][nodeno].out_pro_clash = out_pro_clash;
            dr[drcount][nodeno].display_clash_atoms = display_clash_atoms;
            dr[drcount][nodeno].out_mc = out_mc;
            dr[drcount][nodeno].out_vdw_repuls = out_vdw_repuls;
            dr[drcount][nodeno].mbbr = &g_bbr[0];
            dr[drcount][nodeno].estimated_TDeltaS = g_bbr[0].estimate_DeltaS() * temperature;

            if (dr[drcount][nodeno].ligand_pocket_occlusion < 0.65)
            {
                dr[drcount][nodeno].disqualified = true;
                std::string reason = "Insufficient occlusion ";
                reason += std::to_string(dr[drcount][nodeno].ligand_pocket_occlusion);
                reason += (std::string)". ";
                dr[drcount][nodeno].disqualify_reason += reason;
            }

            n = protein->get_end_resno();
            for (i=1; i<sphres; i++)
            {
                AminoAcid* aa1 = reaches_spheroid[nodeno][i];
                if (!aa1) continue;
                for (j=aa1->get_residue_no()+2; j<=n; j++)
                {
                    AminoAcid* aa2 = protein->get_residue(j);
                    if (!aa2) continue;

                    float f = aa1->get_intermol_clashes(aa2);
                    if (f > clash_limit_per_aa*10)
                    {
                        dr[drcount][nodeno].disqualified = true;
                        std::string reason = (std::string)"Side chain clash " + std::to_string(f) + (std::string)". ";
                        dr[drcount][nodeno].disqualify_reason += reason;
                        i=j=n+2;
                        break;
                    }
                }
            }

            if (nsoftrgn)
            {
                float anomaly = 0, sanomaly = 0;
                for (i=0; i<nsoftrgn; i++)
                {
                    if (!softrgns[i].check_chain_constraints(protein))
                    {
                        dr[drcount][nodeno].disqualified = true;
                        std::string reason = "Soft chain constraint violation: ";
                        if (softrgns[i].prgv) reason += std::to_string(softrgns[i].prge) + (std::string)"~"
                            + std::to_string(softrgns[i].rgn.start) + (std::string)"="
                            + std::to_string(softrgns[i].prgd) + "A. ";
                        if (softrgns[i].nrgv) reason += std::to_string(softrgns[i].rgn.end) + (std::string)"~"
                            + std::to_string(softrgns[i].nrgs) + (std::string)"="
                            + std::to_string(softrgns[i].prgd) + "A. ";
                        dr[drcount][nodeno].disqualify_reason += reason;
                    }
                    if (!softrgns[i].num_contacts())
                    {
                        #if _dbg_soft
                        cout << "No contacts for soft region " << i << endl;
                        #endif
                        continue;
                    }
                    // softrgns[i].optimize_contact(protein, i);
                    float lanom = softrgns[i].contact_anomaly(protein);
                    #if _dbg_soft
                    // TODO:
                    /* float r = softrgns[i].pivaa1->get_CA_location().get_3d_distance(softrgns[i].pivaa2->get_CA_location());
                    cout << softrgns[i].pivaa1->get_letter() << softrgns[i].pivaa1->get_residue_no() << "..."
                        << softrgns[i].pivaa2->get_letter() << softrgns[i].pivaa2->get_residue_no()
                        << " @ " << r << ", expected " << softrgns[i].pivaar << ", anomaly " << lanom << endl; */
                    #endif
                    anomaly += lanom;
                    sanomaly += softrgns[i].contact_distance_anomaly(protein);

                    int nc = softrgns[i].num_contacts();
                    for (j=0; j<nc; j++)
                    {
                        if (softrgns[i].is_contact_paired(j)) continue;
                        float canom = softrgns[i].contact_anomaly(protein, j);
                        if (fabs(canom) < 0.1) continue;
                        AminoAcid *aac1 = softrgns[i].get_local_contact(j, protein), *aac2 = softrgns[i].get_distant_contact(j, protein);
                        if (aac1 && aac2) dr[drcount][nodeno].miscdata += (std::string)"Contact anomaly for "
                            + (std::string)"region " + std::to_string(i) + (std::string)" "
                            + (std::string)aac1->get_name() + (std::string)"..." + (std::string)aac2->get_name()
                            + (std::string)": " + std::to_string(canom)
                            + (std::string)"\n";
                    }
                }
                dr[drcount][nodeno].miscdata += (std::string)"Raw ligand binding energy: " 
                    + std::to_string(dr[drcount][nodeno].kJmol) + (std::string)" kJ/mol.\n";
                dr[drcount][nodeno].miscdata += (std::string)"Soft contact anomaly: " + std::to_string(anomaly) + (std::string)" kJ/mol.\n";
                dr[drcount][nodeno].miscdata += (std::string)"Soft spatial anomaly: " + std::to_string(sanomaly) + (std::string)" A.\n";

                if (dr[drcount][nodeno].kJmol > 0 && !output_something_even_if_it_is_wrong)
                {
                    dr[drcount][nodeno].disqualified = true;
                    dr[drcount][nodeno].disqualify_reason += (std::string)"Unfavorable ligand binding energy. ";
                }
                else if ((anomaly + dr[drcount][nodeno].kJmol) > 0)
                {
                    dr[drcount][nodeno].disqualified = true;
                    dr[drcount][nodeno].disqualify_reason += (std::string)"Soft anomaly greater than ligand binding energy. ";
                }
                else dr[drcount][nodeno].kJmol += anomaly;
            }

            for (i=0; cfmols[i]; i++)
            {
                if (cfmols[i]->is_residue())
                {
                    AminoAcid* aacfmolsi = (AminoAcid*)cfmols[i];
                    if (aacfmolsi->get_letter() == 'P') continue;
                }
                if (cfmols[i]->get_internal_clashes() > clash_limit_per_aa*10)
                {
                    dr[drcount][nodeno].disqualified = true;
                    dr[drcount][nodeno].disqualify_reason += (std::string)cfmols[i]->get_name() + (std::string)" internal clashes too great. ";
                }
            }

            float btot = dr[drcount][nodeno].kJmol;
            float pstot = dr[drcount][nodeno].polsat;
            if (isomers.size()) dr[drcount][nodeno].isomer = ligand->get_name();

            if ((pose==1 && !nodeno) || best_energy > btot) best_energy = btot;
            if ((pose==1 && !nodeno) || best_worst_clash > dr[drcount][nodeno].worst_energy) best_worst_clash = dr[drcount][nodeno].worst_energy;

            #if reuse_best_pose
            if (!aa_best_pose || best_pose_energy > btot)
            {
                float cp = 0;
                for (i=0; i<ncvtys; i++)
                {
                    if (!cvtys[i].count_partials()) continue;
                    cp += cvtys[i].molecule_inside_pocket(ligand);
                }

                if (!ncvtys || cp >= min_cvty_ctnmt/2)
                {
                    n = protein->get_end_resno();
                    l = max(1, (int)isomers.size());
                    if (!aa_best_pose) aa_best_pose = new Pose[n*l+4];
                    for (i=1; i<n; i++)
                    {
                        AminoAcid* aa = protein->get_residue(i);
                        if (!aa) continue;
                        aa_best_pose[i+n*isono].copy_state(aa);
                    }
                    best_pose_energy = btot;
                    ligand_best_pose[isono].copy_state(ligand);
                    best_ligand_stays = ligand->stay_close_mine;
                    best_ligand_stays_other = ligand->stay_close_other;
                }
            }
            #endif

            #if compute_clashdirs
            n = protein->get_end_resno();
            for (i=1; i<=n; i++)
            {
                if (dr[drcount][nodeno].residue_clash[i])
                {
                    for (j=0; regions[j].start; j++)
                    {
                        if (i >= regions[j].start && i <= regions[j].end)
                        {
                            int rgcen, hxno = atoi(&regions[j].name[3]);
                            if (hxno >= 1 && hxno <= 7) rgcen = protein->get_bw50(hxno);
                            else rgcen = (regions[j].start + regions[j].end)/2;
                            k = abs(i - rgcen);
                            if (k < 5) k = 1;
                            else k = (i < rgcen) ? 0 : 2;
                            region_clashes[j][k] = region_clashes[j][k].add(dr[drcount][nodeno].res_clash_dir[i]);
                        }
                    }
                }
            }
            #endif

            dr[drcount][nodeno].proximity = ligand->get_barycenter().get_3d_distance(nodecen);

            if (pdpst == pst_best_binding && out_bb_pairs)
            {
                dr[drcount][nodeno].miscdata += (std::string)"Best-Binding Pairs:\n";
                for (l=0; g_bbr[l].pri_res; l++)
                {
                    if (g_bbr[l].pri_tgt && g_bbr[l].pri_res)
                    {
                        dr[drcount][nodeno].miscdata += (std::string)g_bbr[l].pri_tgt->to_std_string() + (std::string)"..." 
                            + (std::string)g_bbr[l].pri_res->get_name() + (std::string)"\n";
                    }
                    if (g_bbr[l].sec_tgt && g_bbr[l].sec_res)
                    {
                        dr[drcount][nodeno].miscdata += (std::string)g_bbr[l].sec_tgt->to_std_string() + (std::string)"..." 
                            + (std::string)g_bbr[l].sec_res->get_name() + (std::string)"\n";
                    }
                    if (g_bbr[l].tert_tgt && g_bbr[l].tert_res)
                    {
                        dr[drcount][nodeno].miscdata += (std::string)g_bbr[l].tert_tgt->to_std_string() + (std::string)"..." 
                            + (std::string)g_bbr[l].tert_res->get_name() + (std::string)"\n";
                    }
                }
                dr[drcount][nodeno].miscdata += (std::string)"\n";
            }

            #if _DBG_STEPBYSTEP
            if (debug) *debug << "Allocated memory." << endl;
            #endif

            dr[drcount][nodeno].auth = pose;

            if (!nodeno)
            {
                if ((dr[drcount][nodeno].ligand_self + ligand->total_eclipses()) < -clash_limit_per_aa*2)
                {
                    #if _dbg_worst_energy
                    cout << "Internal ligand energy " << -dr[drcount][nodeno].ligand_self << " out of range." << endl << endl;
                    #endif

                    break;          // Exit nodeno loop.
                }
                // else cout << "Internal ligand energy " << -dr[drcount][nodeno].ligand_self << " satisfactory." << endl << endl;

                #if !_dbg_allow_excessive_aa_clashes
                if (dr[drcount][nodeno].worst_energy > l_atom_clash_limit || dr[drcount][nodeno].worst_nrg_aa > clash_limit_per_aa)
                {
                    #if _dbg_worst_energy
                    cout << "Total binding energy " << dr[drcount][nodeno].kJmol
                        << " and worst energy " << dr[drcount][nodeno].worst_energy;
                    if (dr[drcount][nodeno].worst_clash_1 && dr[drcount][nodeno].worst_clash_2)
                    {
                        cout << " (" << dr[drcount][nodeno].worst_clash_1->residue << ":" << dr[drcount][nodeno].worst_clash_1->name
                            << "-" << dr[drcount][nodeno].worst_clash_2->residue << ":" << dr[drcount][nodeno].worst_clash_2->name
                            << ") ";
                    }
                    cout << "; skipping." << endl << endl;
                    #endif

                    // cout << "Least favorable binding energy " << dr[drcount][nodeno].worst_energy << " out of range." << endl << endl;
                    break;          // Exit nodeno loop.
                }
                // else cout << "Least favorable binding energy " << dr[drcount][nodeno].worst_energy << " satisfactory." << endl << endl;
                #endif

                if (pose==1) dr[drcount][nodeno].pose = pose;
                else
                {
                    int bestpose = pose;
                    for (i=0; i<drcount; i++)
                    {
                        if ((	(dr[i][0].kJmol - dr[i][0].ikJmol - dr[i][0].polsat * polar_sat_influence_for_scoring)
                                >
                                (dr[drcount][nodeno].kJmol - dr[drcount][nodeno].ikJmol - dr[drcount][nodeno].polsat * polar_sat_influence_for_scoring)
                            )
                            ||
                            (	(dr[i][0].kJmol - dr[i][0].polsat * polar_sat_influence_for_scoring)
                                >
                                (btot - pstot * polar_sat_influence_for_scoring)
                            ))
                        {
                            if (dr[i][0].pose < bestpose || bestpose < 0) bestpose = dr[i][0].pose;
                            dr[i][0].pose++;
                        }
                    }
                    dr[drcount][nodeno].pose = bestpose;
                    // cout << "Around the posie: "; for (i=0; i<=drcount; i++) cout << dr[i][nodeno].pose << " "; cout << endl << "Best: " << bestpose << endl;
                }
                #if _DBG_STEPBYSTEP
                if (debug) *debug << "Added pose to output array." << endl;
                #endif

                #if _DBG_MAX_CLASHES
                cout << "Pose " << dr[drcount][nodeno].pose << " maxclash " << maxclash << " kJmol " << dr[drcount][nodeno].kJmol << endl;
                #endif
            }

            if (btot <= kJmol_cutoff && !dr[drcount][0].disqualified) success_sofar = true;

            // if (dr[drcount][0].disqualified) cout << dr[drcount][nodeno].disqualify_reason << endl << endl;

            // For performance reasons, once a path node (including #0) fails to meet the binding energy threshold, discontinue further
            // calculations for this pose.
            if (btot > kJmol_cutoff)
            {
                #if _dbg_worst_energy
                cout << "Total binding energy " << -btot << "; skipping." << endl << endl;
                #endif
                drcount++;
                break;
            }
            else if (nodeno == pathnodes) drcount++;
        }	// nodeno loop.
    } // pose loop.

    /////////////////////////////////////////////////////////////////////////////////
    // End main loop.
    /////////////////////////////////////////////////////////////////////////////////

    #if _DBG_STEPBYSTEP
    if (debug) *debug << "Finished poses." << endl;
    #endif

    if (progressbar)
    {
        progb.erase();
    }

    // Output the dr[][] array in order of increasing pose number.
    cout << endl;
    if (output) *output << endl;

    const float energy_mult = kcal ? _kcal_per_kJ : 1;
    pose = 1;
    std::string auths;
    for (i=1; i<=poses; i++)
    {
        for (j=0; j<poses; j++)
        {
            protein = &pose_proteins[j];
            ligand = &pose_ligands[j+1];

            if (dr[j][0].disqualified)
            {
                if (outdisq)
                {
                    cout << "Disqualified pose: " << dr[j][0].disqualify_reason << endl << dr[j][0].miscdata << endl;
                    if (output) *output << "Disqualified pose: " << dr[j][0].disqualify_reason << endl << dr[j][0].miscdata << endl;
                }

                if (!strstr(outdisqrs.c_str(), dr[j][0].disqualify_reason.c_str()))
                    outdisqrs += dr[j][0].disqualify_reason + (std::string)" ";

                continue;
            }
            dr[j][0].ligpos.restore_state(ligand);

            if (ncvtys)
            {
                dr[j][0].disqualified = true;
                int cno;
                float cp = 0;
                for (cno = 0; cno < ncvtys; cno++)
                {
                    if (!cvtys[cno].count_partials()) continue;
                    cp += cvtys[cno].molecule_inside_pocket(ligand);
                }

                if (cp >= min_cvty_ctnmt) dr[j][0].disqualified = false;
                if (dr[j][0].disqualified) dr[j][0].disqualify_reason += (std::string)"Cavity containment. ";
            }

            if (dr[j][0].disqualified)
            {
                if (outdisq)
                {
                    cout << "Disqualified pose: " << dr[j][0].disqualify_reason << endl;
                    if (output) *output << "Disqualified pose: " << dr[j][0].disqualify_reason << endl;
                }

                if (!strstr(outdisqrs.c_str(), dr[j][0].disqualify_reason.c_str()))
                    outdisqrs += dr[j][0].disqualify_reason + (std::string)" ";

                continue;
            }

            if (dr[j][0].pose == i && dr[j][0].pdbdat.length())
            {
                if (dr[j][0].kJmol <= kJmol_cutoff || (output_something_even_if_it_is_wrong && !j))
                {
                    if (dr[j][0].proximity > search_size.magnitude()) continue;
                    if (dr[j][0].worst_nrg_aa > clash_limit_per_aa) continue;
                    protein = &pose_proteins[j];

                    auths += (std::string)" " + std::to_string(dr[j][0].auth);

                    if (!best_acc_energy) best_acc_energy = dr[j][0].kJmol;

                    for (k=0; k<=pathnodes; k++)
                    {
                        // If pathnode is not within kJ/mol cutoff, abandon it and all subsequent pathnodes of the same pose.
                        if (dr[j][k].kJmol > kJmol_cutoff && !output_something_even_if_it_is_wrong)
                        {
                            cout << "Pose " << pose << " node " << k
                                 << " energy " << dr[j][k].kJmol*energy_mult
                                 << " is outside of limit; aborting path nodes." << endl;
                            if (output) *output << "Pose " << pose << " node " << k
                                                << " energy " << -dr[j][k].kJmol*energy_mult
                                                << " is outside of limit; aborting path nodes." << endl;
                            break;
                        }

                        /*if (flex && !dr[j][k].pdbdat.length())
                        {
                            cout << "Pose " << j << " node " << k << " is missing." << endl;
                            if (output) *output << "Pose " << j << " node " << k << " is missing." << endl;
                            continue;
                        }*/

                        if (!k)
                        {
                            if (auditfn.length()) audit = fopen(auditfn.c_str(), "ab");
                            if (audit)
                            {
                                fprintf(audit, "Pose candidate %d accepted as output pose %d.\n", j+1, dr[j][0].pose);
                                fclose(audit);
                            }
                        }
                        do_pose_output(&dr[j][k], k, energy_mult, tmp_pdb_waters[pose], tmp_pdb_metal_locs[pose]);

                        if (!k) found_poses++;
                    }
                    pose++;
                }
                else
                {
                    if (i == 1)
                    {
                        if (kcal)
                        {
                            cout << "No poses found within kcal/mol limit." << endl;
                            if (output) *output << "No poses found within kcal/mol limit." << endl;
                        }
                        else
                        {
                            cout << "No poses found within kJ/mol limit." << endl;
                            if (output) *output << "No poses found within kJ/mol limit." << endl;
                        }
                    }
                    cout << "Exiting." << endl;
                    goto _exitposes;
                }

                break;
            }

            _next_pose:
            ;
        }
    }

_exitposes:
    if (found_poses < poses && triesleft)
    {
        triesleft--;
        goto _try_again;
    }

    if (!found_poses && output_something_even_if_it_is_wrong)
    {
        float leastbad = Avogadro;
        int lbi = 0;
        for (l=0; l<poses; l++)
        {
            if (dr[l][0].kJmol && dr[l][0].kJmol < leastbad)
            {
                leastbad = Avogadro;
                lbi = l;
            }
        }
        pose = 1;
        protein = &pose_proteins[lbi];
        do_pose_output(&dr[lbi][0], 0, energy_mult, tmp_pdb_waters[pose], tmp_pdb_metal_locs[pose]);
    }
    else 
    {
        cout << found_poses << " pose(s) found." << endl;
        if (output) *output << found_poses << " pose(s) found." << endl;
        if (debug) *debug << found_poses << " pose(s) found." << endl;
    }

    cout << "Best candidate pose energy: " << (kcal ? best_energy/_kcal_per_kJ : best_energy) << (kcal ? " kcal/mol." : " kJ/mol.") << endl;
    if (output) *output << "Best candidate pose energy: " << (kcal ? best_energy/_kcal_per_kJ : best_energy) << (kcal ? " kcal/mol." : " kJ/mol.") << endl;
    if (debug) *debug << "Best candidate pose energy: " << (kcal ? best_energy/_kcal_per_kJ : best_energy) << (kcal ? " kcal/mol." : " kJ/mol.") << endl;

    if (found_poses)
    {
        cout << "Best accepted pose energy: " << (kcal ? best_acc_energy/_kcal_per_kJ : best_acc_energy) << (kcal ? " kcal/mol." : " kJ/mol.") << endl;
        if (output) *output << "Best accepted pose energy: " << (kcal ? best_acc_energy/_kcal_per_kJ : best_acc_energy) << (kcal ? " kcal/mol." : " kJ/mol.") << endl;
        if (debug) *debug << "Best accepted pose energy: " << (kcal ? best_acc_energy/_kcal_per_kJ : best_acc_energy) << (kcal ? " kcal/mol." : " kJ/mol.") << endl;
    }

    cout << "Best worst clash: " << (kcal ? best_worst_clash/_kcal_per_kJ : best_worst_clash) << (kcal ? " kcal/mol." : " kJ/mol.") << endl;
    if (output) *output << "Best worst clash: " << (kcal ? best_worst_clash/_kcal_per_kJ : best_worst_clash) << (kcal ? " kcal/mol." : " kJ/mol.") << endl;
    if (debug) *debug << "Best worst clash: " << (kcal ? best_worst_clash/_kcal_per_kJ : best_worst_clash) << (kcal ? " kcal/mol." : " kJ/mol.") << endl;

    if (!found_poses)
    {
        if (outdisqrs.size())
        {
            cout << "Reasons for candidate pose disqualification: " << outdisqrs << endl;
            if (output) *output << "Reasons for candidate pose disqualification: " << outdisqrs << endl;
            if (debug) *debug << "Reasons for candidate pose disqualification: " << outdisqrs << endl;
        }
    }

    #if compute_clashdirs
    if (regions)
    {
        cout << endl;
        if (output) *output << endl;
        for (i=0; regions[i].start; i++)
        {
            for (j=0; j<3; j++) if (region_clashes[i][j].r)
            {
                cout << regions[i].name;
                if (output) *output << regions[i].name;
                if (!j)
                {
                    cout << ".nseg";
                    if (output) *output << ".nseg";
                }
                else if (j==1)
                {
                    cout << ".center";
                    if (output) *output << ".center";
                }
                else if (j==2)
                {
                    cout << ".cseg";
                    if (output) *output << ".cseg";
                }
                cout << ".clashdir = " << (Point)region_clashes[i][j] << endl;
                if (output) *output << ".clashdir = " << (Point)region_clashes[i][j] << endl;
            }
        }
        cout << endl;
        if (output) *output << endl;
    }
    #endif

    if (met) delete met;

    time_t finished = time(NULL);
    int seconds = finished-began;
    int minutes = seconds/60;
    seconds -= 60*minutes;
    int hours = minutes/60;
    minutes -= 60*hours;

    std::string elapsed;
    if (hours)
    {
        elapsed += std::to_string(hours);
        elapsed += (std::string)":";
        if (minutes < 10) elapsed += (std::string)"0";
    }
    elapsed += std::to_string(minutes);
    elapsed += (std::string)":";
    if (seconds < 10) elapsed += (std::string)"0";
    elapsed += std::to_string(seconds);

    cout << "\nCalculation time: " << elapsed << "." << endl;
    if (output) *output << "\nCalculation time: " << elapsed << "." << endl;
    if (debug) *debug << "\nCalculation time: " << elapsed << "." << endl;

    if (output) output->close();
    if (append_pdb)
    {
        if (output)
        {
            pf = fopen(temp_pdb_file.length() ? temp_pdb_file.c_str() : protfname, "r");
            if (!pf)
            {
                cerr << "Error trying to read " << protfname << endl;
                return 0xbadf12e;
            }
            protein->load_pdb(pf);
            fclose(pf);
            apply_protein_specific_settings(protein);
            FILE* pf = fopen(outfname, "ab");
            fprintf(pf, "\nOriginal PDB:\n");
            protein->save_pdb(pf);
            fclose(pf);
            cout << "PDB appended to output file." << endl;
        }
        else cerr << "ERROR: Append PDB can only be used when specifying an output file." << endl;
    }
    if (output_each_iter)
    {
        pf = fopen(temp_pdb_file.length() ? temp_pdb_file.c_str() : protfname, "r");
        if (!pf)
        {
            cerr << "Error trying to read " << protfname << endl;
            return 0xbadf12e;
        }
        protein->load_pdb(pf);
        fclose(pf);
        apply_protein_specific_settings(protein);

        pf = fopen(itersfname.c_str(), "ab");
        if (pf)
        {
            fprintf(pf, "\nOriginal PDB:\n");
            protein->save_pdb(pf);
            fclose(pf);
        }
    }

    if (temp_pdb_file.length()) std::remove(temp_pdb_file.c_str());

    if (debug) debug->close();

    return 0;
}


















