
#include <memory>
#include <algorithm>
#include "protein.h"
#include "cavity.h"

#ifndef _SEARCH
#define _SEARCH

class LigandTarget
{
    public:
    Atom* single_atom = nullptr;
    Conjugation* conjgrp = nullptr;

    float charge();
    float polarity();
    Point barycenter();
    bool is_pi();
    bool has_hb_acceptors();
    bool has_hb_donors();
    float importance(Molecule* mol);
    bool contains(LigandTarget* lt);
    bool contains(Atom* a);
    intera_type best_interaction();
    std::string to_std_string();
};

class BestBindingResult
{
    public:
    float score(Point ligand_center, Cavity* container = nullptr);
    void add_to_candidates();
    int num_assigned();
    bool is_equivalent(BestBindingResult* bbr2);
    float estimate_DeltaS();
    Point barycenter();

    Protein* protein = nullptr;
    Molecule* ligand = nullptr;
    AminoAcid* pri_res = nullptr;
    LigandTarget* pri_tgt = nullptr;
    AminoAcid* sec_res = nullptr;
    LigandTarget* sec_tgt = nullptr;
    AminoAcid* tert_res = nullptr;
    LigandTarget* tert_tgt = nullptr;
    float probability = 0;
    float cached_score = 0;

    protected:
    float entropic_score = 0;
};

class Search
{
    public:
    static void do_tumble_spheres(Protein* protein, Molecule* ligand, Point l_pocket_cen);
    static void copy_ligand_position_from_file(Protein* protein, Molecule* ligand, const char* filename, const char* ligname, int auth_resno);

    static int identify_ligand_pairing_targets(Molecule* ligand, LigandTarget* results, int max_results);
    static void pair_targets(Molecule* ligand, LigandTarget* targets, AminoAcid** pocketres, Point loneliest, 
        BestBindingResult* output, Cavity* container = nullptr, bool allow_thiolation = true);
    static void scan_protein(Protein* prot, Molecule* ligand, LigandTarget* targets, BestBindingResult* results, int max_results, Box limit = Box());
    static void align_targets(Molecule* ligand, Point pocketcen, BestBindingResult* bbr, float amt = 1);
    static bool target_compatibility(float chg1, float chg2, float pol1, float pol2, int pi1, int pi2,
        bool hba1, bool hba2, bool hbd1, bool hbd2);
    static bool target_compatibility(AminoAcid* aa, LigandTarget* lt);
    static void clear_candidates();

    protected:
    static bool any_resnos_priority;
};

extern Point search_size, loneliest;
extern std::vector<int> exclusion;

#define MAX_CS_RES 4096
extern int agqty;
extern AminoAcid* cs_res[MAX_CS_RES];
extern intera_type cs_bt[MAX_CS_RES];
extern int cs_res_qty;
extern int cs_idx;

extern AminoAcid *g_pri_pairres, *g_sec_pairres, *g_tert_pairres;
extern LigandTarget *g_pri_ltarg, *g_sec_ltarg, *g_tert_ltarg;

#define MAX_BBR_CANDIDATES 65536
extern BestBindingResult bbr_candidates[MAX_BBR_CANDIDATES];

std::ostream& operator<<(std::ostream& os, const LigandTarget& lt);
std::ostream& operator<<(std::ostream& os, const BestBindingResult& bbr);

#endif
