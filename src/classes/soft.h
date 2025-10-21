
#include "protein.h"

#ifndef _SOFT
#define _SOFT

#define soft_contact_alloc_block 8

class SoftContact
{
    public:
    int local = 0, distant = 0;
    bool paired = false;
    float energy = 0;
    float CA_distance = 0;
};

class SoftRegion
{
    protected:
    SoftContact* contacts = nullptr;
    int prev_rgn_end = 0;
    int next_rgn_start = 0;
    int allocated = 0;
    float prev_rgn_dist = 0, next_rgn_dist = 0;
    bool prev_rgn_violated = false, next_rgn_violated = false;

    public:
    Region rgn;
    float initclash = 0;

    int num_contacts();
    AminoAcid* get_local_contact(int i, Protein* p);
    AminoAcid* get_distant_contact(int i, Protein* p);
    float get_contact_original_distance(int i);
    bool is_contact_paired(int i) { return contacts[i].paired; }
    Atom* get_pivot_atom_by_contact_idx(int i, Protein* p);
    float contact_anomaly(Protein *p, int i = -1, bool ignore_paired = true);
    float contact_distance_anomaly(Protein *p, int i = -1, bool ignore_paired = true);
    void add_contact(int local, int distant, Protein* p, bool paired = false);
    void link_region(SoftRegion* prev);
    bool check_chain_constraints(Protein* prot);
    void optimize_contact(Protein *p, int i);
    SoftRegion();
    ~SoftRegion();
    SoftRegion(const SoftRegion& sr);
    SoftRegion& operator=(const SoftRegion& sr);
    SoftRegion(SoftRegion&& sr) noexcept;
    SoftRegion& operator=(SoftRegion&& sr) noexcept;

    const int& prge = prev_rgn_end;
    const int& nrgs = next_rgn_start;
    const float& prgd = prev_rgn_dist;
    const float& nrgd = next_rgn_dist;
    const bool& prgv = prev_rgn_violated;
    const bool& nrgv = next_rgn_violated;
};

void soft_docking_iteration(Protein* protein, Molecule* ligand, int nsoftrgn, SoftRegion* softrgns, float softness);

extern float soft_contact_elasticity;

#endif