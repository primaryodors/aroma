
#ifndef _RESHAPE
#define _RESHAPE

#include "protein.h"

#define MAX_ICHX 10
#define MAX_HXIC 10

enum ReshapeType
{
    rshp_GPCR,
};

enum ReshapeMotionType
{
    rshpm_xlate,
    rshpm_pivot,
    rshpm_twist,
    rshpm_wind,
    rshpm_prox,
    rshpm_delete,
    rshpm_flex,
    rshpm_bend,
};

class ReshapeMotion
{
    public:
    void apply(Protein* p);
    void do_measurement(Protein* p);

    ReshapeMotionType rshpmt;
    ResidueAtomPlaceholder rap_start, rap_end, rap_fulcrum, rap_index, rap_target;
    char ba1[16] = {0}, ba2[16] = {0};
    float tgtdist = 0;
    bool fixclash = false;
    bool morethan = false;
    bool entire = false;
    bool tgtligand = false;
    bool measure = false;
    Molecule* ligand = nullptr;

    protected:
    bool get_pt_index_and_tgt(Protein* p, Point* index, Point* target, Atom** a_index = nullptr, Atom** a_target = nullptr);
    float measure_index_tgt_clashes(Protein* p);
    bool fix_clash(Protein* p, int sr, int er, Point pt_fulcrum, int iters=60, float step = 1.5*fiftyseventh);
    bool set_distance(Protein* p, int sr, int er, Point pt_fulcrum, Point pt_index, Point pt_target, int moreorless = -1, float amount=1.0);
    Atom* closest_to_ligand;
};

class Reshape
{
    public:
    void load_rshpm_file(ReshapeType rshpt, Molecule* ligand);
    void load_rshpm_file(const char* infname, Molecule* ligand);
    void apply(Protein* p, bool ones_with_ligands = false);

    protected:
    ReshapeMotion* m_rshpm = nullptr;
    int nrshpm = 0;
};

class ICHelixGroup;
class InternalContact
{
    friend class ICHelixGroup;

    public:
    ResiduePlaceholder res1;
    ResiduePlaceholder res2;
    float r_optimal = 2.5;
    float tolerance = 0.5;

    Vector atom_distance(Protein* prot);

    protected:
    AminoAcid *aa1 = nullptr, *aa2 = nullptr;
    Atom *a1 = nullptr, *a2 = nullptr;
};

class ICHelix
{
    public:
    int hxno = 0;
    ResiduePlaceholder start;
    ResiduePlaceholder end;
    InternalContact ic[MAX_HXIC];
    int n_ic = 0;
    bool hxstatic = false;

    bool contains(InternalContact* ic);
};

class ICHelixGroup
{
    public:
    ICHelix helices[MAX_ICHX];
    int n_helix = 0;
    ResiduePlaceholder deletion_start[MAX_ICHX], deletion_end[MAX_ICHX];
    int n_deletion = 0;

    LocRotation get_motion(InternalContact* ic, ICHelix* ich = nullptr, Protein* prot = nullptr);
    void load_ic_file(const char* filename);
    float optimize_helices(Protein* prot, int iterations = 20);
    float contact_anomaly(Protein* prot);

    protected:
    Protein* m_prot = nullptr;
};

extern bool rshp_verbose;

#endif