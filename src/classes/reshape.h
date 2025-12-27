
#ifndef _RESHAPE
#define _RESHAPE

#include "protein.h"

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

    ReshapeMotionType rshpmt;
    ResidueAtomPlaceholder rap_start, rap_end, rap_fulcrum, rap_index, rap_target;
    char ba1[16] = {0}, ba2[16] = {0};
    float tgtdist = 0;
    bool fixclash = false;
    bool morethan = false;
    bool entire = false;
    bool tgtligand = false;
    Molecule* ligand = nullptr;

    protected:
    bool get_pt_index_and_tgt(Protein* p, Point* index, Point* target);
    bool fix_clash(Protein* p, int sr, int er, Point pt_fulcrum, int iters=60);
    bool set_distance(Protein* p, int sr, int er, Point pt_fulcrum, Point pt_index, Point pt_target, int moreorless = -1, float amount=1.0);
    Atom* closest_to_ligand;
};

class Reshape
{
    public:
    void load_rshpm_file(ReshapeType rshpt, Molecule* ligand);
    void apply(Protein* p, bool ones_with_ligands = false);

    protected:
    ReshapeMotion* m_rshpm = nullptr;
    int nrshpm = 0;
};

#endif