
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