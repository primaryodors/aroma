
#ifndef _ACTIVATION
#define _ACTIVATION

#include "protein.h"

enum AcvType
{
    acv_GPCR,
};

enum AcvMotionType
{
    acvm_xlate,
    acvm_pivot,
    acvm_twist,
    acvm_wind,
    acvm_prox,
    acvm_delete,
    acvm_flex,
    acvm_bend,
};

class ActiveMotion
{
    public:
    void apply(Protein* p);

    AcvMotionType acvmt;
    ResidueAtomPlaceholder rap_start, rap_end, rap_fulcrum, rap_index, rap_target;
    char ba1[16] = {0}, ba2[16] = {0};
    float tgtdist = 0;
    bool fixclash = false;
    bool morethan = false;
    bool entire = false;
    bool tgtligand = false;
    Molecule* ligand = nullptr;
};

class Activation
{
    public:
    void load_acvm_file(AcvType acvt, Molecule* ligand);
    void apply(Protein* p, bool ones_with_ligands = false);

    protected:
    ActiveMotion* m_acvm = nullptr;
    int nacvm = 0;
};

#endif