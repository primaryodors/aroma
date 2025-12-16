
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
};

class ActiveMotion
{
    public:
    void apply(Protein* p);

    AcvMotionType acvmt;
    ResidueAtomPlaceholder rap_start, rap_end, rap_fulcrum, rap_index, rap_target;
    float tgtdist = 0;
};

class Activation
{
    public:
    void load_acvm_file(AcvType acvt);
    void apply(Protein* p);

    protected:
    ActiveMotion* m_acvm = nullptr;
    int nacvm = 0;
};

#endif