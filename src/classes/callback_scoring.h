#ifndef AROMA_CALLBACK_SCORING_H
#define AROMA_CALLBACK_SCORING_H

#include "molecule.h"

template<class Potential>
bool callback_flexion_unfavorable(const Interaction& score, float missed, Potential potential)
{
    float energy = score.clash + score.repulsive - score.attractive;
    if (score.clash > clash_limit_per_aa || (score.clash > 0 && energy > 5) || missed >= 1) return true;
    return energy > -0.3 * potential();
}

#endif
