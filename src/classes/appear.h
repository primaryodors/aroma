
#include "protein.h"

#ifndef _APPEAR
#define _APPEAR

struct Appear               // Should rename this to something better.
{
    ResiduePlaceholder start, end;
    int disappear_iter = -1;
    int reappear_iter = -1;
    float displacement = 0;
    Point dispcen;

    void update(Protein* prot, int effective_iter);
};

#endif
