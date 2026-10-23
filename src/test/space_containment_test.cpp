#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include "../classes/space.h"

void check_falloff(bool paired)
{
    Space cavity;
    SPartial partial;
    partial.s.center = Point(5,0,0);
    partial.s.radius = 4;
    cavity.add_partial(partial);
    if (paired)
    {
        partial.s.center = Point(7,0,0);
        cavity.add_partial(partial);
    }
    Sphere atom;
    atom.radius = 1;
    float inner = partial.s.radius-atom.radius+global_clash_allowance;
    if (paired) inner *= 1.414f;
    float previous = 1;
    for (float offset : {0.f,0.1f,0.5f,0.9f,1.f,1.1f})
    {
        atom.center = paired ? Point(6,inner+offset,0) : Point(5+inner+offset,0,0);
        float containment = cavity.sphere_inside_pocket(atom);
        assert(std::isfinite(containment));
        assert(containment >= 0 && containment <= 1);
        assert(containment <= previous+0.00001f);
        assert(std::fabs(containment-std::max(0.f,1-offset)) < 0.00001f);
        previous = containment;
    }
}

int main()
{
    check_falloff(false);
    check_falloff(true);
    puts("Sphere containment decreases continuously across cavity boundaries.");
}
