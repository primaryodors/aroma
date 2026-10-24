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
    SPartial near, far, middle;
    near.s.center = Point(1,0,0);
    far.s.center = Point(1,20,0);
    middle.s.center = Point(7,3,0);
    near.s.radius = far.s.radius = middle.s.radius = 3;
    Space first, second;
    first.add_partial(near);
    first.add_partial(far);
    first.add_partial(middle);
    second.add_partial(near);
    second.add_partial(middle);
    second.add_partial(far);
    Sphere atom;
    atom.center = Point(1,3,0);
    atom.radius = 1.7;
    float score = first.sphere_inside_pocket(atom);
    assert(score > 0 && score < 1);
    assert(std::fabs(score-second.sphere_inside_pocket(atom)) < 0.00001);
    puts("Sphere containment decreases continuously and ignores sphere insertion order.");
}
