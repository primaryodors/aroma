#include <cassert>
#include <cmath>
#include <cstdio>
#include "../classes/cavity.h"

int main()
{
    SPartial sphere;
    sphere.s.center = Point(1,2,3);
    sphere.s.radius = 4;
    Cavity growing;
    growing.add_partial(sphere);
    assert(growing.get_partial_by_idx(0));
    for (int i=0; i<512; ++i)
    {
        sphere.s.center = Point(i+10,2,3);
        growing.add_partial(sphere);
    }
    Space* base = &growing;
    assert(static_cast<void*>(growing.get_partial_by_idx(0)) ==
           static_cast<void*>(base->get_partial_by_idx(0)));
    assert(growing.get_partial_by_idx(0)->s.center.get_3d_distance(Point(1,2,3)) < 0.00001);
    assert(growing.get_partial_by_idx(-1) == nullptr);
    assert(growing.get_partial_by_idx(growing.count_partials()) == nullptr);

    Cavity first, second, priority, empty;
    sphere.s.center = Point(1,2,3);
    first.add_partial(sphere);
    sphere.s.center = Point(3,2,3);
    second.add_partial(sphere);
    float overlap = first.cavity_intersection(&second);
    assert(overlap > 0);
    assert(std::fabs(overlap-second.cavity_intersection(&first)) < 0.00001);
    sphere.priority = true;
    priority.add_partial(sphere);
    assert(priority.has_priority());
    assert(!empty.has_priority());
    assert(empty.get_partial_by_idx(0) == nullptr);
    puts("Cavity queries use current storage after allocation growth.");
}
