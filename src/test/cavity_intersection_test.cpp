#include <cassert>
#include <cmath>
#include <cstdio>
#include "../classes/cavity.h"

int main()
{
    Protein protein;
    for (int i=1; i<=300; ++i) assert(protein.add_residue(i,'A'));
    SPartial sphere;
    sphere.s.center = Point(1,2,3);
    sphere.s.radius = 10000;
    Cavity first, second;
    first.add_partial(sphere);
    second.add_partial(sphere);
    assert(first.get_partial_by_idx(0));
    assert(second.get_partial_by_idx(0));
    assert(std::fabs(first.residue_intersection(&second,&protein)-1) < 0.00001);
    Cavity distant;
    sphere.s.center = Point(1000000,1000000,1000000);
    sphere.s.radius = 4;
    distant.add_partial(sphere);
    assert(distant.get_partial_by_idx(0));
    assert(distant.residue_intersection(&first,&protein) == 0);
    assert(first.residue_intersection(nullptr,&protein) == 0);
    assert(first.residue_intersection(&second,nullptr) == 0);
    puts("Cavity residue intersections support large and empty residue sets.");
}
