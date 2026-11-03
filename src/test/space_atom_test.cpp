#include <cassert>
#include <cmath>
#include <cstdio>
#include <limits>
#include "../classes/molecule.h"

int main()
{
    Space molecular;
    SPartial atom;
    atom.s.center = Point(1,2,3);
    atom.s.radius = 1.2;
    molecular.add_partial(atom);
    assert(molecular.count_partials() == 1);
    assert(molecular.get_partial_by_idx(0));
    assert(molecular.point_inside_pocket(atom.s.center));
    Sphere probe;
    probe.center = atom.s.center;
    probe.radius = 0.2;
    assert(molecular.sphere_inside_pocket(probe) == 1);

    Space origin;
    atom.s.center = Point(0,0,0);
    origin.add_partial(atom);
    assert(origin.count_partials() == 1);
    assert(origin.point_inside_pocket(Point(0,0,0)));
    atom.s.center = Point(1,2,3);
    atom.s.radius = 0;
    origin.add_partial(atom);
    atom.s.radius = std::numeric_limits<float>::quiet_NaN();
    origin.add_partial(atom);
    atom.s.radius = 1;
    atom.s.center.x = std::numeric_limits<double>::quiet_NaN();
    origin.add_partial(atom);
    assert(origin.count_partials() == 1);
    assert(origin.get_partial_by_idx(1) == nullptr);

    Molecule water("water", "O");
    Space* atoms = water.to_space();
    assert(atoms);
    assert(atoms->count_partials() == water.get_atom_count());
    delete atoms;
    puts("Molecular spaces retain finite atom-sized spheres at every position.");
}
