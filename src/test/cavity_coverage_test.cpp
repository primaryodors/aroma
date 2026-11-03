#include <cassert>
#include <cmath>
#include <cstdio>
#include "../classes/cavity.h"

int main()
{
    Molecule ligand;
    for (float x : {0.f, 1.f, 20.f, 21.f})
    {
        Point location(x,0,0);
        ligand.add_atom("C", "C", &location, nullptr, 0);
    }
    Point hydrogen_location(100,100,100);
    ligand.add_atom("H", "H", &hydrogen_location, nullptr, 0);
    Cavity cavities[3];
    CPartial sphere;
    sphere.s.center = Point(0.5,0,0);
    sphere.s.radius = 6;
    cavities[0].add_partial(sphere);
    cavities[1].add_partial(sphere);
    float half = Cavity::heavy_atom_containment(&ligand, cavities, 1);
    assert(std::fabs(half-0.5f) < 0.00001);
    assert(half < min_cvty_ctnmt);
    assert(Cavity::heavy_atom_containment(&ligand, cavities, 2) == half);
    sphere.s.center = Point(20.5,0,0);
    cavities[2].add_partial(sphere);
    assert(Cavity::heavy_atom_containment(&ligand, cavities, 3) == 1);
    assert(Cavity::heavy_atom_containment(nullptr, cavities, 3) == 0);
    assert(Cavity::heavy_atom_containment(&ligand, nullptr, 3) == 0);
    assert(Cavity::heavy_atom_containment(&ligand, cavities, 0) == 0);
    Molecule empty;
    assert(Cavity::heavy_atom_containment(&empty, cavities, 3) == 0);
    puts("Cavity coverage rejects exposed heavy atoms without counting overlaps twice.");
}
