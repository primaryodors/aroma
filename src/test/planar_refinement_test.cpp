#include <cassert>
#include <cstdio>
#include <fstream>
#include <sstream>
#include "../classes/molecule.h"

static void check_ring(Molecule& ligand, const Point* original)
{
    for (int i=0; i<6; ++i)
        for (int j=0; j<i; ++j)
        {
            float before = original[i].get_3d_distance(original[j]);
            float after = ligand.get_atom(i+2)->loc.get_3d_distance(ligand.get_atom(j+2)->loc);
            assert(fabs(after-before) < 0.001);
        }
}

int main()
{
    std::ifstream input("sdf/methyl_eugenol.sdf");
    assert(input);
    std::stringstream buffer;
    buffer << input.rdbuf();
    Molecule ligand;
    ligand.from_sdf(buffer.str().c_str());
    Point original[6];
    for (int i=0; i<6; ++i) original[i] = ligand.get_atom(i+2)->loc;
    srand(1);
    ligand.refine_structure();
    check_ring(ligand, original);
    ligand.optimize();
    ligand.minimize_internal_clashes();
    check_ring(ligand, original);
    puts("Planar ring geometry survives refinement and docking preparation.");
}
