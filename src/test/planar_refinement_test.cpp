#include <cassert>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <vector>
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

struct LocalDistance
{
    Atom* a;
    Atom* b;
    float distance;
};

static std::vector<LocalDistance> local_geometry(Molecule& ligand)
{
    std::vector<LocalDistance> result;
    for (int i=0; i<ligand.get_atom_count(); ++i)
    {
        Atom* center = ligand.get_atom(i);
        Bond* bonds[32] = {};
        center->fetch_bonds(bonds);
        for (int j=0; bonds[j]; ++j)
        {
            Atom* a = bonds[j]->atom2;
            result.push_back({center, a, center->distance_to(a)});
            for (int k=0; k<j; ++k)
            {
                Atom* b = bonds[k]->atom2;
                result.push_back({a, b, a->distance_to(b)});
            }
        }
    }
    return result;
}

static void check_geometry(const std::vector<LocalDistance>& original)
{
    for (const auto& pair : original)
    {
        float after = pair.a->distance_to(pair.b);
        if (fabs(after-pair.distance) >= 0.001)
            fprintf(stderr, "%s-%s changed from %.6f to %.6f Angstrom.\n",
                    pair.a->name, pair.b->name, pair.distance, after);
        assert(fabs(after-pair.distance) < 0.001);
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

    Molecule prepared;
    prepared.from_sdf(buffer.str().c_str());
    auto geometry = local_geometry(prepared);
    prepared.optimize();
    check_geometry(geometry);
    prepared.minimize_internal_clashes();
    check_geometry(geometry);
    check_ring(prepared, original);
    Bond** bonds = prepared.get_rotatable_bonds();
    for (int step=0; step<20; ++step)
        for (int i=0; bonds[i]; ++i)
        {
            bonds[i]->rotate(0.4);
            check_geometry(geometry);
            bonds[i]->get_reversed()->rotate(-0.3);
            check_geometry(geometry);
        }
    puts("Ring, bond lengths, and bond angles survive docking preparation and torsions.");
}
