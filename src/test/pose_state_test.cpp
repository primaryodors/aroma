#include <cassert>
#include <cstdio>
#include <utility>
#include "../classes/molecule.h"

struct InspectPose : Pose
{
    const Point* storage() const { return saved_atom_locs; }
};

int main()
{
    Molecule m;
    m.from_smiles("CCO", false);
    InspectPose p;
    p.copy_state(&m);
    assert(p.has_data());
    const Point* storage = p.storage();
    for (int i=0; i<100; i++)
    {
        p.copy_state(&m);
        assert(p.storage() == storage);
    }
    Point original = m.get_atom(0)->loc;
    Pose copy(p);
    Pose assigned;
    for (int i=0; i<100; i++) assigned = copy;
    Pose moved(std::move(assigned));
    assert(!assigned.has_data());
    assert(moved.has_data());
    m.get_atom(0)->move(original.add(Point(3,0,0)));
    moved.restore_state(&m);
    assert(m.get_atom(0)->loc.get_3d_distance(original) < 0.001);
    p.copy_state(nullptr);
    assert(!p.has_data());
    p.copy_state(&m);
    p.reset();
    assert(!p.has_data());
    puts("Pose storage passed.");
}
