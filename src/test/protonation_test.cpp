#include <cassert>
#include <cstdio>
#include "../classes/aminoacid.h"

static void check_state(bool present)
{
    AminoAcid histidine('H'), acid('D');
    Atom *bare = nullptr, *donor = nullptr, *hydrogen = nullptr;
    for (int i=0; i<histidine.get_atom_count(); ++i)
    {
        Atom* atom = histidine.get_atom(i);
        if (atom->is_backbone || atom->get_family() != PNICTOGEN) continue;
        Atom* attached = atom->is_bonded_to("H");
        if (attached) { donor = atom; hydrogen = attached; }
        else bare = atom;
    }
    assert(bare && donor && hydrogen);
    if (present)
    {
        histidine.add_atom("H", "extra", bare, 1);
        bare->increment_charge(1);
    }
    Point direction = hydrogen->loc.subtract(donor->loc);
    direction.scale(1.8);
    Point target = hydrogen->loc.add(direction);
    acid.aamove(target.subtract(acid.get_atom("OD1")->loc));
    Molecule* nearby[] = {present ? &acid : nullptr, nullptr};
    const int count = histidine.get_atom_count();
    for (int i=0; i<3; ++i)
    {
        histidine.set_conditional_basicity(nearby);
        assert(histidine.get_atom_count() == count);
        assert(bare->num_bonded_to("H") == (present ? 1 : 0));
        assert(bare->get_charge() == (present ? 1 : 0));
    }
}

int main()
{
    check_state(false);
    check_state(true);
    puts("Conditional protonation is idempotent.");
}
