#include <cassert>
#include <cstdio>
#include "../classes/molecule.h"
#include "../classes/cavity.h"

int main()
{
    Atom named("C");
    named.name = const_cast<char*>("C1");
    named.region = const_cast<char*>("TM6");
    Atom copied(named);
    assert(!strcmp(copied.name, "C1"));
    assert(!strcmp(copied.region, "TM6"));
    assert(copied.name != named.name && copied.region != named.region);
    named.name = named.region = nullptr;
    Atom unnamed(named);
    assert(unnamed.name[0] == 0 && unnamed.region[0] == 0);
    Molecule empty;
    Molecule empty_clone(empty);
    assert(empty_clone.get_atom_count() == 0);
    assert(empty_clone.get_name() == nullptr);
    Molecule m;
    m.from_smiles("CC(=O)O", false);
    m.set_name("copy-test");
    Molecule clone(m);
    assert(!strcmp(clone.get_name(), m.get_name()));
    assert(clone.get_atom_count() == m.get_atom_count());
    assert(clone.get_bond_count(true) == m.get_bond_count(true));
    for (int i=0; i<m.get_atom_count(); i++)
    {
        assert(clone.get_atom(i) != m.get_atom(i));
        assert(clone.get_atom(i)->mol == &clone);
        assert(clone.get_atom(i)->loc.get_3d_distance(m.get_atom(i)->loc) < 0.001);
    }
    Point original = m.get_atom(0)->loc;
    clone.get_atom(0)->move(original.add(Point(1,0,0)));
    assert(m.get_atom(0)->loc.get_3d_distance(original) < 0.001);
    for (int i=0; i<100; i++)
    {
        m.allocate_mandatory_connections(2);
        m.delete_mandatory_connections();
    }
    Cavity a, b;
    SPartial part;
    part.s.center = Point(1,2,3);
    part.s.radius = 3;
    a.add_partial(part);
    b.add_partial(part);
    a.unify(&b);
    assert(a.count_partials() == 2);
    assert(b.count_partials() == 0);
    a.unify(&a);
    assert(a.count_partials() == 2);
    puts("Allocation checks passed.");
}
