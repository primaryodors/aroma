#include <cassert>
#include <cstdio>
#include "../classes/molecule.h"
#include "../classes/cavity.h"

int main()
{
    Molecule m;
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
