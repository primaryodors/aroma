#include <cassert>
#include <cstdio>
#include <vector>
#include "../classes/cavity.h"

int main()
{
    Protein protein;
    for (int i=1; i<=300; i++) assert(protein.add_residue(i, 'A'));
    Cavity empty;
    int absent[] = {-1};
    assert(empty.resnos_as_array(&protein, absent) == 0);
    assert(absent[0] == 0);

    SPartial sphere;
    sphere.s.center = Point(1,2,3);
    sphere.s.radius = 10000;
    Cavity cavity;
    cavity.add_partial(sphere);
    std::vector<int> residues(302, -1);
    const int count = cavity.resnos_as_array(&protein, residues.data());
    assert(count == 300);
    assert(residues[0] == 1);
    assert(residues[count-1] == protein.get_end_resno());
    assert(residues[count] == 0);
    for (int i=0; i<100; i++)
        assert(cavity.resnos_as_array(&protein, nullptr) == count);
    cavity.add_partial(sphere);
    assert(cavity.resnos_as_array(&protein, nullptr) == count);
    assert(cavity.resnos_as_array(&protein, residues.data()) == count);
    puts("Cavity residue counts passed.");
}
