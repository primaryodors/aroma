#include <cassert>
#include <cstdio>
#include "../classes/protein.h"

int main()
{
    Protein protein;
    protein.add_sequence("CAC");
    AminoAcid* first = protein.get_residue(1);
    AminoAcid* second = protein.get_residue(3);
    Atom* sg1 = first->get_atom("SG");
    Atom* sg2 = second->get_atom("SG");
    assert(sg1 && sg2);
    second->movability = MOV_ALL;
    second->aamove(sg1->loc.add(Point(2.07, 0, 0)).subtract(sg2->loc));
    assert(fabs(sg1->distance_to(sg2) - 2.07) < 0.001);
    assert(sg1->is_bonded_to("H") && sg2->is_bonded_to("H"));
    assert(protein.disulfide_bond(1, 3));
    assert(sg1->get_bond_between(sg2));
    assert(!sg1->is_bonded_to("H") && !sg2->is_bonded_to("H"));
    auto links = protein.get_disulfide_residues();
    assert(links.size() == 1);
    assert(links[0].first == 1 && links[0].second == 3);
    puts("Cysteine sulfur atoms form a declared disulfide.");
}
