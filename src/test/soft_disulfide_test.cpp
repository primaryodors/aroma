#include <cassert>
#include <cstdio>
#include "../classes/soft.h"

int main()
{
    Protein protein;
    protein.add_sequence("CAC");
    AminoAcid* first = protein.get_residue(1);
    AminoAcid* second = protein.get_residue(3);
    first->movability = second->movability = MOV_ALL;
    Atom* sg1 = first->get_atom("SG");
    Atom* sg2 = second->get_atom("SG");
    assert(sg1 && sg2);
    second->aamove(sg1->loc.add(Point(2.07,0,0)).subtract(sg2->loc));
    SoftRegion region;
    region.add_contact(1,3,&protein,false,true);
    assert(region.check_disulfide_constraints(&protein));
    first->aamove(Point(1,2,3));
    second->aamove(Point(1,2,3));
    assert(region.check_disulfide_constraints(&protein));
    for (float distance : {2.77f,3.3f,1.5f})
    {
        second->aamove(sg1->loc.add(Point(distance,0,0)).subtract(sg2->loc));
        assert(!region.check_disulfide_constraints(&protein));
    }
    puts("Soft motion preserves disulfides while allowing rigid movement.");
}
