#include <cassert>
#include <cmath>
#include <cstdio>
#include "../classes/protein.h"

int main()
{
    Protein protein;
    protein.add_sequence("ASCA");
    protein.set_clashables();
    const MovabilityType flags[] = {MOV_PINNED,MOV_FLEXONLY,MOV_FLXDESEL,MOV_BKGRND};
    Point positions[4];
    for (int i=0; i<4; ++i)
    {
        AminoAcid* residue = protein.get_residue(i+1);
        residue->movability = flags[i];
        positions[i] = residue->get_CA_location();
    }
    float clashes = protein.get_internal_clashes(1,4,true,0);
    assert(std::isfinite(clashes));
    for (int i=0; i<4; ++i)
    {
        AminoAcid* residue = protein.get_residue(i+1);
        if (residue->movability != flags[i])
        {
            fprintf(stderr,"Repacking changed residue %d motion flags from %x to %x.\n",
                    i+1,flags[i],residue->movability);
            return 1;
        }
        assert(residue->get_CA_location().get_3d_distance(positions[i]) < 0.00001);
    }
    puts("Repacking preserves residue motion constraints.");
}
