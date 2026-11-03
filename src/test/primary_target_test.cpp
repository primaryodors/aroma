#include <cassert>
#include <cstdio>
#include "../classes/search.h"

int main()
{
    Protein protein;
    protein.add_sequence("SA");
    Molecule ligand;
    ligand.from_smiles("O");
    LigandTarget primary;
    primary.single_atom = ligand.get_atom(0);
    BestBindingResult targets;
    targets.protein = &protein;
    targets.ligand = &ligand;
    targets.pri_res = protein.get_residue(1);
    targets.pri_tgt = &primary;
    Point center = targets.pri_res->get_CA_location().add(Point(3,2,1));
    Search::align_targets(&ligand,center,&targets);
    for (int i=0; i<ligand.get_atom_count(); ++i)
    {
        Point location = ligand.get_atom(i)->loc;
        assert(std::isfinite(location.x));
        assert(std::isfinite(location.y));
        assert(std::isfinite(location.z));
    }
    Point position = ligand.get_atom(0)->loc;
    targets.pri_tgt = nullptr;
    Search::align_targets(&ligand,center,&targets);
    assert(ligand.get_atom(0)->loc.get_3d_distance(position) == 0);
    puts("A single ligand target aligns without a secondary pair.");
}
