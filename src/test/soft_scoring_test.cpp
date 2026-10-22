#include <cassert>
#include <cstdio>
#include <vector>
#include "../classes/soft.h"

struct ContactFixture : SoftRegion
{
    void set_reference_energy(float value) { contacts[0].energy = value; }
};

int main()
{
    Protein protein;
    FILE* input = fopen("pdbs/OR1/OR1D2.active.pdb", "rb");
    assert(input);
    protein.load_pdb(input);
    fclose(input);
    AminoAcid* first = protein.get_residue(180);
    AminoAcid* second = protein.get_residue(259);
    assert(first && second);
    first->movability = second->movability = MOV_ALL;
    ContactFixture region;
    region.add_contact(180, 259, &protein);
    // An improved contact must remain unchanged while its penalty is measured.
    region.set_reference_energy(100);
    std::vector<Atom*> atoms;
    std::vector<Point> coordinates;
    for (AminoAcid* residue : {first, second})
        for (int i=0; i<residue->get_atom_count(); ++i)
        {
            atoms.push_back(residue->get_atom(i));
            coordinates.push_back(atoms.back()->loc);
        }
    for (int i=0; i<3; ++i)
    {
        assert(region.contact_anomaly(&protein) < 0.001f);
        for (size_t j=0; j<atoms.size(); ++j)
            assert(atoms[j]->loc.get_3d_distance(coordinates[j]) < 0.000001f);
    }
    puts("Soft scoring preserves coordinates.");
}
