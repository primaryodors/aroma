#include <cassert>
#include <cstdio>
#include "../classes/aminoacid.h"

int main()
{
    AminoAcid residue('K');
    residue.movability = MOV_FORCEFLEX;
    Atom* atom = residue.get_atom("NZ");
    assert(atom);
    const Point target = atom->loc;
    residue.Molecule::conform_atom_to_location("NZ", target, 1);
    assert(atom->loc.get_3d_distance(target) == 0);

    AminoAcid partner('D');
    Atom* other = partner.get_atom("OD1");
    assert(other);
    const float distance = atom->distance_to(other);
    residue.Molecule::conform_atom_to_location(atom, other, 1, distance);
    assert(atom->loc.get_3d_distance(target) == 0);
    residue.conform_atom_to_location("NZ", target, 1, 0, true);
    assert(atom->loc.get_3d_distance(target) == 0);
    residue.conform_atom_to_location(atom, other, 1, distance, true);
    assert(atom->loc.get_3d_distance(target) == 0);
    puts("Atom targeting retains an optimal starting pose.");
}
