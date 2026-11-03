#include <cassert>
#include <cstdio>
#include "../classes/pose_contact.h"

int main()
{
    Protein original, current, missing;
    original.add_sequence("AKA");
    current.add_sequence("AKA");
    missing.add_sequence("AAA");
    Molecule old_ligand("ethanol", "CCO"), new_ligand("ethanol", "CCO");
    old_ligand.stay_close_mine = old_ligand.get_atom("O3");
    old_ligand.stay_close_other = original.get_residue(2)->get_atom("NZ");
    old_ligand.stay_close2_mine = old_ligand.get_atom("C1");
    old_ligand.stay_close2_other = original.get_residue(2)->get_atom("CE");
    assert(old_ligand.stay_close_mine && old_ligand.stay_close_other);
    old_ligand.stay_close_tolerance = 0.4;
    old_ligand.stay_close_optimal = 2.7;
    old_ligand.stay_close2_optimal = 3.2;
    PoseContactSnapshot saved;
    saved.capture(&old_ligand);
    saved.restore(&new_ligand, &current);
    assert(new_ligand.stay_close_mine == new_ligand.get_atom("O3"));
    assert(new_ligand.stay_close_other == current.get_residue(2)->get_atom("NZ"));
    assert(new_ligand.stay_close_mine != old_ligand.stay_close_mine);
    assert(new_ligand.stay_close_other != old_ligand.stay_close_other);
    assert(new_ligand.stay_close_mol == current.get_residue(2));
    assert(new_ligand.stay_close2_mine == new_ligand.get_atom("C1"));
    assert(new_ligand.stay_close2_other == current.get_residue(2)->get_atom("CE"));
    assert(new_ligand.stay_close2_mol == current.get_residue(2));
    assert(new_ligand.stay_close_tolerance == old_ligand.stay_close_tolerance);
    assert(new_ligand.stay_close_optimal == old_ligand.stay_close_optimal);
    assert(new_ligand.stay_close2_optimal == old_ligand.stay_close2_optimal);
    saved.restore(&new_ligand, &missing);
    assert(!new_ligand.stay_close_mine && !new_ligand.stay_close_other && !new_ligand.stay_close_mol);
    assert(!new_ligand.stay_close2_mine && !new_ligand.stay_close2_other && !new_ligand.stay_close2_mol);
    puts("Cached contacts refer to the current ligand and receptor atoms.");
}
