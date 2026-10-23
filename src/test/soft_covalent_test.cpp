#include <cassert>
#include <cstdio>
#include "../classes/soft.h"

int main()
{
    Protein protein;
    protein.add_sequence("AKA");
    AminoAcid* lysine = protein.get_residue(2);
    Molecule ligand("acetaldehyde", "CC=O");
    Molecule* water = ligand.create_Schiff_base(lysine);
    assert(water);
    delete water;
    ligand.check_glued_bond();
    Pose ligand_was(&ligand);
    Point original = ligand.get_atom(0)->loc;
    Vector step = Point(4,0,0);
    protein.move_piece(1,3,step);
    // A nested protein operation can replace the undo snapshot after this move.
    protein.move_piece(1,3,step);
    undo_soft_translation(&protein,&ligand,ligand_was);
    ligand.check_glued_bond();
    assert(ligand.get_atom(0)->loc.get_3d_distance(original.add(Point(4,0,0))) < 0.00001);

    Molecule free_ligand("water", "O");
    Pose free_was(&free_ligand);
    Point free_original = free_ligand.get_atom(0)->loc;
    Point protein_original = lysine->get_CA_location();
    protein.move_piece(1,3,step);
    free_ligand.move(step);
    free_ligand.movability = MOV_PINNED;
    undo_soft_translation(&protein,&free_ligand,free_was);
    assert(free_ligand.get_atom(0)->loc.get_3d_distance(free_original) < 0.00001);
    assert(lysine->get_CA_location().get_3d_distance(protein_original) < 0.00001);
    ligand.check_glued_bond();
    puts("Soft rollback preserves covalent bonds and restores free ligands.");
}
