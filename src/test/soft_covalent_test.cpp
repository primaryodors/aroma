#include <cassert>
#include <cstdio>
#include "../classes/soft.h"
#include "../classes/scoring.h"

struct CovalentLigand : Molecule
{
    using Molecule::Molecule;
    using Molecule::glued_to;
    using Molecule::glued_energy;
};

static void verify_adduct_scoring(CovalentLigand& ligand, AminoAcid* lysine)
{
    assert(ligand.glued_to == lysine);
    assert(fabs(ligand.glued_energy.summed() - intermol_covalent_enthalpy) < 0.001);
    for (bool reverse : {false, true})
    {
        Molecule* first = reverse ? static_cast<Molecule*>(lysine) : &ligand;
        Molecule* second = reverse ? &ligand : static_cast<Molecule*>(lysine);
        ninterall = 0;
        compute_interall = true;
        Interaction score = first->get_intermol_binding(second, false);
        compute_interall = false;
        int covalent_count = 0;
        for (int i=0; i<ninterall; i++)
        {
            assert(interall_a1[i]->mol != interall_a2[i]->mol);
            if (interall_t[i] == covalent)
            {
                covalent_count++;
                assert(fabs(interall[i] - ligand.glued_energy.summed()) < 0.001);
            }
        }
        assert(covalent_count == 1);
        Interaction formation = ligand.glued_energy;
        ligand.glued_energy = 0;
        Interaction without = first->get_intermol_binding(second, false);
        ligand.glued_energy = formation;
        assert(fabs(score.summed() - without.summed() - formation.summed()) < 0.001);
    }
}

int main()
{
    Protein protein;
    protein.add_sequence("AKA");
    AminoAcid* lysine = protein.get_residue(2);
    CovalentLigand ligand("acetaldehyde", "CC=O");
    lysine->been_flexed = ligand.been_flexed = false;
    Molecule* water = ligand.create_Schiff_base(lysine);
    assert(water);
    assert(lysine->been_flexed && ligand.been_flexed);
    delete water;
    ligand.check_glued_bond();
    verify_adduct_scoring(ligand, lysine);
    Pose anchor_placement(lysine);
    Point ligand_position = ligand.get_atom(0)->loc;
    Bond** rotations = lysine->get_rotatable_bonds();
    assert(rotations && rotations[0]);
    rotations[0]->rotate(0.5);
    anchor_placement.restore_state(lysine);
    ligand.check_glued_bond();
    assert(ligand.get_atom(0)->loc.get_3d_distance(ligand_position) < 0.00001);
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
    lysine->been_flexed = false;
    lysine->is_ic_res = false;
    AminoAcid* no_reaches[] = {nullptr};
    DockResult exported(&protein, &ligand, 0, no_reaches, nullptr, 1, nullptr, true);
    assert(lysine->been_flexed);
    assert(exported.pdbdat.find(" NZ ") != std::string::npos);

    Protein reverse_protein;
    reverse_protein.add_sequence("AKA");
    AminoAcid* reverse_lysine = reverse_protein.get_residue(2);
    CovalentLigand reverse_ligand("phenylacetaldehyde", "c1ccccc1CC=O");
    reverse_lysine->been_flexed = reverse_ligand.been_flexed = false;
    water = reverse_lysine->create_Schiff_base(&reverse_ligand);
    assert(water);
    delete water;
    assert(reverse_lysine->been_flexed && reverse_ligand.been_flexed);
    reverse_ligand.check_glued_bond();
    verify_adduct_scoring(reverse_ligand, reverse_lysine);
    puts("Soft rollback preserves covalent bonds and restores free ligands.");
}
