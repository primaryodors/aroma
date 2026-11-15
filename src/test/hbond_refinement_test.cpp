#include <cassert>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <vector>
#include "../classes/protein.h"

static void load(Protein& protein, Molecule& ligand)
{
    FILE* file = fopen("testdata/serine_alcohol_contact.pdb","r");
    assert(file);
    protein.load_pdb(file);
    rewind(file);
    Molecule positions;
    positions.from_pdb(file,true);
    fclose(file);
    std::ifstream sdf("testdata/hexenol_contact.sdf");
    std::stringstream text; text << sdf.rdbuf();
    ligand.from_sdf(text.str().c_str());
    for (int i=0; i<ligand.get_atom_count(); i++)
    {
        Atom* atom = ligand.get_atom(i);
        Atom* source = positions.get_atom(atom->name);
        assert(source);
        atom->move(source->loc);
    }
}

int main()
{
    Protein protein;
    Molecule ligand;
    load(protein,ligand);
    AminoAcid* ser = protein.get_residue(206);
    Atom* oxygen = ligand.get_atom("O7");
    Atom* hydrogen = ligand.get_atom("H19");
    Atom* og = ser->get_atom("OG");
    float radius = InteratomicForce::hydrogen_bond_distance(hydrogen,og);
    assert(fabs(radius-1.97)<.0001);
    float heavy_radius = radius + InteratomicForce::covalent_bond_radius(oxygen,hydrogen,1);
    assert(fabs(heavy_radius-2.94)<.0001);
    assert(InteratomicForce::hydrogen_bond_distance(oxygen,og)==0);
    assert(InteratomicForce::hydrogen_bond_distance(hydrogen,ligand.get_atom("C1"))==0);
    Molecule ammonia;
    ammonia.from_smiles("N",false);
    Atom* nh = ammonia.get_atom("N1")->is_bonded_to("H");
    assert(nh);
    assert(fabs(InteratomicForce::hydrogen_bond_distance(nh,og)-1.66)<.0001);
    assert(oxygen->distance_to(og)>4);
    Point ca = ser->get_atom("CA")->loc;
    Point cb = ser->get_atom("CB")->loc;
    Point ligand_heavy[7];
    for (int i=0; i<7; i++) ligand_heavy[i]=ligand.get_atom(i)->loc;
    Pose pinned(ser);
    ser->movability=MOV_PINNED;
    assert(!protein.refine_hbond_contact(&ligand,ser));
    assert(pinned.total_atom_motions()<.0001);
    ser->movability=MOV_FLEXONLY;
    assert(protein.refine_hbond_contact(&ligand,ser));
    printf("Refined O-O %g H-O %g donor angle %g\n",oxygen->distance_to(og),hydrogen->distance_to(og),
        find_3d_angle(oxygen->loc,og->loc,hydrogen->loc)*fiftyseven);
    assert(oxygen->distance_to(og)<3.3);
    assert(hydrogen->distance_to(og)<2.4);
    assert(find_3d_angle(oxygen->loc,og->loc,hydrogen->loc)*fiftyseven>135);
    assert(ser->get_atom("CA")->loc.get_3d_distance(ca)<.0001);
    assert(ser->get_atom("CB")->loc.get_3d_distance(cb)<.0001);
    for (int i=0; i<7; i++) assert(ligand.get_atom(i)->loc.get_3d_distance(ligand_heavy[i])<.0001);

    // The ligand participates in cached and uncached final hydrogen cleanup.
    Protein reference;
    Molecule reference_ligand;
    load(reference,reference_ligand);
    assert(reference.refine_hbond_contact(&reference_ligand,reference.get_residue(206)));
    protein.optimize_hydrogens(206,206,nullptr,true,&ligand);
    reference.optimize_hydrogens(206,206,nullptr,false,&reference_ligand);
    for (int i=0; i<ser->get_atom_count(); i++)
    {
        Atom* atom=ser->get_atom(i);
        assert(atom->distance_to(reference.get_atom(206,atom->name))<.0001);
    }
    assert(hydrogen->distance_to(og)<2.4);

    // A polar side chain can search a distant rotamer during ordinary docking.
    Protein searching;
    Molecule fixed_ligand;
    load(searching,fixed_ligand);
    AminoAcid* moving = searching.get_residue(206);
    Bond* hydroxyl = fixed_ligand.get_atom("C6")->get_bond_between(fixed_ligand.get_atom("O7"));
    hydroxyl->rotate(-M_PI/4,false,true);
    fixed_ligand.movability = MOV_PINNED;
    Pose fixed_pose(&fixed_ligand);
    // An explicit receptor list must not hide the ligand from the trial score.
    std::vector<Molecule*> background{nullptr};
    moving->mclashables = background.data();
    Bond* chi = moving->get_atom("CA")->get_bond_between(moving->get_atom("CB"));
    float initial_rotations = chi->total_rotations;
    Molecule* molecules[] = {&fixed_ligand,moving,nullptr};
    srand(11973);
    for (int trial=0; trial<10; trial++) Molecule::conform_molecules(molecules,1);
    // Broad exploration can still reject a contact under the docking clash weights.
    // Its angular coverage must exceed the ten bounded Monte Carlo proposals.
    printf("Side-chain angular coverage %g radians\n",chi->total_rotations-initial_rotations);
    assert(fabs(chi->total_rotations-initial_rotations)>2*M_PI);
    assert(fixed_pose.total_atom_motions()<.0001);
    moving->mclashables=nullptr;

    // An isolated Ser must turn its OH toward a ligand acceptor during cleanup.
    Protein isolated;
    assert(isolated.add_residue(1,'S'));
    AminoAcid* donor = isolated.get_residue(1);
    Atom* donor_o = donor->get_atom("OG");
    Atom* donor_h = donor->get_atom("HG");
    Point opposite = rotate3D(donor_h->loc,donor_o->loc,
                             donor_o->loc.subtract(donor->get_atom("CB")->loc),M_PI);
    Vector direction = opposite.subtract(donor_o->loc);
    direction.r = heavy_radius;
    Point destination = donor_o->loc.add(direction);
    Molecule acceptor;
    acceptor.from_smiles("C=O",false);
    Atom* acceptor_o = acceptor.get_atom("O2");
    Atom* acceptor_c = acceptor.get_atom("C1");
    assert(acceptor_o && acceptor_c);
    Rotation rotation = align_points_3d(acceptor_c->loc,acceptor_o->loc.add(direction),acceptor_o->loc);
    acceptor.rotate(rotation,acceptor_o->loc);
    acceptor.move(destination.subtract(acceptor_o->loc));
    float initial_distance = donor_h->distance_to(acceptor_o);
    isolated.optimize_hydrogens(1,1,nullptr,true,&acceptor);
    assert(donor_h->distance_to(acceptor_o)<2.2);
    assert(donor_h->distance_to(acceptor_o)<initial_distance-.5);
    puts("Reference distances, side-chain refinement, and ligand-aware hydrogen cleanup pass.");
}
