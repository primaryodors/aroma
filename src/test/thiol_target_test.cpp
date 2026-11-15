#include "../classes/search.h"
#include <cassert>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <vector>

static void same_bonds(Molecule& ligand,const std::vector<float>& lengths)
{
    int k=0;
    for (int i=0; i<ligand.get_atom_count(); i++) for (int j=i+1; j<ligand.get_atom_count(); j++)
        if (ligand.get_atom(i)->is_bonded_to(ligand.get_atom(j)))
            assert(fabs(ligand.get_atom(i)->distance_to(ligand.get_atom(j))-lengths[k++])<.002);
    assert(k==int(lengths.size()));
}

int main()
{
    Protein protein;
    protein.add_sequence("YQ");
    AminoAcid* tyrosine=protein.get_residue(1);
    AminoAcid* glutamine=protein.get_residue(2);
    Molecule thiol,hydrocarbon,alcohol;
    thiol.from_smiles("CCCCCCCCS",false);
    hydrocarbon.from_smiles("CCCCCCCCC",false);
    alcohol.from_smiles("CCCCCCCCO",false);
    Atom* sulfur=thiol.get_atom("S9");
    Atom* hydrogen=sulfur->is_bonded_to("H");
    assert(sulfur && hydrogen);
    LigandTarget target,hydrogen_target,control;
    target.single_atom=sulfur; hydrogen_target.single_atom=hydrogen;
    control.single_atom=hydrocarbon.get_atom("C9");
    float radius,strength;
    assert(InteratomicForce::thiol_pi_parameters(sulfur,tyrosine->get_atom("CE1"),radius,strength));
    assert(fabs(radius-2.71)<.0001 && fabs(strength-8.3)<.0001);
    assert(fabs(Search::thiol_pi_strength(tyrosine,&target)-strength)<.0001);
    assert(fabs(Search::thiol_pi_strength(tyrosine,&hydrogen_target)-strength)<.0001);
    assert(Search::thiol_pi_strength(tyrosine,&control)==0);
    assert(Search::thiol_pi_strength(glutamine,&target)==0);

    Atom *ligand_atom,*protein_atom;
    float contact_radius;
    assert(!Search::contact_atoms(glutamine,&target,ligand_atom,protein_atom,contact_radius));
    LigandTarget hydroxyl;
    hydroxyl.single_atom=alcohol.get_atom("O9");
    Atom* oh=hydroxyl.single_atom->is_bonded_to("H");
    assert(InteratomicForce::hydrogen_bond_distance(oh,glutamine->get_atom("NE2"))==0);
    assert(Search::contact_atoms(glutamine,&hydroxyl,ligand_atom,protein_atom,contact_radius));
    if (ligand_atom->Z==1) assert(protein_atom==glutamine->get_atom("OE1"));
    else assert(protein_atom->Z==1 && ligand_atom==hydroxyl.single_atom);

    std::vector<float> lengths;
    for (int i=0; i<thiol.get_atom_count(); i++) for (int j=i+1; j<thiol.get_atom_count(); j++)
        if (thiol.get_atom(i)->is_bonded_to(thiol.get_atom(j)))
            lengths.push_back(thiol.get_atom(i)->distance_to(thiol.get_atom(j)));
    BestBindingResult result;
    result.protein=&protein; result.ligand=&thiol; result.pri_res=tyrosine; result.pri_tgt=&target;
    loneliest=tyrosine->get_barycenter();
    assert(result.score(thiol.get_barycenter())>0);
    Search::align_targets(&thiol,loneliest,&result);
    assert(Search::contact_atoms(tyrosine,&target,ligand_atom,protein_atom,contact_radius));
    assert(ligand_atom==hydrogen && protein_atom->Z==6 && protein_atom->is_pi());
    float distance=hydrogen->distance_to(protein_atom);
    float angle=find_3d_angle(sulfur->loc,protein_atom->loc,hydrogen->loc)*fiftyseven;
    printf("Thiol aromatic anchor H-C %g A; S-H-C %g degrees\n",distance,angle);
    assert(fabs(distance-radius)<.02);
    assert(angle>170);
    same_bonds(thiol,lengths);
    Search::align_targets(&thiol,loneliest,&result,.2);
    same_bonds(thiol,lengths);
    // A sulfur anchor in a secondary slot gets the same typed alignment.
    LigandTarget tail; tail.single_atom=thiol.get_atom("C1");
    result.pri_res=glutamine; result.pri_tgt=&tail;
    result.sec_res=tyrosine; result.sec_tgt=&target;
    Search::align_targets(&thiol,loneliest,&result);
    assert(Search::contact_atoms(tyrosine,&target,ligand_atom,protein_atom,contact_radius));
    assert(fabs(ligand_atom->distance_to(protein_atom)-radius)<.02);
    same_bonds(thiol,lengths);
    // An explicit copper site retains precedence over the aromatic alternative.
    Point ion_position=tyrosine->get_barycenter().add(Point(20,20,20));
    Atom copper("Cu",&ion_position,1);
    tyrosine->coordmtl=&copper;
    assert(Search::contact_atoms(tyrosine,&target,ligand_atom,protein_atom,contact_radius));
    assert(ligand_atom==sulfur && protein_atom==&copper && fabs(contact_radius-2.29)<.001);
    result.pri_res=tyrosine; result.pri_tgt=&target;
    result.sec_res=nullptr; result.sec_tgt=nullptr;
    Search::align_targets(&thiol,ion_position,&result);
    assert(fabs(sulfur->distance_to(&copper)-2.29)<.01);
    same_bonds(thiol,lengths);
    tyrosine->coordmtl=nullptr;
    puts("Reference sulfur parameters, aromatic alignment, and amide donor roles pass.");
}
