#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>
#include "../classes/search.h"

int main()
{
    assert(fabs(sphere_intersection(1,2,0)-(4.0f/3.0f)*M_PI)<0.0001f);
    assert(fabs(sphere_intersection(1,1,0)-(4.0f/3.0f)*M_PI)<0.0001f);
    assert(sphere_intersection(1,1,2)==0);
    Protein protein;
    protein.add_sequence("K");
    AminoAcid* lysine = protein.get_residue(1);
    Molecule ligand("octanal","CCCCCCCC=O");
    assert(Search::sample_anchored_conformers(&protein,&ligand)==0);
    Atom *carbon=nullptr,*oxygen=nullptr,*nitrogen=nullptr;
    assert(ligand.identify_Schiff_carbonyl(&carbon,&oxygen));
    assert(lysine->identify_Schiff_amine(&nitrogen,nullptr,nullptr));
    Molecule* water = ligand.create_Schiff_base(lysine);
    assert(water && ligand.glued_to_mol()==lysine);

    std::vector<Atom*> fixed;
    std::vector<Point> positions;
    for (int i=0; i<lysine->get_atom_count(); ++i)
    {
        Atom* atom = lysine->get_atom(i);
        if (atom->mol==&ligand) continue;
        fixed.push_back(atom);
        positions.push_back(atom->loc);
    }
    Point anchor_position = carbon->loc;
    float imine_length = carbon->distance_to(nitrogen);
    struct Length { Atom* a; Atom* b; float distance; };
    std::vector<Length> lengths;
    for (int i=0; i<ligand.get_atom_count(); ++i)
    {
        Atom* atom = ligand.get_atom(i);
        for (int j=0; j<abs(atom->get_geometry()); ++j)
        {
            Bond* bond = atom->get_bond_by_idx(j);
            if (bond && bond->atom2)
                lengths.push_back({atom,bond->atom2,atom->distance_to(bond->atom2)});
        }
    }
    Space pocket;
    for (int i=0; i<ligand.get_atom_count(); ++i)
    {
        Atom* atom = ligand.get_atom(i);
        if (atom->Z<=1) continue;
        SPartial part;
        part.s.center = atom->loc;
        part.s.radius = 2.2;
        pocket.add_partial(part);
    }
    Molecule* nearby[] = {lysine,nullptr};
    // Fold a distal fragment out of the reference-shaped pocket. The sampler
    // must recover a better fit without moving the receptor or the imine.
    Pose reference(&ligand), folded(&ligand);
    float worst_fit = ligand.contained_by_space(&pocket);
    for (int i=0; i<ligand.get_atom_count(); ++i)
    {
        Atom* atom = ligand.get_atom(i);
        if (atom->mol!=&ligand || atom->Z<=1) continue;
        for (int j=0; j<abs(atom->get_geometry()); ++j)
        {
            Bond* bond = atom->get_bond_by_idx(j);
            if (!bond || !bond->atom2 || bond->atom2->mol!=&ligand
                || bond->atom2->Z<=1 || !bond->can_rotate || bond->cardinality!=1) continue;
            std::vector<Atom*> moving(bond->count_moves_with_atom2()+2,nullptr);
            bond->fetch_moves_with_atom2(moving.data());
            bool safe = moving[0]!=nullptr && bond->atom2!=carbon;
            for (Atom* part : moving)
                if (part && (part->mol!=&ligand || part==carbon)) safe = false;
            if (!safe) continue;
            reference.restore_state(&ligand,true);
            if (!bond->rotate(M_PI/3,false,true)) continue;
            float fit = ligand.contained_by_space(&pocket);
            if (fit<worst_fit)
            {
                worst_fit = fit;
                folded.copy_state(&ligand);
            }
        }
    }
    folded.restore_state(&ligand,true);
    float initial = ligand.conformer_binding(nearby,&pocket).summed();
    audit = stderr;
    int samples = Search::sample_anchored_conformers(&protein,&ligand,&pocket);
    audit = nullptr;
    if (samples<=6) fprintf(stderr,"Initial score %g; ligand carbons %d; anchor %s.\n",initial,ligand.get_heavy_atom_count(),carbon->name);
    assert(samples>6);
    float final = ligand.conformer_binding(nearby,&pocket).summed();
    assert(std::isfinite(final) && final<=initial+0.001f);
    assert(final<initial-0.001f);
    assert(carbon->loc.get_3d_distance(anchor_position)<0.0001f);
    assert(fabs(carbon->distance_to(nitrogen)-imine_length)<0.0001f);
    assert(carbon->is_bonded_to(nitrogen)==2);
    for (size_t i=0; i<fixed.size(); ++i)
        assert(fixed[i]->loc.get_3d_distance(positions[i])<0.0001f);
    for (const Length& length : lengths)
        assert(fabs(length.a->distance_to(length.b)-length.distance)<0.0001f);
    assert(ligand.get_atom("O") == nullptr || ligand.get_atom("O") != oxygen);
    printf("Anchored sampling: %d conformers; score %g -> %g; attachment and receptor fixed.\n",samples,initial,final);
    delete water;
}
