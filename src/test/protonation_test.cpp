#include <cassert>
#include <cstdio>
#include "../classes/aminoacid.h"
#include "../classes/protein.h"

static void check_state(bool present)
{
    AminoAcid histidine('H'), acid('D');
    Atom *bare = nullptr, *donor = nullptr, *hydrogen = nullptr;
    for (int i=0; i<histidine.get_atom_count(); ++i)
    {
        Atom* atom = histidine.get_atom(i);
        if (atom->is_backbone || atom->get_family() != PNICTOGEN) continue;
        Atom* attached = atom->is_bonded_to("H");
        if (attached) { donor = atom; hydrogen = attached; }
        else bare = atom;
    }
    assert(bare && donor && hydrogen);
    if (present)
    {
        histidine.add_atom("H", "extra", bare, 1);
        bare->increment_charge(1);
    }
    Point direction = hydrogen->loc.subtract(donor->loc);
    direction.scale(1.8);
    Point target = hydrogen->loc.add(direction);
    acid.aamove(target.subtract(acid.get_atom("OD1")->loc));
    Molecule* nearby[] = {present ? &acid : nullptr, nullptr};
    const int count = histidine.get_atom_count();
    for (int i=0; i<3; ++i)
    {
        histidine.set_conditional_basicity(nearby);
        assert(histidine.get_atom_count() == count);
        assert(bare->num_bonded_to("H") == (present ? 1 : 0));
        assert(bare->get_charge() == (present ? 1 : 0));
    }
}

static void check_steric_protonation(bool blocked, bool internal)
{
    AminoAcid histidine('H', nullptr, false);
    Molecule acid("Acid oxygen");
    histidine.hisflips = nullptr;
    Atom* donor = histidine.get_atom("NE2");
    Atom* hydrogen = histidine.get_atom("HE2");
    Atom* bare = histidine.get_atom("ND1");
    assert(donor && hydrogen && bare);
    Point direction = hydrogen->loc.subtract(donor->loc);
    direction.scale(1.8);
    Point target = hydrogen->loc.add(direction);
    // Isolate the activating oxygen from alternative tautomer contacts.
    acid.add_atom("O", "O", &target, nullptr, 0, -1)->residue = 100;

    Molecule obstacle("Steric obstacle");
    if (blocked)
    {
        Atom test_hydrogen("H");
        Vector vector = bare->get_next_free_geometry(1);
        vector.r = InteratomicForce::covalent_bond_radius(&test_hydrogen, bare, 1);
        Point position = Point(vector).add(bare->loc);
        position.x += 0.5;
        if (internal) histidine.get_atom("HA")->move(position);
        else
        {
            assert(obstacle.from_smiles("C"));
            obstacle.hydrogenate();
            Atom* obstacle_hydrogen = nullptr;
            for (int i = 0; i < obstacle.get_atom_count(); ++i)
                if (obstacle.get_atom(i)->Z == 1) obstacle_hydrogen = obstacle.get_atom(i);
            assert(obstacle_hydrogen);
            obstacle.move(position.subtract(obstacle_hydrogen->loc));
        }
    }
    Molecule* nearby[] = {&acid, &obstacle, nullptr};
    const int count = histidine.get_atom_count();
    const auto histidine_motion = histidine.movability;
    const auto acid_motion = acid.movability;
    for (int repeat = 0; repeat < 3; ++repeat)
    {
        histidine.set_conditional_basicity(nearby);
        assert(histidine.get_atom_count() == count + (blocked ? 0 : 1));
        assert(bare->num_bonded_to("H") == (blocked ? 0 : 1));
        assert(bare->get_charge() == (blocked ? 0 : 1));
        if (blocked)
        {
            assert(histidine.movability == histidine_motion);
            assert(acid.movability == acid_motion);
        }
    }
}

static void check_neutral_protein_obstacle(bool blocked)
{
    Protein protein;
    assert(protein.add_residue(1, 'H'));
    assert(protein.add_residue(2, 'A'));
    AminoAcid* histidine = protein.get_residue(1);
    AminoAcid* neutral = protein.get_residue(2);
    assert(histidine && neutral && histidine->get_parent() == &protein);
    histidine->hisflips = nullptr;
    Atom *donor = nullptr, *hydrogen = nullptr, *bare = nullptr;
    for (int i = 0; i < histidine->get_atom_count(); ++i)
    {
        Atom* atom = histidine->get_atom(i);
        if (atom->is_backbone || atom->get_family() != PNICTOGEN) continue;
        Atom* attached = atom->is_bonded_to("H");
        if (attached) { donor = atom; hydrogen = attached; }
        else bare = atom;
    }
    assert(donor && hydrogen && bare);
    Molecule acid("Acid oxygen");
    Point direction = hydrogen->loc.subtract(donor->loc);
    direction.scale(1.8);
    Point target = hydrogen->loc.add(direction);
    acid.add_atom("O", "O", &target, nullptr, 0, -1)->residue = 100;
    Molecule* acids[] = {&acid, nullptr};
    const int count = histidine->get_atom_count();
    Atom test_hydrogen("H");
    Vector vector = bare->get_next_free_geometry(1);
    vector.r = InteratomicForce::covalent_bond_radius(&test_hydrogen, bare, 1);
    Point obstacle = Point(vector).add(bare->loc);
    obstacle.x += .5;
    if (!blocked) obstacle.x += 100;
    neutral->get_atom("HA")->move(obstacle);
    const auto motion = histidine->movability;
    histidine->set_conditional_basicity(acids);
    assert(histidine->get_atom_count() == count + (blocked ? 0 : 1));
    assert(bare->num_bonded_to("H") == (blocked ? 0 : 1));
    assert(bare->get_charge() == (blocked ? 0 : 1));
    if (blocked) assert(histidine->movability == motion);
}

static void check_flip_after_proton_removal()
{
    AminoAcid histidine('H', nullptr, false);
    assert(histidine.hisflips && histidine.hisflips[0]);
    HistidineFlip* flip = histidine.hisflips[0];
    Atom* original = flip->H;
    Atom* donor = original->get_bond_by_idx(0)->atom2;
    Atom* other = donor == flip->N1 ? flip->N2 : flip->N1;
    const int count = histidine.get_atom_count();

    histidine.delete_atom(original);
    histidine.do_histidine_flip(flip);
    assert(histidine.get_atom_count() == count - 1);
    assert(!donor->is_bonded_to("H") && !other->is_bonded_to("H"));

    Atom* replacement = histidine.add_atom("H", "replacement", other, 1);
    for (int repeat = 0; repeat < 4; ++repeat)
    {
        Atom* before = replacement->get_bond_by_idx(0)->atom2;
        histidine.do_histidine_flip(flip);
        Atom* after = replacement->get_bond_by_idx(0)->atom2;
        assert(after && after != before);
        assert(after == donor || after == other);
        assert(flip->H == replacement);
        assert(after->is_bonded_to(replacement));
        assert(histidine.get_atom_count() == count);
    }
}

int main()
{
    check_flip_after_proton_removal();
    check_state(false);
    check_state(true);
    check_steric_protonation(false, false);
    check_steric_protonation(true, false);
    check_steric_protonation(true, true);
    check_neutral_protein_obstacle(false);
    check_neutral_protein_obstacle(true);
    puts("Conditional protonation preserves state and excludes colliding hydrogens.");
}
