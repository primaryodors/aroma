#include <cassert>
#include <cstdio>
#include "../classes/scoring.h"

int main()
{
    Protein p;
    FILE* f = fopen("pdbs/OR51/OR51E2.active.pdb", "r");
    assert(f);
    p.load_pdb(f);
    fclose(f);
    Molecule ligand, water;
    ligand.from_smiles("CCO", false);
    water.from_smiles("O", false);
    Molecule* waters[] = {&water, nullptr};
    Atom* receptor_atom = p.get_residue(p.get_start_resno())->get_atom(0);
    Atom* ligand_atom = ligand.get_atom(0);
    Atom* water_atom = water.get_atom(0);
    Point receptor_pos = receptor_atom->loc;
    Point ligand_pos = ligand_atom->loc;
    Point water_pos = water_atom->loc;
    DockResult first;
    assert(!first.has_coordinates());
    first.capture_coordinates(nullptr, &ligand);
    assert(!first.has_coordinates());
    first.capture_coordinates(&p, nullptr);
    assert(!first.has_coordinates());
    first.capture_coordinates(&p, &ligand, waters);
    assert(first.has_coordinates());
    receptor_atom->move(receptor_pos.add(Point(1,0,0)));
    ligand_atom->move(ligand_pos.add(Point(2,0,0)));
    water_atom->move(water_pos.add(Point(3,0,0)));
    DockResult second;
    second.capture_coordinates(&p, &ligand, waters);
    first.disqualified = second.disqualified = true;
    first.kJmol = 50;
    second.kJmol = 100;
    assert(first.is_better_fallback_than(&second));
    assert(!second.is_better_fallback_than(&first));
    first.kJmol = 0;
    assert(first.is_better_fallback_than(&second));
    second.kJmol = nan("unset");
    assert(first.is_better_fallback_than(&second));
    assert(!second.is_better_fallback_than(&first));
    DockResult empty;
    assert(!empty.is_better_fallback_than(&first));
    assert(first.is_better_fallback_than(&empty));
    assert(second.is_better_fallback_than(nullptr));
    waters[0] = nullptr;
    first.restore_coordinates();
    assert(first.water_molecules().size() == 1);
    assert(receptor_atom->loc.get_3d_distance(receptor_pos) < 0.001);
    assert(ligand_atom->loc.get_3d_distance(ligand_pos) < 0.001);
    assert(water_atom->loc.get_3d_distance(water_pos) < 0.001);
    second.restore_coordinates();
    assert(receptor_atom->loc.get_3d_distance(receptor_pos.add(Point(1,0,0))) < 0.001);
    assert(ligand_atom->loc.get_3d_distance(ligand_pos.add(Point(2,0,0))) < 0.001);
    assert(water_atom->loc.get_3d_distance(water_pos.add(Point(3,0,0))) < 0.001);
    puts("Result coordinates passed.");
}
