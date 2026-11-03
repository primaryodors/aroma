#include <cassert>
#include <cstdio>
#include "../classes/protein.h"

int main()
{
    Protein p;
    auto load = [&]() {
        FILE* f = fopen("pdbs/OR51/OR51E2.active.pdb", "r");
        assert(f);
        p.load_pdb(f);
        fclose(f);
    };
    load();
    Atom* first = p.get_residue(p.get_start_resno())->get_atom(0);
    assert(p.get_nearest_atom(first->loc)->distance_to(first) < 0.001);
    for (int r=p.get_start_resno(); r<=p.get_end_resno(); r++)
    {
        AminoAcid* aa = p.get_residue(r);
        if (!aa) continue;
        for (int j=0; j<aa->get_atom_count(); j++)
        {
            Atom* a = aa->get_atom(j);
            a->move(a->loc.add(Point(21,0,0)));
        }
    }
    for (int r=p.get_start_resno(); r<=p.get_end_resno(); r++)
    {
        AminoAcid* aa = p.get_residue(r);
        if (!aa) continue;
        Atom* a = aa->get_atom(0);
        assert(p.get_nearest_atom(a->loc)->distance_to(a) < 0.001);
    }
    assert(p.get_nearest_atom(Point(1000,1000,1000)));
    load();
    first = p.get_residue(p.get_start_resno())->get_atom(0);
    assert(p.get_nearest_atom(first->loc)->distance_to(first) < 0.001);
    AminoAcid* residue = p.get_residue(150);
    Atom* sidechain = residue->get_CB();
    assert(sidechain);
    Point removed_position = sidechain->loc;
    int old_count = residue->get_atom_count();
    p.get_nearest_atom(removed_position);
    residue->delete_sidechain();
    assert(residue->get_atom_count() < old_count);
    assert(residue->get_atom(residue->get_atom_count()) == nullptr);
    assert(p.get_nearest_atom(removed_position));
    puts("Spatial queries passed.");
}
