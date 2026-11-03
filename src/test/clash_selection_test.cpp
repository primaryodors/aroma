#include <cassert>
#include <cstdio>
#include "../classes/protein.h"

int main()
{
    Protein p;
    FILE* f = fopen("pdbs/OR51/OR51E2.active.pdb", "r");
    assert(f);
    p.load_pdb(f);
    fclose(f);
    AminoAcid* last = p.get_residue(p.get_end_resno());
    assert(last);
    last->priority = true;
    last->renumber(9999);
    AminoAcid* selected[SPHREACH_MAX+4];
    int count = p.get_residues_can_clash_ligand(selected, nullptr, last->get_CA_location(), Point(8,8,8), nullptr);
    bool found = false;
    for (int i=0; i<count; i++) if (selected[i] == last) found = true;
    assert(found);
    assert(selected[count] == nullptr);
    int excluded[] = {9999, 0};
    count = p.get_residues_can_clash_ligand(selected, nullptr, last->get_CA_location(), Point(8,8,8), nullptr, false, excluded);
    for (int i=0; i<count; i++) assert(selected[i] != last);
    puts("Clash selection passed.");
}
