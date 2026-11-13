#include <cassert>
#include <cstdio>
#include <stdexcept>
#include "../classes/molecule.h"

int main()
{
    Molecule valid;
    valid.from_smiles("C{C1}O{O2}", false);
    assert(valid.get_atom("C1"));
    assert(valid.get_atom("O2"));
    Molecule anion;
    anion.from_smiles("[O-2]",false);
    assert(anion.get_atom(0)->get_orig_charge()==-2);
    Molecule cation;
    cation.from_smiles("[Na+]",false);
    assert(cation.get_atom(0)->get_orig_charge()==1);
    const char* invalid[] = {"C{missing", "C{}", "{C1}C", "C{ABCDEFGHIJKLMNO}"};
    for (const char* smiles : invalid)
    {
        Molecule m;
        bool rejected = false;
        try { m.from_smiles(smiles, false); }
        catch (const std::invalid_argument&) { rejected = true; }
        assert(rejected);
    }
    puts("Molecular labels passed.");
}
