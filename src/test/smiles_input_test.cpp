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
