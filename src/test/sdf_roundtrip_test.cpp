#include "../classes/molecule.h"
#include <cassert>
#include <cstdio>
#include <sstream>

static std::string exported(Molecule& molecule)
{
    FILE* file=tmpfile();assert(file);
    assert(molecule.save_sdf(file));rewind(file);
    std::string text;char buffer[1024];
    while (fgets(buffer,sizeof(buffer),file)) text+=buffer;
    fclose(file);return text;
}

int main()
{
    Molecule charged;
    charged.from_smiles("CC(=O)[O-]",false);
    assert(charged.get_atom(2)->get_orig_charge()==-1);
    std::string sdf=exported(charged);
    const auto end=sdf.find("M  END");assert(end!=std::string::npos);
    sdf.insert(end,"A    2\nC13\n");
    Molecule restored;
    assert(restored.from_sdf(sdf.c_str())==charged.get_atom_count());
    assert(restored.get_atom("C13")==restored.get_atom(1));
    assert(restored.get_atom(2)->get_orig_charge()==-1);
    std::string again=exported(restored);
    assert(again.find("C13\n")!=std::string::npos);
    Molecule final;
    assert(final.from_sdf(again.c_str())==charged.get_atom_count());
    assert(final.get_atom("C13"));
    assert(final.get_atom(2)->get_orig_charge()==-1);

    Molecule benzene;
    benzene.from_smiles("c1ccccc1",false);
    std::string aromatic=exported(benzene);
    assert(final.from_sdf(aromatic.c_str())==benzene.get_atom_count());
    assert(final.get_atom(0)->is_pi());
    Bond* bond=final.get_atom(0)->get_bond_between(final.get_atom(1));
    assert(bond && fabs(bond->cardinality-1.5f)<0.001f);
    const int count=final.get_atom_count();
    std::string truncated=aromatic.substr(0,aromatic.find("M  END"));
    assert(final.from_sdf(truncated.c_str())==0);
    assert(final.get_atom_count()==count);
    assert(final.from_sdf("not an SDF")==0);
    std::string legacy=sdf;
    const auto first_bond=legacy.find("  1  2  1");
    if (first_bond!=std::string::npos)
    {
        legacy.replace(first_bond,9,"  1  2 1.00");
        assert(final.from_sdf(legacy.c_str())==charged.get_atom_count());
    }
    puts("SDF aromatic bonds, charges, atom aliases, reuse, and truncated-input checks passed.");
}
