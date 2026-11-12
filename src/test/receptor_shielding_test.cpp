#include "../classes/scoring.h"
#include <cassert>
#include <cstdio>
#include <cstring>

static float hbond_to_aspartate(Molecule& ligand, AminoAcid* residue, Protein* protein)
{
    for (float& value : total_binding_by_type) value = 0;
    compute_interall = true;
    ninterall = 0;
    ligand.get_intermol_binding(residue,false,false,nullptr,protein);
    float result = total_binding_by_type[hbond];
    compute_interall = false;
    return result;
}

int main()
{
    const char* path = "test/data/OR1A1_ester_shielding.pdb";
    FILE* input = fopen(path,"r");
    assert(input);
    Protein protein;
    protein.load_pdb(input);
    rewind(input);
    Molecule ligand;
    ligand.from_pdb(input,true);
    fclose(input);
    AminoAcid* aspartate = protein.get_residue(111);
    assert(aspartate && ligand.get_atom("O9"));
    Atom* oxygen = ligand.get_atom("O9");
    Atom* hydrogen = aspartate->get_atom("HN");
    assert(hydrogen);
    assert(!ligand.shielded(oxygen,hydrogen));
    assert(!aspartate->shielded(oxygen,hydrogen));
    assert(protein.get_residue(108)->shielded(oxygen,hydrogen)
        || protein.get_residue(109)->shielded(oxygen,hydrogen));
    assert(protein.shielded(oxygen,hydrogen));
    float pair_only = hbond_to_aspartate(ligand,aspartate,nullptr);
    float contextual = hbond_to_aspartate(ligand,aspartate,&protein);
    printf("Asp111 H-bond: pair-only %g, receptor-shielded %g\n",pair_only,contextual);
    assert(pair_only < -0.1f);
    assert(contextual > pair_only+0.1f);
    AminoAcid* selected[] = {aspartate,nullptr};
    DockResult result(&protein,&ligand,1,selected);
    assert(result.bytype[hbond] > pair_only+0.1f);

    // Moving blockers refreshes the atom index and restores the open interaction.
    for (int number : {108,109,110})
    {
        AminoAcid* residue = protein.get_residue(number);
        assert(residue);
        for (int i=0; i<residue->get_atom_count(); i++)
        {
            Atom* atom = residue->get_atom(i);
            atom->move(atom->loc.add(Point(100,0,0)));
        }
    }
    assert(!protein.shielded(oxygen,hydrogen));
    float exposed = hbond_to_aspartate(ligand,aspartate,&protein);
    assert(exposed < -0.1f);
    puts("Receptor shielding, score accounting, and moved-blocker checks passed.");
}
