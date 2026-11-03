#include <cassert>
#include <cstdio>
#include <fstream>
#include <sstream>
#include "../classes/search.h"

int main()
{
    Protein protein;
    FILE* input = fopen("pdbs/OR4/OR4D1.active.pdb", "rb");
    assert(input);
    protein.load_pdb(input);
    fclose(input);
    std::ifstream file("sdf/pi-23472.sdf");
    std::stringstream buffer;
    buffer << file.rdbuf();
    std::string data = buffer.str();
    Molecule ligand;
    ligand.from_sdf(&data[0]);
    LigandTarget primary, secondary;
    primary.single_atom = ligand.get_atom("O17");
    secondary.single_atom = ligand.get_atom("C1");
    assert(primary.single_atom && secondary.single_atom);
    BestBindingResult targets;
    targets.protein = &protein;
    targets.ligand = &ligand;
    targets.pri_res = protein.get_residue(80);
    targets.sec_res = protein.get_residue(103);
    targets.pri_tgt = &primary;
    targets.sec_tgt = &secondary;
    Point center;
    for (const char* position : {"3.33","3.36","3.37","4.56","5.43"})
        center = center.add(protein.get_residue_bw(position)->get_CA_location());
    center.multiply(0.2);
    ligand.recenter(center);
    Point size = ligand.get_bounding_box().size();
    size.multiply(0.5);
    AminoAcid* neighbors[SPHREACH_MAX];
    protein.get_residues_can_clash_ligand(neighbors,&ligand,center,size);
    Search::align_targets(&ligand,center,&targets,1,true);
    const float chosen = ligand.get_intermol_clashes((Molecule**)neighbors);
    float best = chosen;
    LocatedVector axis = (Vector)targets.pri_res->get_barycenter().subtract(center);
    axis.origin = center;
    for (int i=0; i<30; ++i)
    {
        ligand.rotate(axis,M_PI*2/30);
        best = fmin(best,ligand.get_intermol_clashes((Molecule**)neighbors));
    }
    printf("Chosen rotation clash %g; sampled minimum %g.\n",chosen,best);
    assert(chosen <= best + 0.01f + fabs(best)*0.001f);
}
