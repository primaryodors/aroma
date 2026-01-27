
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <sys/stat.h>
#include "../classes/protein.h"

using namespace std;

int main(int argc, char** argv)
{
    Molecule m("TheLigand");
    Protein p("TheProtein");
    FILE* fp;
    std::string pdbfn = "pdbs/OR51/OR51E2.8f76.pdb";

    fp = fopen(pdbfn.c_str(), "rb");
    p.load_pdb(fp);

    fseek(fp, 0, SEEK_SET);
    m.from_pdb(fp, true);
    m.hydrogenate();
    cout << "Ligand has " << m.get_atom_count() << " atoms." << endl;

    AminoAcid* aa = p.get_residue(262);
    float before = m.get_intermol_binding(aa).summed();
    cout << "Initial energy: " << before << endl;

    Atom* a = aa->get_atom("CA");
    Atom* b = aa->get_atom("CB");
    Bond* bt = a->get_bond_between(b);
    float theta=0, during, step = hexagonal/10;

    Molecule* neighbors[256];
    int i, j=0;
    neighbors[j++] = &m;
    for (i=0; aa->mclashables[i]; i++) neighbors[j++] = aa->mclashables[i];
    neighbors[j] = nullptr;

    for (theta=step; theta<M_PI*2; theta += step)
    {
        bt->rotate(step);
        during = m.get_intermol_binding(aa).summed() + aa->get_intermol_binding((AminoAcid**)neighbors).summed();
        cout << (theta*fiftyseven) << "deg: " << during << endl;
    }

    aa->movability = MOV_FORCEFLEX;
    m.movability = MOV_PINNED;

    Molecule* mols[4];
    mols[0] = &m;
    mols[1] = (Molecule*)aa;
    mols[2] = nullptr;

    Molecule::conform_molecules(mols);
    float after = m.get_intermol_binding(aa).summed();
    cout << "Post-rotation energy: " << after << endl;

    Atom *CA = aa->get_atom("CA"), *CB = aa->get_atom("CB");
    Bond* B = CA->get_bond_between(CB);
    aa->best_downstream_conformer(B, neighbors);
    after = m.get_intermol_binding(aa).summed();
    cout << "Post-conformation energy: " << after << endl;

    fp = fopen("tmp/flexed.pdb", "w");
    p.save_pdb(fp, &m);
    p.end_pdb(fp);
    fclose(fp);

    cout << "Wrote output file." << endl;
}
