#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "../classes/protein.h"

using namespace std;

Molecule* mols[4];

int main(int argc, char** argv)
{
    mols[0] = new Molecule("benzaldehyde", "c1ccccc1C=O");
    mols[1] = new Molecule("butylamine", "CCCC[NH3+]");
    mols[2] = mols[3] = nullptr;

    mols[1]->create_Schiff_base(mols[0]);
    cout << "Schiff enthalpy of formation: " << intermol_covalent_enthalpy << endl;

    FILE* fp = fopen("schiff_test.sdf", "w");
    if (!fp)
    {
        cerr << "FAILED to open SDF output file for writing." << endl;
        return 1;
    }
    mols[0]->save_sdf(fp, &mols[1]);
    fclose(fp);
    cout << "Wrote output SDF." << endl;

    Protein p("test");
    p.add_sequence("TIVKRTCQALTILTA");
    p.make_helix(1, p.get_end_resno(), ALPHA_PHI, ALPHA_PSI);
    mols[2] = new Molecule("aldehydeMNA", "O=CC(C)CCCCCCCCC");
    AminoAcid* aa = p.get_residue(4);
    aa->create_Schiff_base(mols[2]);

    Atom *NZ = aa->get_atom("NZ");
    Bond *b = NZ->get_bond_by_idx(2);
    if (!b)
    {
        cerr << "Schiff bond is missing." << endl;
        return 2;
    }
    Atom *C = b->atom2;

    cout << "Schiff atoms are: " << NZ->residue << ":" << NZ->name << ", " << C->residue << ":" << C->name << endl;
    cout << "Bond distance: " << NZ->distance_to(C) << endl;

    Bond **bb = ((Molecule*)aa)->get_rotatable_bonds();
    int i;
    if (!bb)
    {
        cerr << "No rotatable bonds." << endl;
        return 2;
    }
    for (i=0; bb[i]; i++)
    {
        cout << "Rotatable: " << bb[i]->atom1->name << "-" << bb[i]->atom2->name << endl;
    }

    for (i=0; bb[i]; i++)
    {
        bb[i]->rotate(M_PI);
        cout << "Rotated " << bb[i]->atom1->name << "-" << bb[i]->atom2->name
            << " Schiff bond distance: " << NZ->distance_to(C) << endl;
        bb[i]->rotate(M_PI);
    }

    fp = fopen("schiff_test.pdb", "w");
    if (!fp)
    {
        cerr << "FAILED to open PDB output file for writing." << endl;
        return 1;
    }
    p.save_pdb(fp, mols[2]);
    p.end_pdb(fp);
    fclose(fp);
    cout << "Wrote output PDB." << endl;

    delete mols[2];
    mols[2] = new Molecule("isobutyraldehyde", "O=CC(C)C");
    aa = p.get_residue(5);
    aa->identify_Schiff_amine(&NZ, nullptr, nullptr);
    mols[2]->identify_Schiff_carbonyl(&C, nullptr);
    aa->create_Schiff_base(mols[2]);

    b = NZ->get_bond_between(C);
    if (!b)
    {
        cout << "No Schiff bond." << endl;
    }
    else
    {
        cout << "Schiff atoms are: " << NZ->residue << ":" << NZ->name << ", " << C->residue << ":" << C->name << endl;
        cout << "Bond distance: " << NZ->distance_to(C) << endl;
    }

    fp = fopen("schiff_test_1.pdb", "w");
    if (!fp)
    {
        cerr << "FAILED to open PDB output file for writing." << endl;
        return 1;
    }
    p.save_pdb(fp, mols[2]);
    p.end_pdb(fp);
    fclose(fp);
    cout << "Wrote output PDB." << endl;

    aa = p.get_residue(8);
    delete mols[0];
    delete mols[1];
    delete mols[2];
    mols[0] = new Molecule("cadaverine", "[NH3+]CCCCC[NH3+]");
    mols[1] = new Molecule("octanal", "CCCCCCCC=O");
    mols[2] = nullptr;

    if (aa->create_Schiff_base(mols[0]))
    {
        cerr << "ERROR: Amine formed Schiff base with amide." << endl << flush;
        return 1;
    }
    if (aa->create_Schiff_base(mols[1]))
    {
        cerr << "ERROR: Aldehyde formed Schiff base with amide." << endl << flush;
        return 1;
    }

    delete mols[0];
    delete mols[1];
    mols[0] = new Molecule("methyl_anthranilate", "Nc1ccccc1C(=O)OC");
    mols[1] = new Molecule("hydroxycitronellal", "CC(CCCC(C)(C)O)CC=O");
    mols[2] = nullptr;

    mols[1]->create_Schiff_base(mols[0]);
    cout << "Schiff enthalpy of formation: " << intermol_covalent_enthalpy << endl;

    fp = fopen("schiff_test_1.sdf", "w");
    if (!fp)
    {
        cerr << "FAILED to open SDF output file for writing." << endl;
        return 1;
    }
    mols[0]->save_sdf(fp, &mols[1]);
    fclose(fp);
    cout << "Wrote output SDF." << endl;

    delete mols[0];
    delete mols[1];
    mols[0] = new Molecule("arginine_acetylated", "CC(=O)N[C@H](C(=O)O)CCCNC(N)=[NH2+]");
    mols[1] = new Molecule("hydroxycitronellal", "CC(CCCC(C)(C)O)CC=O");
    mols[2] = nullptr;

    mols[1]->create_Schiff_base(mols[0]);
    cout << "Schiff enthalpy of formation: " << intermol_covalent_enthalpy << endl;

    fp = fopen("schiff_test_2.sdf", "w");
    if (!fp)
    {
        cerr << "FAILED to open SDF output file for writing." << endl;
        return 1;
    }
    mols[0]->save_sdf(fp, &mols[1]);
    fclose(fp);
    cout << "Wrote output SDF." << endl;

    return 0;
}