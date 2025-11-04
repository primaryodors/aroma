
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "../classes/molecule.h"

using namespace std;

Molecule* mols[3];

void iteration_callback(int iter, Molecule** mols)
{
    #if _dbg_mol_test
    Interaction e = mols[0]->get_intermol_binding(mols[1]);

    float nodist = 0;
    Atom* O = mols[0]->get_atom("O6");
    Atom* N = mols[1]->get_atom("N1");
    if (N && O) nodist += O->loc.get_3d_distance(N->loc);
    O = mols[1]->get_atom("O6");
    N = mols[0]->get_atom("N1");
    if (N && O) nodist += O->loc.get_3d_distance(N->loc);

    cout << e.summed() << "\t\t" << e.clash << "\t\t" << nodist << endl;
    #endif

    #if _dbg_mol_frames

    char buffer[256];
    sprintf(buffer, "tmp/frame%d.sdf", iter);

    FILE* pf = fopen(buffer, "wb");
    Molecule* ligands[2];
    ligands[0] = mols[1];
    ligands[1] = 0;
    mols[0]->save_sdf(pf, ligands);
    fclose(pf);

    #endif
}

int main(int argc, char** argv)
{
    // TODO: These values are set way too permissively.
    // The intermolecular code should be fine tuned to succeed
    // with threshold values of at worst 15, 3, and 2.5 * 2.
    float energyLevelThreshold = 10;
    float clash_limit = 8;
    float N_O_spacing = 4 * 2;

    Molecule m("nothing");
    cout << "# Created empty molecule named " << m.get_name() << ".\n";

    char buffer[65536];
    char tstname[1024] = "BZN.sdf";

    if (argc > 1) strcpy(tstname, argv[1]);

    Molecule m1(tstname);
    cout << "# Created molecule named " << m1.get_name() << ".\n";
    int nloaded;

    FILE* pf = fopen(tstname, "rb");
    if (!pf)
    {
        m1.from_smiles(argv[1]);
    }
    else
    {
        int got = fread(buffer, 1, 65535, pf);
        fclose(pf);
        cout << "Read data.\n";
        nloaded = m1.from_sdf(buffer);
    }

    cout << "# Loaded " << nloaded << " atoms of " << tstname << " into molecule.\n";

    int rc = m1.get_num_rings();
    if (rc) cout << "# Found " << rc << " ring(s)." << endl;

    int i;
    for (i=0; i<rc; i++)
    {
        cout << "# Ring atoms: ";
        Atom** ra = m1.get_ring_atoms(i);
        Atom::dump_array(ra);
        cout << endl;
        delete ra;

        bool cp = false;
        if (m1.ring_is_aromatic(i))
        {
            cout << "# Ring " << i << " is aromatic." << endl;
            cp = true;
        }
        else if (m1.ring_is_coplanar(i))
        {
            cout << "# Ring " << i << " is coplanar." << endl;
            cp = true;
        }

        Point rc = m1.get_ring_center(i);
        cout << "# Ring centered at [" << rc.x << "," << rc.y << "," << rc.z << "]." << endl;

        if (cp)
        {
            Vector normal = m1.get_ring_normal(i);
            cout << "Ring normal: φ=" << (normal.phi * 180.0/M_PI) << "° θ=" << (normal.theta * 180.0/M_PI) << "°." << endl;
        }
    }

    m1.minimize_internal_clashes();
    float ic1 = m1.get_internal_clashes();
    cout << "# Internal clashes: " << ic1 << " cu. A." << endl;
    if (ic1 > 0.01) cout << "Internal clashes greater than threshold. FAIL." << endl;

    char** anames = m1.get_atom_names();
    if (anames)
    {
        //cout << "Name of first atom " << anames[i] << endl << endl;
        for (i=0; anames[i]; i++)
        {
            Atom* a = m1.get_atom(anames[i]);
            if (a)
            {
                float ab = a->get_acidbase();
                if (ab) cout << "# Atom " << anames[i] << " has basicity " << ab << endl;
            }
        }
        delete anames;
    }
    else cout << "# No names.\n";

    Bond** b = m1.get_rotatable_bonds();
    if (b)
    {
        for (i=0; b[i]; i++)
        {
            cout << "# Bond " << b[i]->atom1->name << " - " << b[i]->atom2->name << " can rotate." << endl;
            b[i]->rotate(30.0 * M_PI / 180);
        }
    }

    Point pt1(0.2,-0.5,-0.3);
    m1.move(pt1);
    Point loc = m1.get_barycenter();
    // cout << "# Molecule moved to [" << loc.x << "," << loc.y << "," << loc.z << "]." << endl;

    Molecule m2("test");

    if (argc > 2)
    {
        m2.from_smiles(argv[2]);
        pt1.x = pt1.y = pt1.z = 0.5;
        while (m1.get_intermol_clashes(&m2) > 10)
        {
            m2.move(pt1);
            // cout << "# " << m2.get_barycenter() << endl;
        }
    }
    else
    {
        Atom* C1 = m2.add_atom("C", "C1", 0, 0);
        cout << "Added a carbon atom. Its location is " << C1->loc.printable()  << "." << endl;

        Atom* C2 = m2.add_atom("C", "C2", C1, 1);
        cout << "Added another carbon atom. Its location is " << C2->loc.printable() << "." << endl;

        Atom* O3 = m2.add_atom("O", "O3", C2, 1);
        cout << "Added an oxygen atom. Its location is " << O3->loc.printable() << "." << endl;

        m2.hydrogenate();
    }

    m2.minimize_internal_clashes();
    float im12 = m1.get_intermol_clashes(&m2);
    cout << "# Loaded test ligand. Intermol clashes: " << im12 << " cu. A." << endl;

    Vector v1(&pt1);
    float rotdeg = -30;
    m2.rotate(&v1, rotdeg * M_PI/180);
    cout << "# Rotated molecule 2 by " << rotdeg << " degrees. Intermol clashes: " << m1.get_intermol_clashes(&m2) << " cu. A." << endl;

    Point pt(0,0,1.0);
    Vector v(&pt);
    float ttlmv = 0;

    Interaction e = m1.get_intermol_binding(&m2);

    cout << "# Initial intermol clashes: " << e.clash << " cu. A." << endl;
    cout << "# Initial intermol energy level: " << e.summed() << " kJ/mol." << endl;

    m1.reset_conformer_momenta();
    m2.reset_conformer_momenta();
    mols[0] = &m1;
    mols[1] = &m2;
    mols[2] = NULL;
    Molecule::conform_molecules(mols, 200, &iteration_callback);
    Vector optimize = m1.motion_to_optimal_contact(&m2);
    cout << "# Optimization moves molecule 1 by " << optimize << endl;
    m1.move(optimize);
    #if _peratom_audit
    interauditing = true;
    #endif

    e = m1.get_intermol_binding(&m2);
    cout << "\n# Post-conformation intermol energy level: " << e.summed() << " kJ/mol." << endl;

    #if _peratom_audit
    cout << endl << "# Interatomic Audit:" << endl;
    int ian = interaudit.size(), iai;
    for (iai=0; iai<ian; iai++) cout << "# " << interaudit[iai] << endl;
    cout << endl << endl;
    interauditing = false;
    #endif

    float nodist = 0;
    Atom* O = m1.get_atom("O6");
    Atom* N = m2.get_atom("N1");
    if (N && O) nodist += O->loc.get_3d_distance(N->loc);
    O = m2.get_atom("O6");
    N = m1.get_atom("N1");
    if (N && O) nodist += O->loc.get_3d_distance(N->loc);

    if (e.summed() <= energyLevelThreshold && e.clash < clash_limit && nodist < N_O_spacing)
        cout << "Energy level below threshold, SUCCESS.\n";
    else
    {
        if (e.summed() > energyLevelThreshold) cout << "Energy level above threshold, FAIL.\n";
        else if (e.clash >= clash_limit) cout << "Intermolecular clashes " << e.clash << " are above threshold, FAIL.\n";
        else if (nodist >= N_O_spacing) cout << "Atoms are too far apart (" << nodist << "A). FAIL.\n";
    }

    const char* tstoutf = "output.sdf";
    pf = fopen(tstoutf, "wb");
    Molecule* ligands[2];
    ligands[0] = &m2;
    ligands[1] = 0;
    if (m1.save_sdf(pf, ligands)) cout << "Saved " << tstoutf << endl;
    else cout << "Failed to save " << tstoutf << endl;
    fclose(pf);

    const char* tstoutf1 = "output.pdb";
    pf = fopen(tstoutf1, "wb");
    m1.save_pdb(pf);
    fclose(pf);



    return 0;
}
