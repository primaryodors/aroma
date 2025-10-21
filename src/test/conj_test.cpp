
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "../classes/molecule.h"
#include "../classes/conj.h"
#include "../classes/aminoacid.h"

using namespace std;

int main(int argc, char** argv)
{
    Molecule m("testmol");
    // m.from_smiles("N[C@H](C(=O)O)CCCNC(N)=[NH2+]", false);
    FILE* fp = fopen("sdf/arginine.sdf", "rb");
    char buffer[8192];
    int gfys = fread(buffer, 1, 8190, fp);
    fclose(fp);
    m.from_sdf(buffer);

    int i, n = m.get_atom_count();
    for (i=0; i<n; i++)
    {
        Atom* a = m.get_atom(i);
        if (!a) continue;
        cout << a->name << " has charge " << a->get_orig_charge() << ", effectively "
            << a->get_charge() << " or " << a->is_conjugated_to_charge();
        if (a->conjugation)
        {
            cout << " from conjugation of " << a->conjugation->count_atoms() << " atoms with net charge " << a->conjugation->get_net_charge();
        }
        cout << "." << endl;
    }

    AminoAcid Arg('R');
    n = Arg.get_atom_count();
    for (i=0; i<n; i++)
    {
        Atom* a = Arg.get_atom(i);
        if (!a) continue;
        cout << a->name << " has charge " << a->get_orig_charge() << ", effectively "
            << a->get_charge() << " or " << a->is_conjugated_to_charge();
        if (a->conjugation)
        {
            cout << " from conjugation of " << a->conjugation->count_atoms() << " atoms with net charge " << a->conjugation->get_net_charge();
        }
        cout << "." << endl;
    }
}