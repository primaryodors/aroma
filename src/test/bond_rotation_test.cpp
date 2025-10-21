
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "../classes/molecule.h"

using namespace std;

Molecule mol("TheMolecule");

int main(int argc, char** argv)
{
    if (argc > 1)
    {
        if (file_exists(argv[1]))
        {
            FILE* fp = fopen(argv[1], "rb");
            char buffer[1024];
            std::string sdfdat = "";
            while (!feof(fp))
            {
                char* got = fgets(buffer, 1022, fp);
                if (!got) break;
                sdfdat += (std::string)buffer;
            }
            fclose(fp);
            mol.from_sdf(sdfdat.c_str());
        }
        else
        {
            mol.from_smiles(argv[1], false);        // obabel is too slow for unit test express.
        }
    }
    else
        mol.from_smiles("c1ccccc1C=O", false);
    
    Bond** b = mol.get_rotatable_bonds();
    if (!b || !b[0])
    {
        cout << "No rotatable bonds." << endl;
        return 0;
    }

    int i, n=0;
    for (i=0; b[i]; i++)
    {
        if (mol.atom_idx_from_ptr(b[i]->atom2) < mol.atom_idx_from_ptr(b[i]->atom1)) b[i] = b[i]->get_reversed();
        if (b[i]->can_rotate)
        {
            cout << "Bond " << *b[i] << " can rotate." << endl;
            n++;
        }
        else if (b[i]->can_flip)
        {
            cout << "Bond " << *b[i] << " can flip." << endl;
            n++;
        }
    }

    if (!n)
    {
        cout << "No rotatable bonds." << endl;
        return 0;
    }

    int l = 2;
    bool save = false;

    if (argc > l)
    {
        if (!strcmp(argv[l], "minc"))
        {
            mol.minimize_internal_clashes();
            l++;
            save = true;
        }
    }

    while (argc > l+2)
    {
        Atom* a = mol.get_atom(argv[l++]);
        Atom* b = mol.get_atom(argv[l++]);
        if (!a)
        {
            cerr << "No such atom " << argv[2] << "." << endl;
            return -1;
        }
        if (!b)
        {
            cerr << "No such atom " << argv[3] << "." << endl;
            return -1;
        }
        Bond* ab = a->get_bond_between(b);
        if (!ab)
        {
            cerr << "No bond between " << a->name << " and " << b->name << "." << endl;
            return -1;
        }
        if (ab->can_rotate)
        {
            ab->rotate(atof(argv[l++])*fiftyseventh, true, true);
            save = true;
        }
        else if (ab->can_flip)
        {
            ab->rotate(M_PI);
            save = true;
        }
        else
        {
            cerr << "Bond between " << a->name << " and " << b->name << " cannot rotate." << endl;
            return -1;
        }
    }

    if (save)
    {
        FILE* fp = fopen("rotated.sdf", "wb");
        if (!fp)
        {
            cerr << "FAILED to save output SDF." << endl;
            return -1;
        }

        mol.save_sdf(fp);
        fclose(fp);
        cout << "Wrote output file." << endl;
    }

    return 0;
}

