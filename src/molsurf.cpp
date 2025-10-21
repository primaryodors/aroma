
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstring>
#include "classes/molecule.h"

using namespace std;

int main(int argc, char** argv)
{
    int i;
    char* SMILES = new char[2];
    strcpy(SMILES, "O");
    for (i=1; i<argc; i++)
    {
        if (argv[i][0] == '-')
        {
            //
        }
        else SMILES = argv[i];
    }

    Molecule mol("mol");
    mol.from_smiles(SMILES);
    if (!global_water.get_atom_count())
    {
        global_water.from_smiles("O{O1}");
        global_water.obtain_vdW_surface(vdw_surface_density);
    }
    Atom* polarH = global_water.get_atom(1);
    Atom* polarO = global_water.get_atom(0);

    mol.obtain_vdW_surface(vdw_surface_density);

    int n = mol.get_atom_count(), nv = mol.get_vdW_vertex_count();
    float nvrat = (float)nv / global_water.get_vdW_vertex_count();
    float wsp = 0;          // water solvation potential
    for (i=0; i<n; i++)
    {
        Atom* a = mol.get_atom(i);
        if (!a) continue;

        cout << a->name << " ";
        float f = mol.get_atom_vdW_vertex_count(a);
        f /= nv;
        cout << (f*100) << "%" << endl;

        float eH = InteratomicForce::potential_binding(a, polarH, false);
        float eO = InteratomicForce::potential_binding(a, polarO, false);

        wsp -= f*(water_surface_H*eH + (1.0-water_surface_H)*eO)*9*nvrat;
    }

    cout << "Potential interaction energy in water: " << wsp << " kJ/mol." << endl;
    cout << "Surface area: " << mol.get_surface_area() << endl;
    cout << "Molecular weight: " << mol.get_molecular_wt() << "." << endl;
    return 0;
}