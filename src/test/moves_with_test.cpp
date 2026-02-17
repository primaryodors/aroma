
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

    m.from_smiles("c1ccccc1CCCO");
    p.add_sequence("MRAFSTVKAACWQNGDAYEQ");

    Atom* buffer[1024];
    AminoAcid* aa;
    Atom *a1, *a2;
    Bond* b;
    int i;

    a1 = m.get_atom("C6");
    a2 = m.get_atom("C7");
    if (!a1 || !a2) cerr << "ATOM ERROR" << endl << flush;
    b = a1->get_bond_between(a2);
    if (!b) cerr << "BOND ERROR" << endl << flush;
    b->fetch_moves_with_atom2(buffer);
    cout << "Moves with C6-C7:";
    for (i=0; buffer[i]; i++) cout << " " << buffer[i]->name;
    cout << endl;
    b->fetch_moves_rigidly_with_atom2(buffer);
    cout << "Moves rigid with C6-C7:";
    for (i=0; buffer[i]; i++) cout << " " << buffer[i]->name;
    cout << endl;

    aa = p.get_residue(18);
    if (!aa) cerr << "RESIDUE ERROR" << endl << flush;
    a1 = aa->get_atom("CA");
    a2 = aa->get_atom("CB");
    b = a1->get_bond_between(a2);
    if (!b) cerr << "BOND ERROR" << endl << flush;
    b->fetch_moves_with_atom2(buffer);
    cout << "Moves with " << aa->get_name() << ":CA-" << aa->get_name() << ":CB:";
    for (i=0; buffer[i]; i++) cout << " " << buffer[i]->name;
    cout << endl;
    b->fetch_moves_rigidly_with_atom2(buffer);
    cout << "Moves rigid with " << aa->get_name() << ":CA-" << aa->get_name() << ":CB:";
    for (i=0; buffer[i]; i++) cout << " " << buffer[i]->name;
    cout << endl;

    return 0;
}