#include <cstring>
#include <iostream>
#include <stdio.h>
#include <string>
#include <math.h>
#include "../classes/molecule.h"

using namespace std;

int main(int argc, char** argv)
{
    Molecule m("Test");

    m.from_smiles("[SH-]");
    m.make_multimer(2);

    Molecule* mon1 = m.get_monomer(0);
    Molecule* mon2 = m.get_monomer(1);

    if (!mon1 || !mon2)
    {
        cout << "Multiplicity failed." << endl;
        throw 0xbadc0de;
    }

    Atom* S1 = mon1->get_atom(0)->get_heavy_atom();
    Atom* S2 = mon2->get_atom(0)->get_heavy_atom();
    Atom* H1 = S1->get_bond_by_idx(0)->atom2;
    Atom* H2 = S2->get_bond_by_idx(0)->atom2;

    if (mon1 == mon2 || S1 == S2 || H1 == H2 || H1 == S1 || H2 == S2)
    {
        cout << "Multiplicity error." << endl;
        throw 0xbadc0de;
    }

    Rotation rot = align_points_3d(H1->loc, S1->loc.add(S1->loc.subtract(S2->loc)), S1->loc);
    mon1->rotate(rot);
    rot = align_points_3d(H2->loc, S2->loc.add(S2->loc.subtract(S1->loc)), S2->loc);
    mon2->rotate(rot);
    SCoord mov = S1->loc.subtract(S2->loc);
    mov.r -= 4;
    mon2->move(mov);

    cout << m.get_intermol_binding(&m) << endl;

    FILE* fp = fopen("tmp/inte.sdf", "wb");
    m.save_sdf(fp);
    fclose(fp);

    return 0;

    Molecule m1;
    m1.from_smiles("CCCCO");
    m1.make_multimer(2);
    mon1 = m1.get_monomer(0);
    mon2 = m1.get_monomer(1);
    Atom *a1, *a2;

    /*mon1->mutual_closest_atoms(mon2, &a1, &a2);
    mov = a1->loc.subtract(a2->loc);
    mov.r -= (a1->vdW_radius+a2->vdW_radius);
    mon2->move(mov);*/

    m1.get_atom("C1")->move(Point(1.855, 7.408, 6.349));
    m1.get_atom("C2")->move(Point(2.257, 7.372, 7.813));
    m1.get_atom("C3")->move(Point(3.108, 8.583, 8.178));
    m1.get_atom("C4")->move(Point(4.431, 8.161, 8.800));
    m1.get_atom("O5")->move(Point(5.076, 9.287, 9.379));
    m1.get_atom("H6")->move(Point(0.774, 7.268, 6.247));
    m1.get_atom("H7")->move(Point(2.357, 6.614, 5.787));
    m1.get_atom("H8")->move(Point(2.118, 8.369, 5.895));
    m1.get_atom("H9")->move(Point(2.815, 6.450, 8.013));
    m1.get_atom("H10")->move(Point(1.355, 7.347, 8.436));
    m1.get_atom("H11")->move(Point(2.568, 9.233, 8.878));
    m1.get_atom("H12")->move(Point(3.300, 9.195, 7.288));
    m1.get_atom("H13")->move(Point(5.100, 7.722, 8.052));
    m1.get_atom("H14")->move(Point(4.269, 7.420, 9.591));
    m1.get_atom("H15")->move(Point(5.156, 9.119, 10.332));
    m1.get_atom("C16")->move(Point(2.635, 9.597, 4.349));
    m1.get_atom("C17")->move(Point(2.006, 10.924, 3.961));
    m1.get_atom("C18")->move(Point(2.789, 12.094, 4.546));
    m1.get_atom("C19")->move(Point(1.964, 12.850, 5.577));
    m1.get_atom("O20")->move(Point(1.592, 11.976, 6.634));
    m1.get_atom("H21")->move(Point(2.112, 8.770, 3.859));
    m1.get_atom("H22")->move(Point(2.581, 9.441, 5.431));
    m1.get_atom("H23")->move(Point(3.687, 9.561, 4.048));
    m1.get_atom("H24")->move(Point(0.970, 10.952, 4.319));
    m1.get_atom("H25")->move(Point(1.974, 11.004, 2.869));
    m1.get_atom("H26")->move(Point(3.097, 12.784, 3.750));
    m1.get_atom("H27")->move(Point(3.718, 11.737, 5.008));
    m1.get_atom("H28")->move(Point(1.054, 13.266, 5.132));
    m1.get_atom("H29")->move(Point(2.542, 13.674, 6.008));
    m1.get_atom("H30")->move(Point(2.133, 11.174, 6.555));

    cout << m1.get_intermol_binding(&m1) << endl;

    fp = fopen("tmp/inte1.sdf", "wb");
    m1.save_sdf(fp);
    fclose(fp);
}