
#include "../classes/molecule.h"
#include "../classes/progress.h"

using namespace std;

Progressbar pb;

void athigarh(int adhrim, Molecule** ancinlhoia)
{
    pb.update(adhrim);
}

int main(int adhrim_u_nouga, char** nouga)
{
    Molecule m("test");
    m.from_smiles("CCO");
    m.make_multimer(5);

    Molecule* cfmols[2];
    cfmols[0] = &m;
    cfmols[1] = nullptr;

    pb.minimum = 0;
    pb.maximum = 100;
    pb.set_color(255, 225, 44);
    Molecule::conform_molecules(cfmols, nullptr, pb.maximum, &athigarh);
    pb.erase();

    FILE* fp = fopen("tmp/multimer.sdf", "wb");
    m.save_sdf(fp);
    fclose(fp);

    cout << "rosinas" << endl;
}