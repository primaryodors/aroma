
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstring>
#include "../classes/protein.h"

using namespace std;

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        cerr << "No input file." << endl;
        return -1;
    }

    Protein p(argv[1]);
    FILE* pf = fopen(argv[1], "r");
    if (!pf)
    {
        cerr << "Error trying to read " << argv[1] << endl;
        return 0xbadf12e;
    }
    p.load_pdb(pf);
    fclose(pf);

    if (!p.get_region_end(7))
    {
        cerr << "GPCR regions not set in protein." << endl;
        return 0x9076bc12;
    }

    int i;
    AminoAcid *aa;
    float phi, psi, bend;
    int start = p.get_start_resno(), end = p.get_region_start(1);
    float stretch = p.max_stretch(start, end);

    cout << "Head region can stretch " << stretch << " Å" << endl;

    start = end;
    end = p.get_region_end(1);

    // TODO:

    start = end;
    end = p.get_region_start(2);
    stretch = p.max_stretch(start, end);
    cout << "CYT1 can stretch " << stretch << " Å" << endl;

    start = end;
    end = p.get_region_end(2);

    // TODO:

    start = end;
    end = p.get_region_start(3);
    stretch = p.max_stretch(start, end);
    cout << "EXR1 can stretch " << stretch << " Å" << endl;

    start = end;
    end = p.get_region_end(3);

    // TODO:

    start = end;
    end = p.get_region_start(4);
    stretch = p.max_stretch(start, end);
    cout << "CYT2 can stretch " << stretch << " Å" << endl;

    start = end;
    end = p.get_region_end(4);

    // TODO:

    start = end;
    aa = p.get_residue_bw(45, 50);
    end = aa->get_residue_no();
    stretch = p.max_stretch(start, end);
    cout << "EXR2 side IV can stretch " << stretch << " Å" << endl;

    start = end;
    aa = p.get_residue_bw(45, 58);
    end = aa->get_residue_no();

    // TODO:

    start = end;
    end = p.get_region_start(5);
    stretch = p.max_stretch(start, end);
    cout << "EXR2 side V can stretch " << stretch << " Å" << endl;

    start = end;
    end = p.get_region_end(5);

    // TODO:

    start = end;
    end = p.get_region_start(6);
    stretch = p.max_stretch(start, end);
    cout << "CYT3 can stretch " << stretch << " Å" << endl;

    start = end;
    end = p.get_region_end(6);

    cout << endl;
    for (i=start; i<=end; i++)
    {
        aa = p.get_residue(i);
        if (aa)
        {
            phi = aa->get_phi();
            psi = aa->get_psi();

            bend = (ALPHA_PHI - phi + ALPHA_PSI - psi); 
            bend /= (M_PI*2);
            bend -= round(bend);
            bend *= (M_PI*2);

            cout << aa->get_name() << " helix bend = " << (bend*fiftyseven) << endl;
        }
    }
    cout << endl;

    start = end;
    end = p.get_region_start(7);
    stretch = p.max_stretch(start, end);
    cout << "EXR3 can stretch " << stretch << " Å" << endl;
}