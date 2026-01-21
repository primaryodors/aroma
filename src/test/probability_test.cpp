#include <cstring>
#include <iostream>
#include <stdio.h>
#include <string>
#include <math.h>
#include "../classes/molecule.h"

using namespace std;

int main(int argc, char** argv)
{
    Interaction e, eref = 0;
    float f;

    for (f=-25.1; f<=25.1; f+=0.5)
    {
        e = f;
        float probs = e.probability(eref);
        cout << e.summed() << ": " << probs << endl;
    }
}