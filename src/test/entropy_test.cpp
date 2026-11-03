#include <cassert>
#include <cmath>
#include <cstdio>
#include "../classes/scoring.h"

int main()
{
    DockResult empty[2];
    DockResult::entropy(empty, 2, 1);
    assert(empty[0].estimated_TDeltaS == 0);
    DockResult::entropy(nullptr, 0, 0);
    Molecule m;
    m.from_smiles("CCO", false);
    DockResult results[48*2];
    for (int i=0; i<48; i++)
    {
        m.get_atom(0)->move(Point(i,0,0));
        results[i*2].capture_coordinates(nullptr, &m);
        results[i*2].kJmol = i % 2 ? -10 : 10;
    }
    DockResult::entropy(results, 48, 2);
    for (int i=0; i<48; i++)
    {
        assert(std::isfinite(results[i*2].estimated_TDeltaS));
        assert(results[i*2].estimated_TDeltaS < 0);
        assert(results[i*2+1].estimated_TDeltaS == 0);
        results[i*2].kJmol = 0;
    }
    DockResult::entropy(results, 48, 2);
    for (int i=0; i<48; i++) assert(std::isfinite(results[i*2].estimated_TDeltaS));
    puts("Entropy limits passed.");
}
