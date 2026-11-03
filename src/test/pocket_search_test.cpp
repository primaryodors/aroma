#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdio>
#include "../classes/protein.h"

int main()
{
    Protein p;
    FILE* f = fopen("pdbs/OR51/OR51E2.active.pdb", "r");
    assert(f);
    p.load_pdb(f);
    fclose(f);
    Point center = p.get_residue(150)->get_CA_location();
    Point expected = center;
    float best = 0;
    auto began = std::chrono::steady_clock::now();
    for (float x=-1; x<=1; x+=0.25)
    for (float y=-1; y<=1; y+=0.25)
    for (float z=-1; z<=1; z+=0.25)
    {
        if (sqrt(x*x+y*y+z*z) > 1) continue;
        Point sample(static_cast<float>(center.x+x), static_cast<float>(center.y+y),
                     static_cast<float>(center.z+z));
        float nearest = Avogadro;
        for (int r=p.get_start_resno(); r<=p.get_end_resno(); r++)
        {
            AminoAcid* aa = p.get_residue(r);
            if (!aa) continue;
            Atom* a = aa->get_nearest_atom(sample);
            if (a) nearest = fmin(nearest, a->loc.get_3d_distance(sample));
        }
        if (nearest > best && nearest < 1e9) { best = nearest; expected = sample; }
    }
    auto reference_end = std::chrono::steady_clock::now();
    float clearance = 0;
    Point actual = p.find_loneliest_point(center, Point(2,2,2), &clearance);
    auto grid_end = std::chrono::steady_clock::now();
    assert(actual.get_3d_distance(expected) < 0.001);
    assert(fabs(clearance-best) < 0.001);
    double reference_ms = std::chrono::duration<double,std::milli>(reference_end-began).count();
    double grid_ms = std::chrono::duration<double,std::milli>(grid_end-reference_end).count();
    printf("Pocket search passed: reference %.2f ms; grid %.2f ms.\n", reference_ms, grid_ms);
}
