#include "../classes/configtext.h"
#include <cassert>
#include <cstdio>

int main()
{
    char line[128];
    for (auto expected : {std::array<float,3>{{10.253f,10.131f,8.020f}},
                          std::array<float,3>{{-123.456f,1234.567f,-98.765f}}})
    {
        snprintf(line,sizeof(line),"%-30s%8.3f%8.3f%8.3f", "ATOM      1   C1 LIG A   1",
            expected[0],expected[1],expected[2]);
        auto actual = parse_pdb_coordinates(line);
        for (int axis=0; axis<3; ++axis) assert(fabs(actual[axis]-expected[axis]) < 0.0001f);
    }
    for (const char* invalid : {"ATOM", "ATOM      1   C1 LIG A   1         bad  10.131   8.020"})
    {
        bool rejected = false;
        try { parse_pdb_coordinates(invalid); }
        catch (const std::invalid_argument&) { rejected = true; }
        assert(rejected);
    }
    puts("External pose coordinates preserve all three PDB columns.");
}
