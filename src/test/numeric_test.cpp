#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstdio>

int main()
{
    volatile double invalid = std::strtod("nan", nullptr);
    volatile double unbounded = std::strtod("inf", nullptr);
    assert(std::isnan(invalid));
    assert(!std::isfinite(invalid));
    assert(std::isinf(unbounded));
    assert(!std::isfinite(unbounded));
    puts("Floating-point guards passed.");
}
