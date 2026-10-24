#include <cassert>
#include <cstdio>
#include "../classes/space.h"

int main()
{
    Space space;
    SPartial sphere;
    sphere.s.center = Point(1,2,3);
    sphere.s.radius = 3;
    space.add_partial(sphere);
    Box bounds = space.boundingbox();
    assert(bounds.x1 == -2 && bounds.x2 == 4);
    assert(bounds.y1 == -1 && bounds.y2 == 5);
    assert(bounds.z1 == 0 && bounds.z2 == 6);
    float volume = space.get_volume();
    assert(volume > 100 && volume < 125);

    Space shifted;
    sphere.s.center = Point(11,22,33);
    shifted.add_partial(sphere);
    assert(shifted.get_volume() == volume);
    sphere.s.center = Point(-8,-6,2);
    sphere.s.radius = 2.5;
    space.add_partial(sphere);
    bounds = space.boundingbox();
    assert(bounds.x1 == -10.5 && bounds.x2 == 4);
    assert(bounds.y1 == -8.5 && bounds.y2 == 5);
    assert(bounds.z1 == -0.5 && bounds.z2 == 6);
    assert(space.get_volume() > volume);
    puts("Space bounds include both extents and preserve translated volumes.");
}
