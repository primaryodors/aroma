
#include <iostream>
#include <math.h>
#include <time.h>
#include "../classes/point.h"

using namespace std;

int main (int argc, char** argv)
{
    Point pt(2,-5,3);
    cout << "Point pt is located at " << pt.x << "," << pt.y << "," << pt.z << ".\n";

    Point pt1(1,1,1);
    Point pt2 = pt1.negate();
    cout << "Points pt1 and pt2 are " << pt1.get_3d_distance(&pt2) << " Angstroms apart.\n";

    Point pt0;
    Point pt3 = pt.multiply_3d_distance(&pt0, 0.1);

    if (pt.pt_in_bounding_box(&pt1, &pt2))
        cout << "Point pt is in bounding box [pt1, pt2].\n";
    else
        cout << "Point pt is not in bounding box [pt1, pt2].\n";

    if (pt3.pt_in_bounding_box(&pt1, &pt2))
        cout << "Point pt3 is in bounding box [pt1, pt2].\n";
    else
        cout << "Point pt3 is not in bounding box [pt1, pt2].\n";

    cout << "Points pt and pt1 are "
         << find_3d_angle(&pt, &pt1, 0) * 180.0 / M_PI
         << " degrees apart relative to [0,0,0].\n";

    Point pt5(-1.6,3.1,-0.7);
    Vector v(&pt5);
    Point pt6(&v);

    cout << "Vector v from point pt5 has r=" << v.r
         << ", theta=" << v.theta * 180.0 / M_PI
         << ", phi=" << v.phi * 180.0 / M_PI
         << ".\n";

    cout << "Point pt6 from Vector v is " << (0.01 * round(pt6.get_3d_distance(&pt5)*100)) << " Angstroms from pt5.\n";

    Point ptrot8, ptref(pt5.x + 0, pt5.y + 0, pt5.z + 1);

    float f, f57, f_57, test;
    for (f=0; f<M_PI*2; f+=0.1)
    {
        ptrot8.x = pt5.x + sin(f);
        ptrot8.y = pt5.y + 0;
        ptrot8.z = pt5.z + cos(f);

        f57 = f*fiftyseven;
        f_57 = 360 - f57;
        test = find_3d_angle(ptrot8, ptref, pt5)*fiftyseven;

        cout << ((fabs(test-f57) > 0.01 && fabs(test-f_57) > 0.01) ? "" : "# ")
             << "Rotation " << f57 << " points are "
             << test << " degrees apart."
             << endl;
    }

    Vector axis(1.0, fiftyseventh*-53, fiftyseventh*281);
    for (f=0; f<M_PI*2; f+=0.1)
    {
        ptrot8 = rotate3D(ptref, pt5, axis, f);
        float th3d = find_3d_angle(ptrot8, ptref, pt5);
        float thvec = find_angle_along_vector(ptrot8, ptref, pt5, axis);
        Point ptrot9 = rotate3D(ptrot8, pt5, axis, thvec);

        f57 = f*fiftyseven;
        f_57 = 360 - f57;

        float th3d57 = th3d*fiftyseven;
        float thvec57 = thvec*fiftyseven;
        test = ptrot9.get_3d_distance(ptref);

        cout << ((test >= 0.0001) ? "" : "# ")
            << "Rotation " << f57
            << " find_3d_angle() gives " << th3d57 << "deg"
            << " find_angle_along_vector() gives " << thvec57 << "deg"
            << " putitback anomaly = " << test
            << endl;
    }
}











