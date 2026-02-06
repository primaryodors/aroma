
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <math.h>

#include "misc.h"

#ifndef _POINT
#define _POINT

using namespace std;

struct Vector;

struct Point
{
    double x=0, y=0, z=0;
    float weight = 1;

    Point();
    Point(Vector v);
    Point(Vector* v);
    Point(double x, double y, double z);
    Point add(Point add_to) const
    {
        return add(&add_to);
    }
    Point add(Point* add_to) const;
    Point add(Vector* add_to) const;
    Point subtract(const Point subtracted) const;
    Point subtract(const Point* subtracted) const;
    Point subtract(Vector* subtracted) const
    {
        Point pt(subtracted);
        return subtract(pt);
    }
    Point negate();
    Point randomize(float amt) const;
    float get_3d_distance(const Point reference) const;
    float get_3d_distance(const Point* reference) const;
    float get_distance_to_line(const Point a, const Point b) const;         // Where a and b are the termini of the line.
    Point multiply_3d_distance(const Point* reference, float r_mult);
    bool pt_in_bounding_box(const Point* corner1, const Point* corner2);
    float magnitude() const;
    void scale(float new_magn);
    void multiply(float multiplier);
    bool fits_inside(Point container);

    std::string printable() const;

    Point& operator=(Vector v);
    friend std::ostream& operator<<(std::ostream& os, const Vector& v);
};

struct Box
{
    double x1=0, y1=0, z1=0;
    double x2=0, y2=0, z2=0;

    Box();
    Box(Point pt1, Point pt2);
    Box(double X1, double X2, double Y1, double Y2, double Z1, double Z2);
    Point size();
    Point center();
    float width() const { return fabs(x2-x1); }
    float height() const { return fabs(y2-y1); }
    float depth() const { return fabs(z2-z1); }
    Box outer(Box other);
    void accommodate(Point p);

    double volume();
    bool contains(Point pt);
};

struct Sphere
{
    Point center;
    float radius = 0;
    float area() { return 4.0 * M_PI*radius*radius; }
    float volume() { return 4.0/3 * M_PI*radius*radius*radius; }
};

struct Vector
{
    double r=0, theta=0, phi=0;

    Vector()
    {
        r = theta = phi = 0;
    }
    Vector(const Point from);
    Vector(const Point* from);
    Vector(double r, double theta, double phi);
    Vector negate()
    {
        Point pt(this);
        pt = pt.negate();
        Vector v(&pt);
        return v;
    }
    Vector add(Vector v)
    {
        return add(&v);
    };
    Vector add(Vector* v);

    std::string printable() const;
    int octant_idx();           // Assign octants of relative 3D space to indices for an array. The opposite of octant i always equals 7-i.

    Vector& operator=(Point p);
    friend std::ostream& operator<<(std::ostream& os, const Point& p);
};

struct LocatedVector : public Vector
{
    LocatedVector()
    {
        r = theta = phi = origin.x = origin.y = origin.z = 0;
    }
    LocatedVector(Vector sc)
    {
        copy(sc);
    }
    LocatedVector add(Vector v)
    {
        LocatedVector result = Vector::add(&v);
        result.origin = origin;
        return result;
    }
    void copy(Vector sc)
    {
        r=sc.r;
        theta=sc.theta;
        phi=sc.phi;
    }
    Point to_point();
    Point origin;
};

struct Rotation
{
    Vector v;
    float a=0;
    Rotation()
    {
        a = 0;
    };

    Rotation add(Rotation rot)
    {
        return add(&rot);
    };
    Rotation add(Rotation* rot);

    char* printable();
};

struct LocRotation : public Rotation
{
    Point origin;

    LocRotation();
    LocRotation(Rotation rot);
    LocRotation(LocatedVector lv);
    LocatedVector get_lv();
};

class Atom;

struct Tug
{
    Vector vec;
    Rotation rot;

    Tug() { ; }
    Tug(Atom* atom, Point molcen, Vector pull);

    Tug add(Tug t);
};

Point average_of_points(Point* points, int count);
Point size_of_point_space(Point* points, int count);
float equidistance_anomaly(Point point, Point* refs, int count);
Point find_equidistant_point(Point* points, int count, Point* bias = nullptr);

float find_angle(float dx, float dy);
float find_angle_delta(float a1, float a2);
float find_3d_angle(Point* A, Point* B, Point* source);
float find_3d_angle(Point A, Point B, Point source);
float find_angle_along_vector(Point* pt1, Point* pt2, Point* source, Vector* v);
float find_angle_along_vector(Point pt1, Point pt2, Point source, Vector v);

Point rotate3D(Point* point, Point* source, Vector* axis, float theta);
Point rotate3D(Point point, Point source, Vector axis, float theta);
Point rotate3D(Point* point, Point* source, Rotation* rot);
Point rotate3D(Point point, Point source, Rotation rot);

Rotation align_points_3d(Point* point, Point* align, Point* center);
Rotation align_points_3d(Point point, Point align, Point center);
Rotation* align_2points_3d(Point* point1, Point* align1, Point* point2, Point* align2, Point* center);

Vector compute_normal(Point* pt1, Point* pt2, Point* pt3);
Vector compute_normal(Point pt1, Point pt2, Point pt3);

float are_points_planar(Point p1, Point p2, Point p3, Point p4);
float polygon_radius(float side_length, int num_sides);

float cap_volume(float r1, float r2, float d);                  // TODO:
float cap_area(float r1, float r2, float d);                    // Important: this is the OUTER (i.e. NON-intersecting) cap of the sphere.
float sphere_intersection(float r1, float r2, float d);         // Volume of the lens composed of the two caps.
float sphere_inter_area(float r1, float r2, float d);           // Area of the circle formed by the spheres' intersection.

Vector v_from_pt_sub(Point distal, Point reference);

std::ostream& operator<<(std::ostream& os, const Point& p);
std::ostream& operator<<(std::ostream& os, const Vector& v);
std::ostream& operator<<(std::ostream& os, const Rotation& r);

#endif

