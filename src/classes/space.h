#include "atom.h"
#include "progress.h"

#ifndef _SPACE
#define _SPACE

#define min_partial_radius 0.7
#define min_dist_bounding_box 11
#define cav_360_step fiftyseventh*8
#define cav_xyz_step 1.1
#define cav_min_partials 4
#define cav_linking_threshold 2.8

class SPartial
{
    public:
    Sphere s;
    bool chargedp = false;
    bool chargedn = false;
    bool metallic = false;
    bool polar = false;
    bool thio = false;
    bool pi = false;
    bool priority = false;

    float atom_match_score(Atom* a);
};

class Space
{
    public:
    SPartial* get_nearest_partial(Point pt);
    float partial_intersects_cavity(SPartial p);
    void add_partial(SPartial p);
    void output_ngl_js(FILE* fp);
    int count_partials();
    SPartial* get_partial_by_idx(int idx) { return &spartials[idx]; }
    Point get_center();
    float get_volume();
    SPartial* point_inside_pocket(Point pt);
    float sphere_inside_pocket(Sphere s, SPartial** partial = nullptr);
    float atom_inside_pocket(Atom* a, bool match_attributes = false);
    float space_intersection(Space* other);
    Box boundingbox();

    protected:
    void compute_vdW_surface(float d);
    Point nearest_surface_vertex(Point pt);

    SPartial* spartials = nullptr;
    int pallocd = 0;
    bool priority = false;
    Point* vdw_surface = nullptr;
    SPartial** vdw_vertex_partial = nullptr;
    int vdw_vertex_count = 0;
    float cached_volume = -1;
};

#endif