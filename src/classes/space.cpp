#include "space.h"


Box Space::boundingbox()
{
    Box result(0,0,0,0,0,0);
    if (pallocd)
    {
        int i;
        for (i=0; i<pallocd; i++)
        {
            if (spartials[i].s.radius < min_partial_radius) break;
            Point size = Point(spartials[i].s.radius, spartials[i].s.radius, spartials[i].s.radius);
            Point corner1 = spartials[i].s.center.subtract(size);
            Point corner2 = spartials[i].s.center.subtract(size);
            if (!i) result = Box(corner1, corner2);
            else result = result.outer(Box(corner1, corner2));
        }
    }

    return result; 
}

float Space::partial_intersects_cavity(SPartial p)
{
    if (!pallocd) return 0;
    int i;
    float result = 0;
    for (i=0; i<pallocd; i++)
    {
        if (spartials[i].s.radius < min_partial_radius) break;
        // float inter = sphere_intersection(spartials[i].s.radius, p.s.radius, p.s.center.get_3d_distance(spartials[i].s.center));
        float inter = spartials[i].s.radius + p.s.radius - p.s.center.get_3d_distance(spartials[i].s.center);
        // inter /= p.s.volume();
        if (inter > result) result = inter;
        // result += inter;
    }

    return result;
}

void Space::add_partial(SPartial p)
{
    if (!pallocd)
    {
        pallocd = 256;
        spartials = new SPartial[pallocd+4];
        priority = false;
    }

    if (!p.s.center.x && !p.s.center.y && !p.s.center.z) return;

    int i, j;
    for (i=0; i<pallocd; i++) if (spartials[i].s.radius < min_partial_radius) break;             // Get count.

    if (i >= pallocd-4)
    {
        SPartial* old = spartials;
        spartials = new SPartial[pallocd+260];
        for (j=0; j<i; j++) spartials[j] = old[j];
        delete[] old;
        pallocd += 256;
    }

    spartials[i] = p;
    spartials[i+1].s.radius = 0;
    if (p.priority) this->priority = true;
    cached_volume = -1;
    // cout << "Cavity " << this << " added partial at " << spartials[i].s.center << " radius " << spartials[i].s.radius << endl << flush;
}

void Space::output_ngl_js(FILE* fp)
{
    if (!fp) return;
    if (!pallocd) return;

    int i;
    if (priority) fprintf(fp, "// PRIORITY:\n");
    fprintf(fp, "var shape = new NGL.Shape( \"shape\" );\n");
    fprintf(fp, "var sphereBuffer = new NGL.SphereBuffer(\n{\n");
    fprintf(fp, "    position: new Float32Array( [ ");
    for (i=0; i<pallocd; i++)
    {
        if (spartials[i].s.radius < min_partial_radius) break;
        if (i) fprintf(fp, ", ");
        fprintf(fp, "%f, %f, %f", spartials[i].s.center.x, spartials[i].s.center.y, spartials[i].s.center.z);
        if (!(i%5)) fprintf(fp, "\n\t\t");
    }
    fprintf(fp, " ] ),\n");
    fprintf(fp, "    color: new Float32Array( [ ");
    for (i=0; i<pallocd; i++)
    {
        if (spartials[i].s.radius < min_partial_radius) break;
        if (i) fprintf(fp, ", ");
        if (spartials[i].priority) fprintf(fp, "0.0, 1.0, 0.0");
        else if (spartials[i].metallic) fprintf(fp, "0.7, 0.5, 0.3");
        else if (spartials[i].chargedp && spartials[i].chargedn) fprintf(fp, "1, 0, 1");
        else if (spartials[i].chargedp) fprintf(fp, "0.1, 0.1, 1");
        else if (spartials[i].chargedn) fprintf(fp, "1, 0.1, 0.1");
        else if (spartials[i].polar) fprintf(fp, "0.1, 1, 1");
        else if (spartials[i].thio) fprintf(fp, "1, 0.8, 0.1");
        else if (spartials[i].pi) fprintf(fp, "0.8, 0.6, 0.8");
        else fprintf(fp, "0.6, 0.6, 0.6");
        if (!(i%5)) fprintf(fp, "\n\t\t");
    }
    fprintf(fp, " ] ),\n");
    fprintf(fp, "    radius: new Float32Array( [ ");
    for (i=0; i<pallocd; i++)
    {
        if (spartials[i].s.radius < min_partial_radius) break;
        if (i) fprintf(fp, ", ");
        fprintf(fp, "%f", spartials[i].s.radius);
        if (!(i%5)) fprintf(fp, "\n\t\t");
    }
    fprintf(fp, " ] ),\n");
    fprintf(fp, "} );\n");
    fprintf(fp, "shape.addBuffer( sphereBuffer );\n");
    fprintf(fp, "var shapeComp = stage.addComponentFromObject( shape );\n");
    fprintf(fp, "shapeComp.addRepresentation( \"buffer\", { opacity: 0.3 } );\n");
    // fprintf(fp, "shapeComp.autoView();\n");
    fprintf(fp, "\n");
}

int Space::count_partials()
{
    if (!pallocd) return 0;
    int i;
    for (i=0; i<pallocd; i++) if (spartials[i].s.radius < min_partial_radius) return i;
    return 0;
}

Point Space::get_center()
{
    if (!pallocd) return Point(0,0,0);
    int i, j=0;
    Point foravg[pallocd+4];
    for (i=0; i<pallocd; i++)
    {
        if (spartials[i].s.radius < min_partial_radius) break;
        foravg[j] = spartials[i].s.center;
        foravg[j].weight = spartials[i].s.radius;
        j++;
    }

    return average_of_points(foravg, j);
}

float Space::get_volume()
{
    if (!pallocd) return 0;
    if (cached_volume> 0) return cached_volume;

    Box bounds = boundingbox();
    float resolution = 0.5;
    float x, y, z;
    int present = 0;
    for (x=bounds.x1; x<=bounds.x2; x+= resolution)
    {
        for (y=bounds.y1; y<=bounds.y2; y+= resolution)
        {
            for (z=bounds.z1; z<=bounds.z2; z+= resolution)
            {
                Point pt(x,y,z);
                if (this->point_inside_pocket(pt)) present++;
            }
        }
    }
    cached_volume = (double)present * (resolution*resolution*resolution);
    return cached_volume;

    #if 0
    int i, j;
    float result = 0;
    for (i=0; i<pallocd; i++)
    {
        if (spartials[i].s.radius < min_partial_radius) break;
        result += spartials[i].s.volume();
    }

    for (i=0; i<pallocd; i++)
    {
        if (spartials[i].s.radius < min_partial_radius) break;
        for (j=i+1; j<pallocd; j++)
        {
            if (spartials[j].s.radius < min_partial_radius) break;
            result -= sphere_intersection(spartials[i].s.radius, spartials[j].s.radius, spartials[i].s.center.get_3d_distance(spartials[j].s.center));
        }
    }

    return result;
    #endif
}

Point Space::nearest_surface_vertex(Point pt)
{
    if (!vdw_vertex_count) return Point(0,0,0);
    int i;
    float bestr;
    Point result;
    for (i=0; i<vdw_vertex_count; i++)
    {
        float r = pt.get_3d_distance(vdw_surface[i]);
        if (!i || r < bestr)
        {
            bestr = r;
            result = vdw_surface[i];
        }
    }

    return result;
}

SPartial* Space::get_nearest_partial(Point pt)
{
    if (!pallocd || !spartials) return nullptr;

    int i;
    float bestr;
    SPartial* result;
    for (i=0; i<pallocd && spartials[i].s.radius; i++)
    {
        float r = spartials[i].s.center.get_3d_distance(pt);
        if (!i || r < bestr)
        {
            bestr = r;
            result = &spartials[i];
        }
    }

    return result;
}

void Space::compute_vdW_surface(float d)
{
    if (!pallocd || !spartials) return;

    int maxpoints = pallocd * d * d / 3 + 256;
    if (!vdw_surface)
    {
        vdw_surface = new Point[maxpoints];
        vdw_vertex_partial = new SPartial*[maxpoints];
    }

    float halfstep = M_PI / d;
    float step = halfstep * 2;

    int i, ivdW = 0;
    Vector v;
    for (i=0; i<pallocd && spartials[i].s.radius; i++)
    {
        Point ploc = spartials[i].s.center;
        v.r = spartials[i].s.radius;
        float ystep = step / v.r / v.r;
        for (v.theta = -square; v.theta <= square; v.theta += step)
        {
            float xstep = step / v.r / fmax(cos(v.theta), 0.000001);
            float end = M_PI*2-xstep/2;
            for (v.phi = 0; v.phi < end; v.phi += xstep)
            {
                Point pt = ploc.add(v);
                SPartial* np = this->get_nearest_partial(pt);
                if (np != &spartials[i] && pt.get_3d_distance(np->s.center) < np->s.radius) continue;
                if (!pt.x && !pt.y && !pt.z) pt = Point(-0.001, 0.001, -0.001);
                vdw_vertex_partial[ivdW] = &spartials[i];
                vdw_surface[ivdW++] = pt;
                if (ivdW >= maxpoints)
                {
                    cout << "Too many cavity surface ligand_vertices. Please increase limit in code." << endl;
                    throw 0xbadc0de;
                }
            }
        }
    }
    vdw_vertex_count = ivdW;
}

SPartial* Space::point_inside_pocket(Point pt)
{
    if (!pallocd) return nullptr;
    int i, j=-1;
    float rbest = Avogadro;
    for (i=0; i<pallocd; i++)
    {
        if (spartials[i].s.radius < min_partial_radius) break;
        float r = spartials[i].s.center.get_3d_distance(pt);
        if (r < spartials[i].s.radius && r < rbest)
        {
            rbest = r;
            j = i;
        }
    }
    if (j >= 0) return &spartials[j];

    return nullptr;
}

float Space::sphere_inside_pocket(Sphere s, SPartial** p)
{
    if (p) *p = nullptr;
    if (!pallocd) return 0;

    int i;
    float result = 0, br = Avogadro, br1 = Avogadro;
    SPartial *p0 = nullptr, *p1 = nullptr;
    for (i=0; i<pallocd; i++)
    {
        if (spartials[i].s.radius < min_partial_radius) break;
        float r = spartials[i].s.center.get_3d_distance(s.center);
        if (!p0 || r < br)
        {
            p1 = p0;
            br1 = br;
            p0 = &spartials[i];
            if (p) *p = p0;
            br = r;
            // cout << r << endl;
        }
        else if (!p1 || r < br1)
        {
            p1 = &spartials[i];
            br1 = br;
        }
    }

    if (p0 && p1)
    {
        float r = s.center.get_distance_to_line(p0->s.center, p1->s.center);
        float rlim = 1.414*(fmax(p0->s.radius, p1->s.radius) - s.radius + global_clash_allowance);
        if (r <= rlim) result = 1;
        else if ((r - rlim) < s.radius) result = (r - rlim) / s.radius;
        else result = 0;
    }
    else if (p0)
    {
        float r = s.center.get_3d_distance(p0->s.center);
        float rlim = p0->s.radius - s.radius + global_clash_allowance;
        if (r <= rlim) result = 1;
        else if ((r - rlim) < s.radius) result = (r - rlim) / s.radius;
        else result = 0;
    }

    return result;
}

float Space::atom_inside_pocket(Atom *a, bool match_attributes)
{
    float best = 0;
    int i;
    Sphere s = a->get_sphere();
    for (i=0; i<pallocd; i++)
    {
        float f = s.intersection(spartials[i].s);
        if (f >= 0.99999) return f;
        if (f > best) best = f;
    }
    return best;
}

float SPartial::atom_match_score(Atom* a)
{
    float result = 0;
    if (!a) return result;
    if (chargedp || chargedn)
    {
        float achg = a->get_charge();
        if (achg > 0 && chargedn) result += 0.3;
        else if (achg < 0 && chargedp) result += 0.3;
    }

    if (metallic)
    {
        int fam = a->get_family();
        int Z = a->Z;

        if (fam == PNICTOGEN) result += 1.5 / sqrt(Z/8);
        else if (fam == CHALCOGEN)
        {
            if (Z < 10) result += 0.5;
            else result += 2.5 / sqrt(Z/16);
        }
    }

    if (polar)
    {
        float apol = a->is_polar();
        if (fabs(apol) >= hydrophilicity_cutoff) result += 0.2 * fabs(apol);
    }

    if (thio && a->is_thio()) result += 0.05;

    if (this->pi && a->is_pi()) result += 0.12;

    if (priority) result *= 5;

    return result;
}
