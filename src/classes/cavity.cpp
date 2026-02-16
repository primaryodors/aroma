
#include "cavity.h"

float cav_xmax = Avogadro, cav_xmin = -Avogadro, cav_ymax = Avogadro, cav_ymin = -Avogadro, cav_zmax = Avogadro, cav_zmin = -Avogadro;
float cav_xyrlim = Avogadro, cav_xzrlim = Avogadro, cav_yzrlim = Avogadro;
int cav_resmin = -99999, cav_resmax = 99999;

float Cavity::cavity_intersection(Cavity* other)
{
    float result = 0;
    if (!partials || !pallocd || !other->partials || !other->pallocd) return result;
    int i, j;

    for (i=0; i<pallocd && partials[i].s.radius >= min_partial_radius; i++)
    {
        float r1 = partials[i].s.radius;
        for (j=0; j<other->pallocd && other->partials[j].s.radius >= min_partial_radius; j++)
        {
            float r2 = other->partials[j].s.radius;
            float d = partials[i].s.center.get_3d_distance(other->partials[j].s.center);
            if (d > r1+r2) continue;
            float si = sphere_intersection(r1, r2, d);

            /*if (d < 2)
                cout << "overlap: " << partials[i].s.center << other->partials[j].s.center << " " << r1 << " " << r2 << " " << d << " " << si << endl;*/

            result += si;
            // if (d < 2) cout << " " << d;
        }
    }

    return result;
}

void Cavity::unify(Cavity* cavfrom)
{
    int i, j, n;

    n = 0;
    if (partials) for (; n<pallocd && partials[n].s.radius >= min_partial_radius; n++);
    j = n;
    if (cavfrom->partials)
    {
        for (i=0; i<cavfrom->pallocd && cavfrom->partials[i].s.radius >= min_partial_radius; i++);
        n += i;
    }
    if (n > pallocd)
    {
        pallocd += n + 256;
        CPartial* oldpart = partials;
        partials = new CPartial[pallocd+4];
        j=0;
        if (oldpart)
        {
            for (i=0; oldpart[i].s.radius >= min_partial_radius; i++) partials[i] = oldpart[i];
            j = i;
            delete oldpart;
        }
    }

    if (cavfrom->partials)
    {
        for (i=0; cavfrom->partials[i].s.radius >= min_partial_radius; i++) partials[i+j] = cavfrom->partials[i];
        delete cavfrom->partials;
    }
}

Box Cavity::boundingbox()
{
    Box result(0,0,0,0,0,0);
    if (pallocd)
    {
        int i;
        for (i=0; i<pallocd; i++)
        {
            if (partials[i].s.radius < min_partial_radius) break;
            Point size = Point(partials[i].s.radius, partials[i].s.radius, partials[i].s.radius);
            Point corner1 = partials[i].s.center.subtract(size);
            Point corner2 = partials[i].s.center.subtract(size);
            if (!i) result = Box(corner1, corner2);
            else result = result.outer(Box(corner1, corner2));
        }
    }

    return result; 
}

int Cavity::scan_in_protein(Protein* p, Cavity* cavs, int cmax, Progressbar* pgb)
{
    if (!p || !cavs) return 0;
    if (cmax < 1) return 0;

    int i, j, l, n, sr = max(p->get_start_resno(), cav_resmin), er = min(p->get_end_resno(), cav_resmax);
    float x, y, z, step, yoff = 0, zoff = 0;
    Point pcen = p->get_region_center(sr, er), pbox = p->get_region_bounds(sr, er);

    int priorities[1024], pqty;
    priorities[0] = pqty = 0;
    l = p->get_end_resno();
    for (i=1; i<=l; i++)
    {
        AminoAcid* a = p->get_residue(i);
        if (a && a->priority) priorities[pqty++] = i;
    }
    priorities[pqty] = 0;

    step = cav_xyz_step;
    AminoAcid* can_clash[SPHREACH_MAX+4];
    Molecule dummy("DUMMY");
    dummy.add_atom("H", "H1", nullptr, 0);
    Point size(_INTERA_R_CUTOFF, _INTERA_R_CUTOFF, _INTERA_R_CUTOFF);
    CPartial parts[65536];
    j=0;
    bool any_priority = false;

    float xmin = pcen.x - pbox.x + min_dist_bounding_box, xmax = pcen.x + pbox.x - min_dist_bounding_box;
    if (pgb)
    {
        pgb->minimum = xmin;
        pgb->maximum = xmax;
        cout << endl;
    }
    else cout << "Cavity search in progress..." << flush;

    for (x = xmin; x <= xmax; x += step)
    {
        if (pgb) pgb->update(x);
        else cout << "." << flush;
        yoff = yoff ? 0 : step/2;
        for (y = pcen.y - pbox.y + min_dist_bounding_box - yoff; y <= pcen.y + pbox.y - min_dist_bounding_box; y += step)
        {
            zoff = zoff ? 0 : step/2;
            for (z = pcen.z - pbox.z + min_dist_bounding_box - zoff; z <= pcen.z + pbox.z - min_dist_bounding_box; z += step)
            {
                Point pt(x,y,z);
                dummy.recenter(pt);
                int sphres = p->get_residues_can_clash_ligand(can_clash, &dummy, pt, size, priorities, true);
                if (sphres < 8+pqty) continue;          // Too isolated.

                float occltot = 0;
                for (i=0; i<sphres; i++)
                {
                    float f = can_clash[i]->octant_occlusion();
                    occltot += f;
                }
                float occlavg = occltot / i;
                // cout << occlavg << endl;
                if (occlavg < 0.7) continue;            // Too isolated.

                float rmin;
                CPartial working;
                for (i=0; i<sphres; i++)
                {
                    Atom* a = can_clash[i]->get_nearest_atom(pt);
                    float r = a->loc.get_3d_distance(pt) - a->vdW_radius + global_clash_allowance*1.5;

                    if (!a->is_backbone)
                    {
                        // Add a slight "give" for side chains of aminos with greater flexional probability.
                        float aafp = can_clash[i]->get_aa_definition()->flexion_probability;
                        if (aafp >= 0.075) r += aafp * 2.9 * max(0, a->get_Greek()-2);
                    }

                    if (r < rmin || !i)
                    {
                        rmin = r;
                    }
                    if (r < min_partial_radius) break;
                }

                for (i=0; i<sphres; i++)
                {
                    Atom* a = can_clash[i]->get_nearest_atom(pt);
                    float r = a->loc.get_3d_distance(pt)-a->vdW_radius;
                    if (r > rmin /* && !can_clash[i]->priority */) continue;

                    Atom* CB = can_clash[i]->get_atom("CB");
                    if (CB
                        && a->loc.get_3d_distance(CB->loc) < a->loc.get_3d_distance(can_clash[i]->get_CA_location())
                        && (!pqty || can_clash[i]->priority)
                        )
                    {
                        if (can_clash[i]->priority)
                        {
                            // cout << pt << can_clash[i]->get_name() << " distance " << r << " has priority." << endl;
                            any_priority = /*working.priority =*/ true;
                        }
                        if (can_clash[i]->coordmtl) working.metallic = true;
                        if (can_clash[i]->get_charge() < -hydrophilicity_cutoff) working.chargedn = true;
                        if (can_clash[i]->get_charge() >  hydrophilicity_cutoff) working.chargedp = true;
                        if (can_clash[i]->pi_stackability() >= 0.2) working.pi = true;
                        if (fabs(can_clash[i]->hydrophilicity()) > hydrophilicity_cutoff) working.polar = true;
                        if (can_clash[i]->count_atoms_by_element("S")) working.thio = true;
                    }
                }

                if (rmin >= min_partial_radius)
                {
                    working.s.center = pt;
                    working.s.radius = rmin;
                    // cout << "Found partial at " << pt << " radius " << rmin << endl << flush;
                    parts[j++] = working;
                }
            }
        }
    }
    if (pgb) pgb->erase();
    else cout << endl;

    // Now consolidate all partials into glommed cavities.
    Cavity tmpcav[4096];
    int pmax = j;
    j=0;
    for (i=0; i<pmax; i++)
    {
        if (parts[i].s.center.magnitude() == 0) break;
        // if (parts[i].resno && (parts[i].resno < sr || parts[i].resno > er)) continue;
        bool glommed = false;
        for (l=0; l<j; l++)
        {
            float inter = tmpcav[l].partial_intersects_cavity(parts[i]);
            if (/*(parts[i].priority && tmpcav[l].priority) ||*/ inter >= cav_linking_threshold)
            {
                tmpcav[l].add_partial(parts[i]);
                glommed = true;
                /*cout << "Partial at " << parts[i].s.center << " radius " << parts[i].s.radius << " belongs to cavity #" << l
                    << " with intersection " << inter << endl << flush;*/
                break;
            }
        }
        if (!glommed)
        {
            tmpcav[j++].add_partial(parts[i]);
        }
        if (j >= 4090) break;
    }
    l=j;

    // Any cavities that intersect more than a threshold amount, unify them.
    for (i=0; i<l; i++)
    {
        Point cen = tmpcav[i].get_center();
        if (cen.x < cav_xmin || cen.x > cav_xmax) continue;
        if (cen.y < cav_ymin || cen.y > cav_ymax) continue;
        if (cen.z < cav_zmin || cen.z > cav_zmax) continue;

        for (j=i+1; j<l; j++)
        {
            cen = tmpcav[j].get_center();
            if (cen.x < cav_xmin || cen.x > cav_xmax) continue;
            if (cen.y < cav_ymin || cen.y > cav_ymax) continue;
            if (cen.z < cav_zmin || cen.z > cav_zmax) continue;
            float u = tmpcav[i].cavity_intersection(&tmpcav[j]);
            if (u >= cavity_intersect_threshold)
            {
                // cout << "Cavities " << i << " and " << j << " intersect by " << u << endl;
                tmpcav[i].unify(&tmpcav[j]);
                for (n=j+1; n<l; n++) tmpcav[n-1] = tmpcav[n];
                l--;
                j = i;
            }
        }
    }

    // cout << "l = " << l << endl;

    // if (any_priority) cout << "Priority residues found." << endl;
    j=0;
    for (i=0; i<l; i++)
    {
        tmpcav[i].prot = p;
        Point cen = tmpcav[i].get_center();
        if (cen.x < cav_xmin || cen.x > cav_xmax) continue;
        if (cen.y < cav_ymin || cen.y > cav_ymax) continue;
        if (cen.z < cav_zmin || cen.z > cav_zmax) continue;

        float r = sqrt(pow(cen.x - pcen.x, 2) + pow(cen.y - pcen.y, 2));
        if (r > cav_xyrlim) continue;
        r = sqrt(pow(cen.x - pcen.x, 2) + pow(cen.z - pcen.z, 2));
        if (r > cav_xzrlim) continue;
        r = sqrt(pow(cen.y - pcen.y, 2) + pow(cen.z - pcen.z, 2));
        if (r > cav_yzrlim) continue;

        // cout << "i = " << i << " (" << tmpcav[i].count_partials() << ")" << endl;

        int x;
        for (x=0; x<pqty; x++)
        {
            AminoAcid* aa = p->get_residue(priorities[x]);
            if (!aa) continue;
            Atom* a = aa->get_nearest_atom(tmpcav[i].get_center());
            if (!a) continue;
            CPartial* part = tmpcav[i].get_nearest_partial(a->loc);
            if (!part) continue;
            r = part->s.center.get_3d_distance(a->loc);
            // if (r > _INTERA_R_CUTOFF+part->s.radius+a->vdW_radius) continue;
            tmpcav[i].priority = part->priority = true;
            if (aa->coordmtl) part->metallic = true;
            if (aa->get_charge() < -hydrophilicity_cutoff) part->chargedn = true;
            if (aa->get_charge() >  hydrophilicity_cutoff) part->chargedp = true;
            if (aa->pi_stackability() >= 0.2) part->pi = true;
            if (fabs(aa->hydrophilicity()) > hydrophilicity_cutoff) part->polar = true;
            if (aa->count_atoms_by_element("S")) part->thio = true;
            break;
        }

        if (tmpcav[i].count_partials() >= cav_min_partials
            && (!any_priority || tmpcav[i].priority)
            )
        {
            cavs[j++] = tmpcav[i];
            // cout << "Accepted " << i << endl;
        }
        if (j >= cmax-1) break;
    }

    // cout << cavs[4].cavity_intersection(&cavs[6]) << endl;

    return j;
}

float Cavity::partial_intersects_cavity(CPartial p)
{
    if (!pallocd) return 0;
    int i;
    float result = 0;
    for (i=0; i<pallocd; i++)
    {
        if (partials[i].s.radius < min_partial_radius) break;
        // float inter = sphere_intersection(partials[i].s.radius, p.s.radius, p.s.center.get_3d_distance(partials[i].s.center));
        float inter = partials[i].s.radius + p.s.radius - p.s.center.get_3d_distance(partials[i].s.center);
        // inter /= p.s.volume();
        if (inter > result) result = inter;
        // result += inter;
    }

    return result;
}

void Cavity::add_partial(CPartial p)
{
    if (!pallocd)
    {
        pallocd = 256;
        partials = new CPartial[pallocd+4];
        priority = false;
    }

    if (!p.s.center.x && !p.s.center.y && !p.s.center.z) return;

    int i, j;
    for (i=0; i<pallocd; i++) if (partials[i].s.radius < min_partial_radius) break;             // Get count.

    if (i >= pallocd-4)
    {
        CPartial* old = partials;
        partials = new CPartial[pallocd+260];
        for (j=0; j<i; j++) partials[j] = old[j];
        delete[] old;
        pallocd += 256;
    }

    partials[i] = p;
    partials[i+1].s.radius = 0;
    if (p.priority) this->priority = true;
    cached_volume = -1;
    // cout << "Cavity " << this << " added partial at " << partials[i].s.center << " radius " << partials[i].s.radius << endl << flush;
}

void Cavity::output_ngl_js(FILE* fp)
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
        if (partials[i].s.radius < min_partial_radius) break;
        if (i) fprintf(fp, ", ");
        fprintf(fp, "%f, %f, %f", partials[i].s.center.x, partials[i].s.center.y, partials[i].s.center.z);
        if (!(i%5)) fprintf(fp, "\n\t\t");
    }
    fprintf(fp, " ] ),\n");
    fprintf(fp, "    color: new Float32Array( [ ");
    for (i=0; i<pallocd; i++)
    {
        if (partials[i].s.radius < min_partial_radius) break;
        if (i) fprintf(fp, ", ");
        if (partials[i].priority) fprintf(fp, "0.0, 1.0, 0.0");
        else if (partials[i].metallic) fprintf(fp, "0.7, 0.5, 0.3");
        else if (partials[i].chargedp && partials[i].chargedn) fprintf(fp, "1, 0, 1");
        else if (partials[i].chargedp) fprintf(fp, "0.1, 0.1, 1");
        else if (partials[i].chargedn) fprintf(fp, "1, 0.1, 0.1");
        else if (partials[i].polar) fprintf(fp, "0.1, 1, 1");
        else if (partials[i].thio) fprintf(fp, "1, 0.8, 0.1");
        else if (partials[i].pi) fprintf(fp, "0.8, 0.6, 0.8");
        else fprintf(fp, "0.6, 0.6, 0.6");
        if (!(i%5)) fprintf(fp, "\n\t\t");
    }
    fprintf(fp, " ] ),\n");
    fprintf(fp, "    radius: new Float32Array( [ ");
    for (i=0; i<pallocd; i++)
    {
        if (partials[i].s.radius < min_partial_radius) break;
        if (i) fprintf(fp, ", ");
        fprintf(fp, "%f", partials[i].s.radius);
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

int Cavity::count_partials()
{
    if (!pallocd) return 0;
    int i;
    for (i=0; i<pallocd; i++) if (partials[i].s.radius < min_partial_radius) return i;
    return 0;
}

Point Cavity::get_center()
{
    if (!pallocd) return Point(0,0,0);
    int i, j=0;
    Point foravg[pallocd+4];
    for (i=0; i<pallocd; i++)
    {
        if (partials[i].s.radius < min_partial_radius) break;
        foravg[j] = partials[i].s.center;
        foravg[j].weight = partials[i].s.radius;
        j++;
    }

    return average_of_points(foravg, j);
}

float Cavity::get_volume()
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
        if (partials[i].s.radius < min_partial_radius) break;
        result += partials[i].s.volume();
    }

    for (i=0; i<pallocd; i++)
    {
        if (partials[i].s.radius < min_partial_radius) break;
        for (j=i+1; j<pallocd; j++)
        {
            if (partials[j].s.radius < min_partial_radius) break;
            result -= sphere_intersection(partials[i].s.radius, partials[j].s.radius, partials[i].s.center.get_3d_distance(partials[j].s.center));
        }
    }

    return result;
    #endif
}

Point Cavity::nearest_surface_vertex(Point pt)
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

CPartial* Cavity::get_nearest_partial(Point pt)
{
    if (!pallocd || !partials) return nullptr;

    int i;
    float bestr;
    CPartial* result;
    for (i=0; i<pallocd && partials[i].s.radius; i++)
    {
        float r = partials[i].s.center.get_3d_distance(pt);
        if (!i || r < bestr)
        {
            bestr = r;
            result = &partials[i];
        }
    }

    return result;
}

void Cavity::compute_vdW_surface(float d)
{
    if (!pallocd || !partials) return;

    int maxpoints = pallocd * d * d / 3 + 256;
    if (!vdw_surface)
    {
        vdw_surface = new Point[maxpoints];
        vdw_vertex_partial = new CPartial*[maxpoints];
    }

    float halfstep = M_PI / d;
    float step = halfstep * 2;

    int i, ivdW = 0;
    Vector v;
    for (i=0; i<pallocd && partials[i].s.radius; i++)
    {
        Point ploc = partials[i].s.center;
        v.r = partials[i].s.radius;
        float ystep = step / v.r / v.r;
        for (v.theta = -square; v.theta <= square; v.theta += step)
        {
            float xstep = step / v.r / fmax(cos(v.theta), 0.000001);
            float end = M_PI*2-xstep/2;
            for (v.phi = 0; v.phi < end; v.phi += xstep)
            {
                Point pt = ploc.add(v);
                CPartial* np = this->get_nearest_partial(pt);
                if (np != &partials[i] && pt.get_3d_distance(np->s.center) < np->s.radius) continue;
                if (!pt.x && !pt.y && !pt.z) pt = Point(-0.001, 0.001, -0.001);
                vdw_vertex_partial[ivdW] = &partials[i];
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

CPartial* Cavity::point_inside_pocket(Point pt)
{
    if (!pallocd) return nullptr;
    int i, j=-1;
    float rbest = Avogadro;
    for (i=0; i<pallocd; i++)
    {
        if (partials[i].s.radius < min_partial_radius) break;
        float r = partials[i].s.center.get_3d_distance(pt);
        if (r < partials[i].s.radius && r < rbest)
        {
            rbest = r;
            j = i;
        }
    }
    if (j >= 0) return &partials[j];

    return nullptr;
}

float Cavity::sphere_inside_pocket(Sphere s, CPartial** p)
{
    if (p) *p = nullptr;
    if (!pallocd) return 0;

    int i;
    float result = 0, br = Avogadro, br1 = Avogadro;
    CPartial *p0 = nullptr, *p1 = nullptr;
    for (i=0; i<pallocd; i++)
    {
        if (partials[i].s.radius < min_partial_radius) break;
        float r = partials[i].s.center.get_3d_distance(s.center);
        if (!p0 || r < br)
        {
            p1 = p0;
            br1 = br;
            p0 = &partials[i];
            if (p) *p = p0;
            br = r;
            // cout << r << endl;
        }
        else if (!p1 || r < br1)
        {
            p1 = &partials[i];
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

float Cavity::molecule_inside_pocket(Molecule* m, bool mattr)
{
    int i, j, n = m->get_atom_count();
    float result = 0;
    for (i=0; i<n; i++)
    {
        Atom* a = m->get_atom(i);
        if (!a) continue;

        Sphere s = a->get_sphere();
        int Z = a->Z;
        if (Z < 2) s.radius /= 2;

        CPartial* p;
        float partial = sphere_inside_pocket(s, &p);

        if (!partial) continue;

        if (mattr)
        {
            bool apol = fabs(a->is_polar()) >= hydrophilicity_cutoff;
            float achg = a->get_charge();
            bool apos = achg >= hydrophilicity_cutoff;
            bool aneg = achg <= -hydrophilicity_cutoff;
            bool api = a->is_pi();
            bool amet = a->is_metal();
            int afam = a->get_family();
            bool acm = (afam == CHALCOGEN || a->get_family() == PNICTOGEN) && !amet;

            // if (p->thio) cout << a->name << " partial was " << partial;

            if (p->metallic && acm) partial *= 5;
            else if (p->thio && amet) partial *= 5;
            else if (p->chargedn && apos) partial *= 1.3;
            else if (p->chargedn && aneg) partial *= 0.6;
            else if (p->chargedp && aneg) partial *= 1.3;
            else if (p->chargedp && apos) partial *= 0.6;
            else if (p->thio && acm) partial *= 1.1;
            // else if (!p->thio && !acm && p->polar != apol) partial *= 0.5;

            if (p->pi && api) partial *= 1.3;
            if (p->polar && apol) partial *= 1.7;

            // if (p->thio) cout << " now " << partial << endl;
        }

        result += partial;
    }

    if (n) result /= n;
    return result;
}

float Cavity::cavity_filling(Molecule *m)
{
    float empties = 0.0f;
    int filleds = 0;
    float x, y, z;
    float step = 0.5;

    Box b = boundingbox();

    for (x=b.x1; x<=b.x2; x+=step)
        for (y=b.y1; y<=b.y2; y+=step)
            for (z=b.z1; z<=b.z2; z+=step)
            {
                Point pt(x,y,z);
                CPartial* part = point_inside_pocket(pt);
                if (part)
                {
                    Atom* a = m->get_nearest_atom(pt);
                    if (a && a->loc.get_3d_distance(pt) <= a->vdW_radius)
                        filleds++;
                    else
                    {
                        if (part->chargedn || part->chargedp) empties += 0.1;
                        else if (part->pi && part->polar) empties += 0.333;
                        else if (part->polar) empties += 0.5;
                        else if (part->thio) empties += 0.8;
                        else if (part->pi) empties += 0.9;
                        else empties += 1;
                    }
                }
            }

    if (!empties && !filleds) return 0;
    return (float)filleds / (empties+filleds);
}

const Point* ligand_vertices;
float Cavity::containment_violations(Molecule* m, float simt)
{
    int i, n = m->get_atom_count();
    float viol = 0;
    for (i=0; i<n; i++)
    {
        Atom* a = m->get_atom(i);
        float f = sphere_inside_pocket(a->get_sphere());
        int Z = a->Z;

        viol += ((Z > 1) ? 1 : 0.5) * (1.0 - f);
        if ((simt >= 0) && (viol > simt)) return viol;
    }

    return viol;
}

float Cavity::find_best_containment(Molecule* m, bool mbt)
{
    ligand_vertices = m->obtain_vdW_surface(10);
    Point cen = get_center();
    m->recenter(cen);
    Pose best(m);
    float bestviol = Avogadro;
    LocatedVector lv;

    Vector axisx = Point(1,0,0);
    Vector axisy = Point(0,1,0);
    Vector axisz = Point(0,0,1);
    int i, j, k, l, n = m->get_atom_count();
    float thx, thy, thz;

    for (thx=0; thx < M_PI*2; thx += cav_360_step)
    {
        for (thy=0; thy < M_PI*2; thy += cav_360_step)
        {
            for (thz=0; thz < M_PI*2; thz += cav_360_step)
            {
                float atoms_outside_cavity = 0;
                std::string ldbg = "";
                for (i=0; i<n; i++)
                {
                    Atom* a = m->get_atom(i);
                    if (a->Z < 2) continue;

                    CPartial* inside;
                    float f = sphere_inside_pocket(a->get_sphere(), &inside);
                    /* if (inside && a->get_charge() && a->distance_to(prot->get_residue(264)->get_atom("HH2")) < 5)
                        cout << a->distance_to(prot->get_residue(264)->get_atom("HH2")) << endl << flush; */
                    if (mbt && inside)
                    {
                        float e = 0;
                        float achg = a->get_charge();
                        if (!achg) achg = a->is_conjugated_to_charge();
                        if (inside->chargedp && achg < -hydrophilicity_cutoff) e = 0.6;
                        else if (inside->chargedn && achg > hydrophilicity_cutoff) e = 0.6;
                        else if (inside->metallic && (a->get_family() == CHALCOGEN || a->get_family() == PNICTOGEN) && a->Z != 8) e = 0.8;
                        else if (inside->polar && fabs(a->is_polar()) > hydrophilicity_cutoff) e = 0.3;
                        else if (inside->thio && a->get_family() == CHALCOGEN && a->Z != 8) e = 0.15;
                        else if (inside->pi && a->is_pi()) e = 0.12;
                        if (e && inside->priority)
                        {
                            /* if (inside->chargedp) cout << inside->s.center << " "
                                << inside->resnos_as_string(prot) << (inside->chargedn ? "-" : "") << (inside->chargedp ? "+" : "")
                                << " " << a->name << " " << achg << " " << e << endl; */
                            f = e;
                        }
                        f += e;
                    }
                    atoms_outside_cavity += (1.0-f);
                }   // for i

                if (atoms_outside_cavity < bestviol)
                {
                    best.copy_state(m);
                    bestviol = atoms_outside_cavity;
                }

                lv = axisz;
                lv.origin = cen;
                m->rotate(lv, cav_360_step);
            }   // for thz

            lv = axisy;
            lv.origin = cen;
            m->rotate(lv, cav_360_step);
        }   // for thy

        lv = axisx;
        lv.origin = cen;
        m->rotate(lv, cav_360_step);
    }   // for thx

    best.restore_state(m);

    for (j=0; j<=26; j++)
    {
        Point maybe = cen;
        k=0;
        _retry_linear_motion:
        maybe.x += 0.5 * (j%3-1);
        maybe.y += 0.5 * ((j/3)%3-1);
        maybe.z += 0.5 * j/9;

        m->recenter(maybe);
        float viol = containment_violations(m) + m->total_eclipses()*33;
        if (viol < bestviol)
        {
            best.copy_state(m);
            bestviol = viol;
            k++;
            if (k < 5) goto _retry_linear_motion;
        }
    }
    best.restore_state(m);

    return bestviol;
}

float CPartial::atom_match_score(Atom* a)
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

float Cavity::match_ligand(Molecule* ligand, Atom** match_atom, CPartial** match_partial, Protein* prot)
{
    CPartial* mpart[10];
    Atom*     matom[10];
    float     score[10];
    int matches = 0;

    int i, j, l, m, n;
    float f, f0;

    for (i=0; i<10; i++) score[i] = 0;

    n = ligand->get_atom_count();
    for (i=0; i<n; i++)
    {
        Atom* a = ligand->get_atom(i);
        for (j=0; j<pallocd; j++)
        {
            CPartial* p = &partials[j];
            if (p->s.radius < min_partial_radius) break;

            f = p->atom_match_score(a);

            for (l=0; l<10 && l<=matches; l++)
            {
                if (f > score[l])
                {
                    for (m=9; m>l; m--)
                    {
                        mpart[m] = mpart[m-1];
                        matom[m] = matom[m-1];
                        score[m] = score[m-1];
                    }
                    mpart[l] = p;
                    matom[l] = a;
                    score[l] = f;
                    matches = max(matches, l+1);
                    break;
                }
            }
        }
    }

    if (!matches)
    {
        if (match_atom) *match_atom = nullptr;
        if (match_partial) *match_partial = nullptr;
        return -Avogadro;
    }

    Pose best_of_all(ligand);
    for (m=0; m<matches; m++)
    {
        // Each match, try centering the atom inside the partial.
        Vector mov = mpart[m]->s.center.subtract(matom[m]->loc);
        ligand->move(mov);

        AminoAcid* reaches[SPHREACH_MAX+4];
        for (i=0; i<SPHREACH_MAX; i++) reaches[i] = nullptr;
        if (prot) resnos(prot, reaches);

        // Do 3-axis rotations to maximize containment.
        Vector axisx = Point(1,0,0);
        Vector axisy = Point(0,1,0);
        Vector axisz = Point(0,0,1);
        Pose best(ligand);
        float bestc = 0;
        LocatedVector lv;
        float thx, thy, thz;

        for (thx=0; thx < M_PI*2; thx += cav_360_step)
        {
            for (thy=0; thy < M_PI*2; thy += cav_360_step)
            {
                for (thz=0; thz < M_PI*2; thz += cav_360_step)
                {
                    f = 0;
                    for (i=0; i<n; i++)
                    {
                        f += sphere_inside_pocket(ligand->get_atom(i)->get_sphere());
                    }

                    if (prot)
                    {
                        f /= n;
                        if (f > 0.5)
                        {
                            f = 0;
                            for (i=0; i<SPHREACH_MAX && reaches[i]; i++)
                            {
                                float f1 = reaches[i]->get_intermol_potential(ligand);
                                if (f1 > f) f = f1;
                            }

                            if (f > bestc)
                            {
                                bestc = f;
                                best.copy_state(ligand);
                            }
                        }
                    }
                    else if (f > bestc)
                    {
                        bestc = f;
                        best.copy_state(ligand);
                    }

                    lv = axisz;
                    lv.origin = matom[m]->loc;
                    ligand->rotate(lv, cav_360_step);
                }   // for thz

                lv = axisy;
                lv.origin = matom[m]->loc;
                ligand->rotate(lv, cav_360_step);
            }   // for thy

            lv = axisx;
            lv.origin = matom[m]->loc;
            ligand->rotate(lv, cav_360_step);
        }   // for thx
        best.restore_state(ligand);

        // Do xyz wiggle to improve containment.
        for (j=0; j<=26; j++)
        {
            Point maybe = ligand->get_barycenter();
            l=0;
            _retry__linear_motion:
            maybe.x += 0.5 * (j%3-1);
            maybe.y += 0.5 * ((j/3)%3-1);
            maybe.z += 0.5 * j/9;

            ligand->recenter(maybe);
            f = 0;
            for (i=0; i<n; i++)
            {
                f += sphere_inside_pocket(ligand->get_atom(i)->get_sphere());
            }

            if (f > bestc)
            {
                best.copy_state(ligand);
                bestc = f;
                l++;
                if (l < 5) goto _retry__linear_motion;
            }
        }
        best.restore_state(ligand);
        if (!m) best_of_all.copy_state(ligand);

        // If no satisfactory containment, go on to next match.
        // Return true if good match found.
        if (bestc >= min_cavmatch_ctainmt)
        {
            if (match_atom) *match_atom = matom[m];
            if (match_partial) *match_partial = mpart[m];
            return bestc;
        }
        if (!m) f0 = bestc;
    }

    if (match_atom) *match_atom = matom[0];
    if (match_partial) *match_partial = mpart[0];
    best_of_all.restore_state(ligand);
    return f0;
}

int Cavity::resnos(Protein* p, AminoAcid** result)
{
    if (!partials || !pallocd) return 0;
    int i, j, n=p->get_end_resno();
    bool resincluded[n+4];
    for (i=0; i<=n; i++) resincluded[i] = false;
    for (i=0; i<pallocd; i++)
    {
        if (!partials[i].s.radius) break;
        std::string pres = partials[i].resnos_as_string(p);             // This is a kludge.
        char buffer[1024];
        strcpy(buffer, pres.c_str());
        char** words = chop_spaced_words(buffer);
        if (!words) continue;
        for (j=0; words[j]; j++) resincluded[atoi(words[j])] = true;
    }

    j = 0;
    for (i=1; i<=n; i++)
    {
        if (resincluded[i])
        {
            result[j++] = p->get_residue(i);
            if (j >= SPHREACH_MAX-1) break;
        }
    }

    return j;
}

std::string Cavity::resnos_as_string(Protein* p)
{
    if (!partials || !pallocd) return (std::string)"";
    int i, j, n=p->get_end_resno();
    bool resincluded[n+4];
    for (i=0; i<=n; i++) resincluded[i] = false;
    for (i=0; i<pallocd; i++)
    {
        if (!partials[i].s.radius) break;
        std::string pres = partials[i].resnos_as_string(p);
        char buffer[1024];
        strcpy(buffer, pres.c_str());
        char** words = chop_spaced_words(buffer);
        if (!words) continue;
        for (j=0; words[j]; j++) resincluded[atoi(words[j])] = true;
    }

    std::string result = "";
    for (i=1; i<=n; i++)
    {
        if (resincluded[i])
        {
            if (result.length()) result += " ";
            result += std::to_string(i);
        }
    }

    return result;
}

int Cavity::resnos_as_array(Protein *p, int *output)
{
    int i, j, l, m, n;

    n = p->get_end_resno();
    bool included[n+1];
    int tmp[256];

    for (i=0; i<n; i++) included[i] = false;

    for (i=0; i<pallocd; i++)
    {
        if (partials[i].s.radius)
        {
            m = partials[i].resnos_as_array(p, tmp);
            for (j=0; j<m; j++)
            {
                included[tmp[j]] = true;
            }
        }
    }

    l=0;
    for (i=1; i<n; i++)
    {
        if (included[i]) output[l++] = i;
    }
    return l;
}

std::string CPartial::resnos_as_string(Protein* p)
{
    int i, j, n = p->get_end_resno();
    bool intersect[n+4];
    for (i=1; i<=n; i++)
    {
        AminoAcid* aa = p->get_residue(i);
        if (!aa)
        {
            intersect[i] = false;
            continue;
        }
        j = aa->atoms_inside_sphere(s, nullptr, aa->priority ? 4.0 : 1.1);
        intersect[i] = (j>0);
    }

    std::string result;
    for (i=1; i<=n; i++) if (intersect[i]) result += (std::string)(result.length() ? " " : "") + std::to_string(i);
    return result;
}

int CPartial::resnos_as_array(Protein *p, int *output)
{
    int i, j, l, n = p->get_end_resno();
    bool intersect[n+4];
    for (i=1; i<=n; i++)
    {
        AminoAcid* aa = p->get_residue(i);
        if (!aa)
        {
            intersect[i] = false;
            continue;
        }
        j = aa->atoms_inside_sphere(s, nullptr, aa->priority ? 4.0 : 1.1);
        intersect[i] = (j>0);
    }

    l=0;
    for (i=1; i<=n; i++) if (intersect[i]) output[l++] = i;
    output[l] = 0;
    return l;
}

int CPartial::from_cvty_line(char* lndata)
{
    int cno;

    //   0    2.212   14.679    9.921   2.180 M-+HSP! 161 200 201 204 264
    char buffer[65536];
    strcpy(buffer, lndata);
    char** words = chop_spaced_words(buffer);

    if (!words || !words[0] || !words[1] || !words[2] || !words[3] || !words[4] || !words[5]) return -1;

    cno = atoi(words[0]);
    s.center.x = atof(words[1]);
    s.center.y = atof(words[2]);
    s.center.z = atof(words[3]);
    s.radius = atof(words[4]);
    chargedn = strchr(words[5], '-');
    chargedp = strchr(words[5], '+');
    polar    = strchr(words[5], 'H');
    thio     = strchr(words[5], 'S');
    pi       = strchr(words[5], 'P');
    priority = strchr(words[5], '!');

    delete[] words;

    return cno;
}

void CPartial::write_cvty_line(char* outdata, int cno, Protein* protein)
{
    sprintf(outdata, "%4d %8.3f %8.3f %8.3f %7.3f %c%c%c%c%c%c%c %s\n", cno, 
        s.center.x,
        s.center.y,
        s.center.z,
        s.radius,
        metallic ? 'M' : '_',
        chargedn ? '-' : '_',
        chargedp ? '+' : '_',
        polar    ? 'H' : '_',
        thio     ? 'S' : '_',
        pi       ? 'P' : '_',
        priority ? '!' : '_',
        protein ? resnos_as_string(protein).c_str() : ""
        );
}

int Cavity::estimate_multiplicity(Molecule* ligand)
{
    // Later we'll examine partial radii and molecular dimensions but for now here's a quick and dirty.

    float lv = ligand->get_volume();
    return this->get_volume() / lv / 3;
}




