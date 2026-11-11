#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>
#include "../classes/protein.h"

// Reference calculation from the original shielding path, including angle side effects.
static bool reference_shielded(const Molecule& m, Atom* a, Atom* b)
{
    if (!a || !b) return false;
    float r = a->distance_to(b);
    float va = a->vdW_radius, vb = b->vdW_radius;
    if (r < va + vb || r < 2) return false;
    float r6 = r*1.26f, r125 = 1.25f*r, r6_2 = r6*r6;
    a->shielding_angle = b->shielding_angle = 0;
    Point aloc = a->loc, bloc = b->loc;
    float low[] = {float(fmin(aloc.x, bloc.x)-r6), float(fmin(aloc.y, bloc.y)-r6), float(fmin(aloc.z, bloc.z)-r6)};
    float high[] = {float(fmax(aloc.x, bloc.x)+r6), float(fmax(aloc.y, bloc.y)+r6), float(fmax(aloc.z, bloc.z)+r6)};
    for (int i=0; i<m.get_atom_count(); i++)
    {
        Atom* ai = m.get_atom(i);
        if (!ai) break;
        if (ai == a || ai == b) continue;
        if (ai->loc.x < low[0] || ai->loc.x > high[0] ||
            ai->loc.y < low[1] || ai->loc.y > high[1] ||
            ai->loc.z < low[2] || ai->loc.z > high[2]) continue;
        float rai2 = ai->loc.get_3d_distance_squared(aloc);
        if (rai2 > r6_2) continue;
        float rbi2 = ai->loc.get_3d_distance_squared(bloc);
        if (rbi2 > r6_2) continue;
        float rai = sqrt(rai2), rbi = sqrt(rbi2);
        if (rai+rbi > r125) continue;
        float vs = ai->vdW_radius;
        if (rai < va+vs || rbi < vb+vs) return true;
        Point sloc = ai->loc;
        float angle = find_3d_angle(&aloc, &bloc, &sloc);
        if (angle > a->shielding_angle) a->shielding_angle = b->shielding_angle = angle;
        if (angle > _shield_angle) return true;
    }
    return false;
}

static void compare(const Molecule& m, Atom* a, Atom* b)
{
    a->shielding_angle = 0.25f;
    b->shielding_angle = 0.75f;
    bool expected = reference_shielded(m, a, b);
    float first = a->shielding_angle, second = b->shielding_angle;
    a->shielding_angle = 0.25f;
    b->shielding_angle = 0.75f;
    bool actual = m.shielded(a, b);
    assert(actual == expected);
    assert(fabs(a->shielding_angle-first) < 0.000001f);
    assert(fabs(b->shielding_angle-second) < 0.000001f);
}

struct Query
{
    const Molecule* molecule;
    Atom* a;
    Atom* b;
};

static std::vector<Query> receptor_queries(Protein& protein, Molecule& ligand)
{
    AminoAcid* center = protein.get_residue(108);
    assert(center);
    ligand.recenter(center->get_CA_location().add(Point(4,0,0)));
    std::vector<Query> queries;
    for (int r=protein.get_start_resno(); r<=protein.get_end_resno(); r++)
    {
        AminoAcid* residue = protein.get_residue(r);
        if (!residue) continue;
        for (int i=0; i<ligand.get_atom_count(); i++)
        {
            Atom* a = ligand.get_atom(i);
            for (int j=0; j<residue->get_atom_count(); j++)
            {
                Atom* b = residue->get_atom(j);
                if (a->distance_to(b) >= _INTERA_R_CUTOFF) continue;
                queries.push_back({&ligand,a,b});
                queries.push_back({residue,a,b});
            }
        }
    }
    assert(!queries.empty());
    return queries;
}

static volatile double benchmark_total = 0;

static double measure(const std::vector<Query>& queries, int rounds, bool optimized, bool refresh = false)
{
    double total = 0;
    auto start = std::chrono::steady_clock::now();
    for (int round=0; round<rounds; round++)
    {
        if (refresh) queries[0].a->move(queries[0].a->loc);
        for (const Query& q : queries)
        {
            bool result = optimized ? q.molecule->shielded(q.a,q.b) : reference_shielded(*q.molecule,q.a,q.b);
            total += result + q.a->shielding_angle + q.b->shielding_angle;
        }
    }
    benchmark_total = total;
    return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
}

int main(int argc, char** argv)
{
    srand(731);
    Molecule endpoints;
    endpoints.from_smiles("CO", false);
    Atom* a = endpoints.get_atom(0);
    Atom* b = endpoints.get_atom(1);
    std::mt19937 random(731);
    std::uniform_real_distribution<double> coordinate(-20,20);
    for (int count : {1,8,32,63,64,65,256,1024})
    {
        Molecule cloud;
        for (int i=0; i<count; i++)
        {
            char name[24];
            snprintf(name, sizeof(name), "X%d", i);
            Point point(coordinate(random),coordinate(random),coordinate(random));
            cloud.add_atom(i%3 ? "C" : "H", name, &point, nullptr, 0);
        }
        for (int i=0; i<2000; i++)
        {
            Point first(coordinate(random),coordinate(random),coordinate(random));
            a->move(first);
            b->move(first.add(Point(coordinate(random)*0.5,coordinate(random)*0.5,coordinate(random)*0.5)));
            compare(cloud,a,b);
            compare(cloud,b,a);
        }
        Pose original(&cloud);
        Point midpoint((a->loc.x+b->loc.x)/2,(a->loc.y+b->loc.y)/2,(a->loc.z+b->loc.z)/2);
        cloud.get_atom(0)->move(midpoint);
        compare(cloud,a,b);
        original.restore_state(&cloud);
        compare(cloud,a,b);
        Atom* added = cloud.add_atom("C", "New", &midpoint, nullptr, 0);
        compare(cloud,a,b);
        cloud.delete_atom(added);
        compare(cloud,a,b);
        compare(cloud,cloud.get_atom(0),b);
        Atom* removed = cloud.get_atom(count/2);
        cloud.delete_atom(removed);
        compare(cloud,a,b);
        delete removed;
        assert(!cloud.shielded(nullptr,b));
    }

    // Distance-sum boundaries and large coordinate offsets must survive spatial culling.
    Molecule edge;
    Point point;
    Atom* blocker = edge.add_atom("C", "Edge", &point, nullptr, 0);
    for (double offset : {0.0,10000.0,-10000.0})
    {
        a->move(Point(offset-4,offset,offset));
        b->move(Point(offset+4,offset,offset));
        for (double delta : {-0.00001,0.0,0.00001})
        {
            blocker->move(Point(offset,offset+3+delta,offset));
            compare(edge,a,b);
        }
    }
    b->move(a->loc);
    compare(edge,a,b);

    Molecule chain;
    chain.from_smiles("CCCCCCCCCCCCCCCCCCCCCCCC",false);
    a->move(Point(0,0,0));
    b->move(Point(8,0,0));
    compare(chain,a,b);
    chain.dehydrogenate();
    assert(chain.get_atom_count() == 24);
    compare(chain,a,b);

    Protein protein;
    const char* path = argc>2 ? argv[2] : "pdbs/OR51/OR51E2.active.pdb";
    FILE* input = fopen(path,"r");
    assert(input);
    protein.load_pdb(input);
    fclose(input);
    Molecule ligand;
    ligand.from_smiles("CC1=CC[C@@H](CC1)C(=C)C", false);
    auto queries = receptor_queries(protein,ligand);
    for (const Query& q : queries) compare(*q.molecule,q.a,q.b);
    Pose pose(&ligand);
    ligand.move(Point(0.5,0.25,-0.5));
    for (const Query& q : queries) compare(*q.molecule,q.a,q.b);
    pose.restore_state(&ligand);
    for (const Query& q : queries) compare(*q.molecule,q.a,q.b);
    Bond** bonds = ligand.get_rotatable_bonds();
    assert(bonds);
    for (int i=0; bonds[i]; i++) if (bonds[i]->can_rotate && !bonds[i]->is_ring_bond())
    {
        bonds[i]->rotate(0.23f);
        for (const Query& q : queries) compare(*q.molecule,q.a,q.b);
        pose.restore_state(&ligand);
        for (const Query& q : queries) compare(*q.molecule,q.a,q.b);
        break;
    }

    Molecule* neighbors[SPHREACH_MAX+4];
    protein.get_residues_can_clash_ligand(reinterpret_cast<AminoAcid**>(neighbors), &ligand,
        ligand.get_barycenter(), Point(8,8,8), nullptr);
    Interaction original_energy = ligand.get_intermol_binding(neighbors);
    ligand.move(Point(0.5,0.25,-0.5));
    ligand.get_intermol_binding(neighbors);
    pose.restore_state(&ligand);
    Interaction restored_energy = ligand.get_intermol_binding(neighbors);
    assert(std::isfinite(restored_energy.summed()));
    assert(fabs(original_energy.summed()-restored_energy.summed()) < 0.00001f*(1+fabs(original_energy.summed())));
    puts("Shielding decisions and angles match the reference.");
    if (argc>1 && !strcmp(argv[1],"--bench"))
    {
        measure(queries,10,true);
        double reference = measure(queries,1000,false);
        double optimized = measure(queries,1000,true);
        printf("Receptor shielding: %zu queries/round; reference %.3fs, optimized %.3fs, %.2fx\n",
            queries.size(),reference,optimized,reference/optimized);
        reference = measure(queries,1000,false,true);
        optimized = measure(queries,1000,true,true);
        printf("With geometry refresh each round: reference %.3fs, optimized %.3fs, %.2fx\n",
            reference,optimized,reference/optimized);
        auto start = std::chrono::steady_clock::now();
        for (int round=0; round<1000; round++)
            benchmark_total = ligand.get_intermol_binding(neighbors).summed();
        printf("Binding smoke benchmark: 1000 evaluations in %.3fs\n",
            std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count());
    }
}
