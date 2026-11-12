#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <vector>
#include "../classes/protein.h"

static std::uint64_t fingerprint = 1469598103934665603ULL;
static void record(int number)
{
    fingerprint ^= static_cast<unsigned>(number);
    fingerprint *= 1099511628211ULL;
}

static void verify(Protein& protein)
{
    protein.set_clashables();
    for (int r=protein.get_start_resno(); r<=protein.get_end_resno(); r++)
    {
        AminoAcid* residue = protein.get_residue(r);
        if (!residue) continue;
        std::vector<Molecule*> expected;
        Molecule** actual = residue->mclashables;
        for (int other=protein.get_start_resno(); other<=protein.get_end_resno(); other++)
        {
            AminoAcid* neighbor = protein.get_residue(other);
            if (neighbor && neighbor != residue && residue->can_reach(neighbor)) expected.push_back(neighbor);
        }
        for (size_t i=0; i<expected.size(); i++)
        {
            assert(actual[i] == expected[i]);
            record(expected[i]->is_residue());
        }
        assert(!actual[expected.size()]);
        expected.push_back(nullptr);
        float reference = residue->get_intermol_clashes(expected.data());
        assert(fabs(reference-residue->get_base_mclashes()) < 0.0001f*(1+fabs(reference)));
    }
    for (int r : {108,150,200})
    for (float radius : {4.0f,8.0f,12.0f})
    for (bool facing : {false,true})
    {
        Point point = protein.get_residue(r)->get_CA_location().add(Point(1,2,-1));
        AminoAcid* actual[1024];
        int count = protein.fetch_residues_near(point,radius,actual,facing);
        int expected = 0;
        for (int other=protein.get_start_resno(); other<=protein.get_end_resno(); other++)
        {
            AminoAcid* residue = protein.get_residue(other);
            if (!residue) continue;
            float reach = radius+residue->get_reach();
            if (residue->get_CA_location().get_3d_distance_squared(point) > reach*reach) continue;
            Atom* nearest = residue->get_nearest_atom(point);
            if (!nearest) continue;
            float distance2 = point.get_3d_distance_squared(nearest->loc);
            if (distance2 > radius*radius) continue;
            Atom* cb = residue->get_CB();
            if (facing && cb)
            {
                float tolerance = residue->get_CA_location().get_3d_distance(cb->loc)*sin(hexagonal);
                if (point.get_3d_distance(cb->loc) > sqrt(distance2)+tolerance) continue;
            }
            assert(actual[expected++] == residue);
            record(other);
        }
        assert(count == expected && !actual[count]);
    }
}

int main()
{
    Protein protein;
    FILE* file = fopen("pdbs/OR1/OR1A1.active.pdb","r");
    assert(file);
    protein.load_pdb(file);
    fclose(file);
    verify(protein);
    Atom* cb = protein.get_residue(150)->get_CB();
    assert(cb);
    cb->move(cb->loc.add(Point(0.4,0.1,-0.2)));
    verify(protein);
    Atom* ca = protein.get_residue(150)->get_CA();
    ca->move(ca->loc.add(Point(25,-5,0)));
    verify(protein);
    Molecule ligand;
    ligand.from_smiles("CCO",false);
    AminoAcid* distant = protein.get_residue(protein.get_end_resno());
    distant->priority = true;
    int additional[] = {distant->get_residue_no(),0};
    int excluded[] = {108,0};
    for (bool ignore_priority : {false,true})
    for (Point point : {protein.get_residue(108)->get_CA_location(),Point(1000,1000,1000)})
    {
        ligand.recenter(point);
        AminoAcid* selected[SPHREACH_MAX+4];
        int count = protein.get_residues_can_clash_ligand(selected,&ligand,point,Point(8,8,8),
            additional,ignore_priority,excluded);
        for (int i=0; i<count; i++)
        {
            assert(selected[i]->get_residue_no() != 108);
            record(selected[i]->get_residue_no());
        }
        assert(!selected[count]);
    }
    printf("Residue grid comparisons passed; fingerprint %llu\n",static_cast<unsigned long long>(fingerprint));
}
