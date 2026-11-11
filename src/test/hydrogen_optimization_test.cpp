#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdio>
#include "../classes/protein.h"

static void load(Protein& protein, const char* path)
{
    FILE* file = fopen(path,"r");
    assert(file);
    protein.load_pdb(file);
    fclose(file);
}

static void same_positions(Protein& reference, Protein& optimized)
{
    for (int r=reference.get_start_resno(); r<=reference.get_end_resno(); r++)
    {
        AminoAcid* first = reference.get_residue(r);
        AminoAcid* second = optimized.get_residue(r);
        assert(bool(first) == bool(second));
        if (!first) continue;
        assert(first->get_atom_count() == second->get_atom_count());
        assert(first->movability == second->movability);
        for (int i=0; i<first->get_atom_count(); i++)
        {
            Atom* a = first->get_atom(i);
            Atom* b = second->get_atom(a->name);
            assert(b);
            float error = a->distance_to(b);
            if (error >= 0.0001f)
                fprintf(stderr,"Hydrogen optimization diverged at residue %d atom %s: %g A\n",r,a->name,error);
            assert(error < 0.0001f);
        }
    }
}

static double optimize(Protein& protein, bool cached, int* forced = nullptr)
{
    srand(11973);
    auto start = std::chrono::steady_clock::now();
    protein.optimize_hydrogens(1, INT_MAX, forced, cached);
    return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
}

static void close_value(double actual, double expected)
{
    assert(std::isfinite(actual));
    assert(fabs(actual-expected) < 0.0001*(1+fabs(expected)));
}

static void compare_sweeps(Protein& reference, Protein& optimized)
{
    const char* axes[][2] = {{"CB","OG"},{"CB","OG1"},{"CB","SG"},{"CZ","OH"},{"CE","NZ"},{"CG","CD1"}};
    int checked = 0;
    for (const auto& names : axes)
    {
        for (int r=reference.get_start_resno(); r<=reference.get_end_resno(); r++)
        {
            AminoAcid* first = reference.get_residue(r);
            AminoAcid* second = optimized.get_residue(r);
            if (!first || !second) continue;
            Atom* a = first->get_atom(names[0]);
            Atom* b = first->get_atom(names[1]);
            if (!a || !b) continue;
            Bond* axis = a->get_bond_between(b);
            if (!axis || !axis->can_rotate) continue;
            Bond* other_axis = second->get_atom(names[0])->get_bond_between(second->get_atom(names[1]));
            Pose saved_first(first), saved_second(second);
            AminoAcid* near_first[SPHREACH_MAX+4];
            AminoAcid* near_second[SPHREACH_MAX+4];
            reference.get_residues_can_clash_ligand(near_first,first,first->get_barycenter(),Point(8,8,8),nullptr);
            optimized.get_residues_can_clash_ligand(near_second,second,second->get_barycenter(),Point(8,8,8),nullptr);
            first->get_intermol_binding(near_first);
            second->get_intermol_binding(near_second);
            {
                HydrogenRotationCache cache(*second,other_axis,near_second);
                for (int angle=0; angle<10; angle++)
                {
                    axis->rotate(M_PI/5);
                    other_axis->rotate(M_PI/5);
                    for (float& value : total_binding_by_type) value = 0;
                    missed_connection = Vector();
                    mc_bpotential = 0;
                    Interaction expected = first->get_intermol_binding(near_first);
                    float types[_INTER_TYPES_LIMIT];
                    for (int i=0; i<_INTER_TYPES_LIMIT; i++) types[i] = total_binding_by_type[i];
                    Vector missed = missed_connection;
                    float potential = mc_bpotential;
                    for (float& value : total_binding_by_type) value = 0;
                    missed_connection = Vector();
                    mc_bpotential = 0;
                    Interaction actual = second->get_intermol_binding(near_second);
                    close_value(actual.attractive,expected.attractive);
                    close_value(actual.repulsive,expected.repulsive);
                    close_value(actual.clash,expected.clash);
                    close_value(actual.worst_atom_clash,expected.worst_atom_clash);
                    assert(actual.stays_met == expected.stays_met);
                    for (int i=0; i<_INTER_TYPES_LIMIT; i++) close_value(total_binding_by_type[i],types[i]);
                    close_value(missed_connection.r,missed.r);
                    close_value(missed_connection.theta,missed.theta);
                    close_value(missed_connection.phi,missed.phi);
                    close_value(mc_bpotential,potential);
                }
            }
            saved_first.restore_state(first,true);
            saved_second.restore_state(second,true);
            checked++;
            break;
        }
    }
    assert(checked >= 3);
    puts("Rotation scores and force accounting match the reference.");
}

int main(int argc, char** argv)
{
    Protein reference, optimized;
    const char* path = argc>1 ? argv[1] : "pdbs/OR51/OR51E2.active.pdb";
    srand(771);
    load(reference,path);
    srand(771);
    load(optimized,path);
    same_positions(reference,optimized);
    compare_sweeps(reference,optimized);
    reference.get_residue(150)->movability = MOV_PINNED;
    optimized.get_residue(150)->movability = MOV_PINNED;
    Pose pinned(reference.get_residue(150));
    double baseline = optimize(reference,false);
    double accelerated = optimize(optimized,true);
    same_positions(reference,optimized);
    assert(pinned.total_atom_motions() < 0.0001f);
    printf("Hydrogen optimization matches reference: %.3fs vs %.3fs (%.2fx)\n",
        baseline,accelerated,baseline/accelerated);

    // Forced side-chain rotations and residue flips retain the full scoring path.
    int forced[] = {108,150,155,0};
    optimize(reference,false,forced);
    optimize(optimized,true,forced);
    same_positions(reference,optimized);
    assert(pinned.total_atom_motions() < 0.0001f);
    puts("Forced rotations and pinned residue checks passed.");
}
