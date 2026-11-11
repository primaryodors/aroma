#include <cassert>
#include <cmath>
#include <cstdio>
#include "../classes/protein.h"

struct ScoringMolecule : Molecule
{
    using Molecule::cfmol_multibind;
    using Molecule::intermol_bind_for_multimol_dock;
};

struct ScoringResidue : AminoAcid
{
    ScoringResidue() : AminoAcid('K') {}
    using Molecule::intermol_bind_for_multimol_dock;
};

static void close_score(float actual, float expected)
{
    assert(std::isfinite(actual));
    if (std::fabs(actual - expected) >= 0.001f * (1.0f + std::fabs(expected)))
        fprintf(stderr, "Score mismatch: actual %g, expected %g\n", actual, expected);
    assert(std::fabs(actual - expected) < 0.001f * (1.0f + std::fabs(expected)));
}

template<class T>
static Interaction verify_score(T& m, Molecule** neighbors)
{
    Interaction score = ScoringMolecule::cfmol_multibind(&m, neighbors);
    Interaction pairs;
    Molecule** selected = m.mclashables ? m.mclashables : neighbors;
    if (selected) for (int i=0; selected[i]; i++)
        pairs += m.intermol_bind_for_multimol_dock(selected[i], nullptr, false, false);
    float strain = m.get_internal_clashes() + m.total_eclipses();
    if (m.is_residue()) strain -= reinterpret_cast<AminoAcid&>(m).initial_eclipses;
    close_score(score.clash - pairs.clash, strain * iteration_internal_clash_coefficient);
    close_score(score.attractive, pairs.attractive);
    close_score(score.repulsive, pairs.repulsive);
    return score;
}

int main()
{
    srand(12345);
    ScoringMolecule m;
    Molecule far_a, far_b;
    m.from_smiles("CCCC", false);
    far_a.from_smiles("C", false);
    far_b.from_smiles("C", false);
    far_a.recenter(Point(1000, 0, 0));
    far_b.recenter(Point(-1000, 0, 0));
    Molecule* zero[] = {nullptr};
    Molecule* one[] = {&far_a, nullptr};
    Molecule* two[] = {&far_a, &far_b, nullptr};

    Pose original(&m);
    Interaction before = verify_score(m, zero);
    verify_score(m, nullptr);
    verify_score(m, one);
    verify_score(m, two);
    m.mclashables = two;
    verify_score(m, one);
    m.mclashables = nullptr;

    // A changed conformation must replace the cached strain, even with no neighbors.
    m.get_atom("C4")->move(m.get_atom("C1")->loc.add(Point(0.2, 0, 0)));
    Interaction clashing = verify_score(m, zero);
    assert(clashing.clash > before.clash + 1.0f);
    assert(!clashing.improved(before));
    assert(clashing.probability(before) < 0.5f);
    verify_score(m, two);
    original.restore_state(&m);
    close_score(verify_score(m, zero).summed(), before.summed());

    // Bond rotation must refresh eclipse strain and pose restoration must recover it.
    ScoringMolecule ethane;
    ethane.from_smiles("CC", false);
    Pose staggered(&ethane);
    Interaction initial = verify_score(ethane, zero);
    Bond* axis = ethane.get_atom("C1")->get_bond_between(ethane.get_atom("C2"));
    assert(axis);
    axis->rotate(hexagonal);
    Interaction rotated = verify_score(ethane, zero);
    assert(std::fabs(rotated.clash - initial.clash) > 0.001f);
    staggered.restore_state(&ethane);
    close_score(verify_score(ethane, zero).summed(), initial.summed());

    ScoringResidue residue;
    residue.initial_eclipses = 2.0f;
    verify_score(residue, zero);
    verify_score(residue, two);
    puts("Internal scoring passed.");
}
