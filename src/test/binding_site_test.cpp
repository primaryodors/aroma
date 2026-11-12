#include <cassert>
#include <cstdio>
#include <fstream>
#include <sstream>
#include "../classes/search.h"

static void check_receptor(const char* pdb, const char* cvty, const char* sdf,
                           int buried, int exposed, int alternative)
{
    Protein protein;
    FILE* input = fopen(pdb, "rb");
    assert(input);
    protein.load_pdb(input);
    fclose(input);
    std::vector<Cavity> cavities(16);
    int cavity_count = 0;
    input = fopen(cvty, "rb");
    assert(input);
    char line[4096];
    while (fgets(line, sizeof(line), input))
    {
        CPartial partial;
        const int index = partial.from_cvty_line(line);
        assert(index >= 0 && index < int(cavities.size()));
        cavities[index].add_partial(partial);
        cavity_count = std::max(cavity_count, index+1);
    }
    fclose(input);
    Cavity& cavity = cavities[0];
    const std::vector<int> pocket_anchors = buried == 285
        ? std::vector<int>{111,116,120,123,213,269,285}
        : std::vector<int>{110,204,256};
    std::vector<bool> eligible(cavity_count, true);
    const auto choices = Search::pocket_cavity_candidates(&protein, cavities.data(), cavity_count, pocket_anchors, eligible);
    assert(std::find(choices.begin(), choices.end(), 0) != choices.end());
    if (buried == 285 && cavity_count > 1)
    {
        assert(choices.size() == 1);
        eligible[0] = false;
        const auto constrained = Search::pocket_cavity_candidates(&protein, cavities.data(), cavity_count, pocket_anchors, eligible);
        assert(std::find(constrained.begin(), constrained.end(), 0) == constrained.end());
        assert(!constrained.empty());
    }
    const auto unconstrained = Search::pocket_cavity_candidates(&protein, cavities.data(), cavity_count, {}, {});
    assert(int(unconstrained.size()) == cavity_count);
    Point center = cavity.get_center();
    const float blocked = Search::polar_site_availability(&protein, protein.get_residue(buried), center);
    const float available = Search::polar_site_availability(&protein, protein.get_residue(exposed), center);
    const float alternate = Search::polar_site_availability(&protein, protein.get_residue(alternative), center);
    printf("%s: site %d %.5f; site %d %.5f; site %d %.5f\n",
           pdb, buried, blocked, exposed, available, alternative, alternate);
    assert(available > blocked);
    assert(available > 0.1f);
    assert(alternate > blocked);
    if (buried == 285) assert(available > 25*blocked);
    if (buried == 110)
    {
        AminoAcid* residue = protein.get_residue(buried);
        Atom* oxygen = nullptr;
        for (int i=0; i<residue->get_atom_count(); ++i)
            if (!residue->get_atom(i)->is_backbone && residue->get_atom(i)->Z == 8)
                oxygen = residue->get_atom(i);
        assert(oxygen);
        Atom* hydrogen = oxygen->is_bonded_to("H");
        assert(hydrogen);
        const Point original = hydrogen->loc;
        hydrogen->move(Point(100,100,100));
        const float freed = Search::polar_site_availability(&protein, residue, center);
        hydrogen->move(original);
        assert(freed > blocked*1.3f);
    }

    std::ifstream stream(sdf);
    assert(stream);
    std::ostringstream contents;
    contents << stream.rdbuf();
    Molecule ligand;
    assert(ligand.from_sdf(contents.str().c_str()));
    std::vector<LigandTarget> targets(ligand.get_heavy_atom_count()+8);
    Search::identify_ligand_pairing_targets(&ligand, targets.data(), targets.size());
    // Include the competing polar sites and several plausible nonpolar anchors.
    const std::vector<int> anchors = buried == 285
        ? std::vector<int>{buried, exposed, alternative, 166, 120, 170}
        : std::vector<int>{buried, exposed, alternative, 182, 156, 256};
    std::vector<AminoAcid*> candidates;
    for (int residue : anchors)
    {
        AminoAcid* aa = protein.get_residue(residue);
        assert(aa);
        candidates.push_back(aa);
    }
    candidates.push_back(nullptr);
    // Preserve the problematic priority, so geometry must overcome that bias.
    protein.get_residue(buried)->priority = true;
    protein.get_residue(exposed)->priority = buried == 110;
    assert(Search::pocket_cavity_candidates(&protein, cavities.data(), cavity_count, pocket_anchors, {}) == choices);
    for (int seed=1; seed<=4; ++seed)
    {
        srand(seed);
        BestBindingResult result;
        Search::clear_candidates();
        assert(Search::pair_targets(&protein, &ligand, targets.data(), candidates.data(), center, &result, &cavity));
        assert(result.pri_res);
        printf("Seed %d primary: %s\n", seed, result.pri_res->get_name());
        AminoAcid* residues[] = {result.pri_res, result.sec_res, result.tert_res};
        LigandTarget* assigned[] = {result.pri_tgt, result.sec_tgt, result.tert_tgt};
        bool polar_contact = false;
        for (int pair=0; pair<3; ++pair)
            if (residues[pair] && assigned[pair]
                && fabs(assigned[pair]->polarity()) >= hydrophilicity_cutoff)
            {
                polar_contact = true;
                assert(residues[pair] != protein.get_residue(buried));
            }
        assert(polar_contact);
    }
}

static void check_generic_geometry()
{
    Protein protein;
    assert(protein.add_sequence("SK"));
    AminoAcid* site = protein.get_residue(1);
    AminoAcid* wall = protein.get_residue(2);
    for (int i=0; i<wall->get_atom_count(); ++i)
        wall->get_atom(i)->move(Point(1000+i*10,1000,1000));
    Atom* oxygen = site->get_atom("OG");
    assert(oxygen);
    Point outward = oxygen->loc.subtract(site->get_atom("CB")->loc);
    outward.multiply(5/outward.magnitude());
    Point center = oxygen->loc.add(outward);
    const float open = Search::polar_site_availability(&protein, site, center);
    assert(open > 0.1f);

    // Absolute world coordinates cannot decide eligibility.
    const Point translation(31,-27,19);
    for (int i=0; i<site->get_atom_count(); ++i)
        site->get_atom(i)->move(site->get_atom(i)->loc.add(translation));
    center = center.add(translation);
    assert(fabs(Search::polar_site_availability(&protein, site, center)-open) < 0.001f);

    Point offsets[] = {Point(2.8,0,0), Point(-2.8,0,0), Point(0,2.8,0),
                       Point(0,-2.8,0), Point(0,0,2.8), Point(0,0,-2.8)};
    int heavy = 0;
    for (int i=0; i<wall->get_atom_count(); ++i)
    {
        Atom* atom = wall->get_atom(i);
        atom->move(Point(1000+i*10,1000,1000));
        if (atom->Z > 1 && heavy < 6) atom->move(oxygen->loc.add(offsets[heavy++]));
    }
    assert(heavy == 6);
    site->priority = true;
    assert(Search::polar_site_availability(&protein, site, center) == 0);
    for (int i=0; i<wall->get_atom_count(); ++i)
        wall->get_atom(i)->move(Point(1000+i*10,1000,1000));
    assert(fabs(Search::polar_site_availability(&protein, site, center)-open) < 0.001f);

    Cavity separate;
    SPartial sphere;
    outward.multiply(2);
    sphere.s.center = oxygen->loc.add(outward);
    sphere.s.radius = 3;
    separate.add_partial(sphere);
    site->priority = false;
    assert(separate.resnos_as_array(&protein, nullptr) == 0);
    site->priority = true;
    assert(separate.resnos_as_array(&protein, nullptr) == 0);
}

static void check_nonpolar_fallback()
{
    Protein protein;
    assert(protein.add_sequence("AAA"));
    Molecule ligand;
    std::ifstream input("sdf/beta-ionone.sdf");
    assert(input);
    std::ostringstream contents;
    contents << input.rdbuf();
    assert(ligand.from_sdf(contents.str().c_str()));
    std::vector<LigandTarget> targets(ligand.get_heavy_atom_count()+8);
    assert(Search::identify_ligand_pairing_targets(&ligand, targets.data(), targets.size()) >= 2);
    AminoAcid* residues[] = {protein.get_residue(1), protein.get_residue(2), protein.get_residue(3), nullptr};
    BestBindingResult result;
    // An acceptor-only ligand can still use nonpolar contacts without a receptor donor.
    assert(Search::pair_targets(&protein, &ligand, targets.data(), residues,
        protein.get_region_center(1,3), &result));
    assert(result.pri_res && result.pri_tgt);
}

static void check_declared_center()
{
    Cavity cavity;
    CPartial partial;
    partial.s.radius = 2.5;
    partial.s.center = Point(6,11,6);
    cavity.add_partial(partial);
    // An adjacent connected cleft must not shift an explicitly declared site.
    for (int i=0; i<12; ++i)
    {
        partial.s.center = Point(2,14+i*0.1f,3);
        cavity.add_partial(partial);
    }
    Point declared(6,10.5,5.8);
    Point chosen = Search::pocket_starting_center(&cavity,declared);
    assert(chosen.get_3d_distance(declared) < 0.0001f);
    assert(chosen.get_3d_distance(cavity.get_center()) > 3);
    Point outside(9,9,7);
    chosen = Search::pocket_starting_center(&cavity,outside);
    assert(chosen.get_3d_distance(Point(6,11,6)) < 0.0001f);
    assert(cavity.point_inside_pocket(chosen));
    Cavity empty;
    assert(Search::pocket_starting_center(&empty,declared).get_3d_distance(declared) < 0.0001f);
    assert(Search::pocket_starting_center(nullptr,declared).get_3d_distance(declared) < 0.0001f);
}

int main()
{
    check_declared_center();
    check_generic_geometry();
    check_nonpolar_fallback();
    check_receptor("pdbs/OR10/OR10J1.active.pdb", "pdbs/OR10/OR10J1.active.cvty", "sdf/dimetol.sdf", 285, 213, 209);
    check_receptor("pdbs/OR11/OR11A1.active.pdb", "pdbs/OR11/OR11A1.active.cvty", "sdf/geosmin.sdf", 110, 204, 200);
    puts("Pocket-facing polar contacts outrank buried and occupied sites.");
}
