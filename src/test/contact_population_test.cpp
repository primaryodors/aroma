#include "../classes/scoring.h"
#include <cassert>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>

int main()
{
    assert(fabs(Protein::contact_population(0)-.5)<1e-6);
    assert(Protein::contact_population(1)>.40 && Protein::contact_population(1)<.41);
    assert(Protein::contact_population(1000)<1e-20);
    assert(Protein::contact_population(-1000)>.999);
    Protein protein("OR52I2");
    Molecule ligand("benzyl_salicylate"),positions;
    FILE* file=fopen("testdata/serine_competing_contact.pdb","r"); assert(file);
    protein.load_pdb(file); rewind(file); positions.from_pdb(file,true); fclose(file);
    std::ifstream sdf("testdata/benzyl_salicylate_contact.sdf");
    std::stringstream text; text<<sdf.rdbuf(); ligand.from_sdf(text.str().c_str());
    ligand.set_name("benzyl_salicylate");
    std::vector<Point> original;
    for (int i=0;i<ligand.get_atom_count();i++)
    {
        Atom* atom=ligand.get_atom(i); Atom* source=positions.get_atom(atom->name); assert(source);
        atom->move(source->loc); original.push_back(atom->loc);
    }
    AminoAcid* ser=protein.get_residue(233); assert(ser);
    ser->movability=MOV_FLEXONLY;
    Point ca=ser->get_atom("CA")->loc, cb=ser->get_atom("CB")->loc;
    float co=ser->get_atom("CB")->distance_to(ser->get_atom("OG"));
    printf("Typed Ser HG to O10 target %.3f A\n",InteratomicForce::hydrogen_bond_distance(ser->get_atom("HG"),ligand.get_atom("O10")));
    assert(protein.refine_hbond_contact(&ligand,ser,false));
    Atom* h=ser->get_atom("HG"); Atom* oxygen=ligand.get_atom("O10");
    printf("Selected contact H-O %.3f A, angle %.1f, population %.4f, gap %.3f kJ/mol\n",
        h->distance_to(oxygen),find_3d_angle(ser->get_atom("OG")->loc,oxygen->loc,h->loc)*fiftyseven,
        ser->contact_population,ser->contact_energy_gap);
    assert(h->distance_to(oxygen)<2.4);
    assert(find_3d_angle(ser->get_atom("OG")->loc,oxygen->loc,h->loc)*fiftyseven>135);
    assert(ser->contact_partner==&ligand && ser->contact_population>=.1 && ser->contact_population<.95);
    assert(ser->get_atom("CA")->loc.get_3d_distance(ca)<1e-5);
    assert(ser->get_atom("CB")->loc.get_3d_distance(cb)<1e-5);
    assert(fabs(ser->get_atom("CB")->distance_to(ser->get_atom("OG"))-co)<1e-5);
    for (int i=0;i<ligand.get_atom_count();i++) if (ligand.get_atom(i)->Z>1)
        assert(ligand.get_atom(i)->loc.get_3d_distance(original[i])<1e-5);
    AminoAcid* residues[]={ser,nullptr};
    float fraction=ser->contact_population;
    ser->contact_population=1;
    DockResult full(&protein,&ligand,1,residues);
    ser->contact_population=fraction;
    DockResult weighted(&protein,&ligand,1,residues);
    assert(full.bytype[hbond]<-1);
    assert(fabs(weighted.bytype[hbond]-fraction*full.bytype[hbond])<.01);
    assert(weighted.contact_population_report.find("Sampled contact population")!=std::string::npos);
    assert(weighted.kJmol_raw>full.kJmol_raw);
    assert(hbond_population_matches(&ligand,ser));
    Atom* moved=ligand.get_atom("C2"); Point saved=moved->loc;
    moved->move(saved.add(Point(.1,0,0)));
    assert(!hbond_population_matches(&ligand,ser));
    moved->move(saved);
    Protein acid_protein("OR52E8");
    Molecule acid("hydroxy_acid"),acid_positions;
    file=fopen("testdata/asparagine_inward_contact.pdb","r"); assert(file);
    acid_protein.load_pdb(file); rewind(file); acid_positions.from_pdb(file,true); fclose(file);
    std::ifstream acid_sdf("testdata/hydroxy_acid_contact.sdf");
    std::stringstream acid_text; acid_text<<acid_sdf.rdbuf(); acid.from_sdf(acid_text.str().c_str());
    for (int i=0;i<acid.get_atom_count();i++)
        acid.get_atom(i)->move(acid_positions.get_atom(acid.get_atom(i)->name)->loc);
    Atom* anchor=acid.get_atom("O9"); Atom* arginine=acid_protein.get_atom(268,"NH1"); assert(arginine);
    acid.stay_close_mine=anchor; acid.stay_close_other=arginine;
    float anchor_distance=anchor->distance_to(arginine);
    std::vector<float> lengths;
    std::vector<std::pair<Atom*,Atom*>> pairs;
    for (int i=0;i<acid.get_atom_count();i++) for (int j=i+1;j<acid.get_atom_count();j++)
        if (acid.get_atom(i)->is_bonded_to(acid.get_atom(j)))
        { pairs.push_back({acid.get_atom(i),acid.get_atom(j)}); lengths.push_back(pairs.back().first->distance_to(pairs.back().second)); }
    AminoAcid* asn=acid_protein.get_residue(208); asn->movability=MOV_FLEXONLY;
    acid.movability=MOV_FLEXONLY;
    acid_protein.refine_hbond_contact(&acid,asn,true);
    assert(anchor->distance_to(arginine)<=anchor_distance+.351);
    for (size_t i=0;i<pairs.size();i++) assert(fabs(pairs[i].first->distance_to(pairs[i].second)-lengths[i])<.002);
    puts("Contact selection, geometry, Boltzmann fractions, and scored partial credit pass.");
}
