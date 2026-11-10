#include "../classes/protein.h"
#include <cassert>
#include <cmath>
#include <map>
#include <string>
#include <stdexcept>
#include <vector>
#include <fstream>

int main(int argc,char** argv)
{
    assert(argc==3);
    Protein protein;
    FILE* input=fopen(argv[1],"r");assert(input);
    protein.load_pdb(input);fclose(input);
    std::map<Atom*,Point> backbone;
    struct Length{Atom* a;Atom* b;float distance;};
    std::vector<Length> lengths;
    for(int r=1;r<=protein.get_end_resno();++r)
    {
        AminoAcid* aa=protein.get_residue(r);
        if(!aa)continue;
        for(int j=0;j<aa->get_atom_count();++j)
        {
            Atom* a=aa->get_atom(j);
            if(a->is_backbone)backbone[a]=a->loc;
            if(a->Z<2)continue;
            for(int k=0;k<a->get_bonded_atoms_count();++k)
            {
                Bond* bond=a->get_bond_by_idx(k);
                if(bond && bond->atom2 && bond->atom2->Z>1)
                    lengths.push_back({a,bond->atom2,a->distance_to(bond->atom2)});
            }
        }
    }
    MCoord site;
    site.Z=29;site.charge=1;site.ncoordres=3;
    site.coordres[0].set("C5.42");site.coordres[1].set("C5.43");site.coordres[2].set("M5.46");
    protein.coordinate_metal(&site,1);
    assert(site.mtl && protein.get_metals_count()==1);
    assert(protein._mcoords[0].mtl==site.mtl);
    assert(protein._mcoords[0].mtl_original_location.get_3d_distance(site.mtl->loc)<1e-5);
    for(auto item:backbone)assert(item.first->loc.get_3d_distance(item.second)<1e-5);
    for(auto length:lengths)assert(std::fabs(length.a->distance_to(length.b)-length.distance)<.001);
    for(int j=0;j<3;++j)
    {
        Atom* donor=protein.get_atom(site.coordres[j].resno,j==2?"SD":"SG");
        float target=InteratomicForce::coordination_distance(donor,site.mtl);
        assert(std::fabs(target-(j==2?2.29:2.25))<.001);
        assert(std::fabs(donor->distance_to(site.mtl)-target)<.15);
        if(j<2){assert(!donor->is_bonded_to("H"));assert(donor->get_charge()<-.5);}
    }
    Atom* first=site.mtl;
    protein.coordinate_metal(&site,1);
    assert(site.mtl==first && protein.get_metals_count()==1);
    assert(protein._mcoords[1].mtl==nullptr);
    MCoord invalid=site;
    invalid.coordres[2]=invalid.coordres[1];
    bool rejected=false;
    try{protein.coordinate_metal(&invalid,1);}catch(const std::invalid_argument&){rejected=true;}
    assert(rejected && protein.get_metals_count()==1);
    FILE* output=fopen(argv[2],"w");assert(output);
    protein.save_pdb(output);protein.end_pdb(output);fclose(output);
    Protein loaded;
    input=fopen(argv[2],"r");assert(input);loaded.load_pdb(input);fclose(input);
    assert(loaded.get_metals_count()==1);
    assert(std::fabs(loaded.get_metal(0)->get_charge()-1)<.001);
    Atom* serialized_ion=loaded.get_metal(0);
    MCoord again=site; // Its pointer belongs to the previous Protein instance.
    loaded.coordinate_metal(&again,1);
    assert(again.mtl==serialized_ion);
    assert(loaded.get_metals_count()==1 && loaded._mcoords[0].mtl==again.mtl);
    output=fopen(argv[2],"w");assert(output);loaded.save_pdb(output);loaded.end_pdb(output);fclose(output);
    Protein saved;
    input=fopen(argv[2],"r");assert(input);saved.load_pdb(input);fclose(input);
    output=fopen(argv[2],"w");assert(output);saved.save_pdb(output);saved.end_pdb(output);fclose(output);
    assert(saved.get_metals_count()==1);
    // Repeated ions in separate MODEL blocks remain separate structures.
    std::ifstream saved_file(argv[2]);std::vector<std::string> records;std::string line;
    while(std::getline(saved_file,line))if(line!="TER" && line!="END")records.push_back(line);
    saved_file.close();
    FILE* frames=tmpfile();assert(frames);
    for(int frame=1;frame<=2;++frame)
    {
        fprintf(frames,"MODEL %d\n",frame);
        for(const auto& record:records)fprintf(frames,"%s\n",record.c_str());
        fprintf(frames,"ENDMDL\n");
    }
    rewind(frames);
    Protein framed;
    for(int frame=0;frame<2;++frame)
    {
        assert(framed.load_pdb(frames)>0);
        assert(framed.get_metals_count()==1 && framed.get_residue(199));
    }
    fclose(frames);
    input=fopen(argv[1],"r");assert(input);saved.load_pdb(input);fclose(input);
    assert(saved.get_metals_count()==0 && saved._mcoords[0].mtl==nullptr);
    assert(saved.mcoord_resnos[0]==0);
    Protein impossible;
    input=fopen(argv[1],"r");assert(input);impossible.load_pdb(input);fclose(input);
    AminoAcid* distant=impossible.get_residue(203);
    for(int i=0;i<distant->get_atom_count();++i)
    {
        Atom* a=distant->get_atom(i);a->move(a->loc.add(Point(30,0,0)));
    }
    Atom* sulfur=impossible.get_atom(199,"SG");
    float charge=sulfur->get_charge();Point position=sulfur->loc;
    auto mobility=impossible.get_residue(199)->movability;
    MCoord unreachable=site;unreachable.mtl=nullptr;
    rejected=false;
    try{impossible.coordinate_metal(&unreachable,1);}catch(const std::runtime_error&){rejected=true;}
    assert(rejected && impossible.get_metals_count()==0 && impossible._mcoords[0].mtl==nullptr);
    assert(sulfur->loc.get_3d_distance(position)<1e-5 && std::fabs(sulfur->get_charge()-charge)<.001);
    assert(sulfur->is_bonded_to("H") && impossible.get_residue(199)->movability==mobility);
    std::cout<<"Metal fit preserves backbone, covalent geometry, metadata, and serialized ions.\n";
}
