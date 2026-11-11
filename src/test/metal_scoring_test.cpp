#include "../classes/scoring.h"
#include <cassert>
#include <cmath>
#include <cstring>
#include <memory>

int main(int argc,char** argv)
{
    assert(argc<=2);
    Protein protein("metal scoring");
    FILE* input=fopen(argc==2?argv[1]:"testdata/metal-site-OR2T11.pdb","r");assert(input);
    protein.load_pdb(input);fclose(input);
    MCoord site;site.Z=29;site.charge=1;site.ncoordres=3;
    site.coordres[0].resno=199;site.coordres[1].resno=200;site.coordres[2].resno=203;
    protein.coordinate_metal(&site,1);
    std::unique_ptr<Molecule> ions(protein.metals_as_molecule());
    assert(ions->get_atom_count()==1 && ions->get_atom(0)==site.mtl);
    Vector* geometry=site.mtl->get_geometry_aligned_to_bonds();
    Point third=protein.get_atom(203,"SD")->loc.subtract(site.mtl->loc);
    Point origin;
    assert(find_3d_angle(third,Point(geometry[2]),origin)
        <=find_3d_angle(third,Point(geometry[3]),origin)+1e-6);
    Molecule disulfide("metal search probe");
    disulfide.from_smiles("CSSC");
    LigandTarget primary;primary.single_atom=disulfide.get_atom(1);
    assert(primary.single_atom->Z==16);
    BestBindingResult targets;
    targets.protein=&protein;targets.ligand=&disulfide;
    targets.pri_res=protein.get_residue(199);targets.pri_tgt=&primary;
    Point ion_position=site.mtl->loc;
    Search::align_targets(&disulfide,ion_position.add(Point(8,4,2)),&targets);
    assert(std::fabs(primary.single_atom->distance_to(site.mtl)-2.29)<.01);
    for(int i=0;i<_INTER_TYPES_LIMIT;++i)total_binding_by_type[i]=0;
    disulfide.get_intermol_binding(ions.get(),false);
    assert(total_binding_by_type[mcoord-covalent]<-1);
    assert(site.mtl->loc.get_3d_distance(ion_position)<1e-6);
    Vector away=primary.single_atom->loc.subtract(ion_position);away.r=2;
    disulfide.move(away);
    float displaced=primary.single_atom->distance_to(site.mtl);
    float bond_length=disulfide.get_atom(0)->distance_to(primary.single_atom);
    Search::align_targets(&disulfide,ion_position,&targets,.3);
    assert(std::fabs(primary.single_atom->distance_to(site.mtl)-(displaced-.3*(displaced-2.29)))<.01);
    assert(std::fabs(disulfide.get_atom(0)->distance_to(primary.single_atom)-bond_length)<.001);
    Molecule ligand("sulfur probe");
    Point position=site.mtl->loc.add(Point(2.29,0,0));
    Atom* sulfur=ligand.add_atom("S","S1",&position,nullptr,0,0);
    Point best=position;float lowest=0;
    for(int x=-1;x<=1;++x)for(int y=-1;y<=1;++y)for(int z=-1;z<=1;++z)
    {
        if(!x && !y && !z)continue;
        Point direction(x,y,z);direction.scale(2.29);
        sulfur->move(site.mtl->loc.add(direction));
        float energy=ligand.get_intermol_binding(ions.get(),false).summed();
        if(energy<lowest){lowest=energy;best=sulfur->loc;}
    }
    assert(lowest<-1);
    sulfur->move(best);
    for(int i=0;i<_INTER_TYPES_LIMIT;++i)total_binding_by_type[i]=0;
    float expected=ligand.get_intermol_binding(ions.get(),false).summed();
    float expected_coord=total_binding_by_type[mcoord-covalent];
    assert(expected_coord<-1);
    AminoAcid* donors[]={protein.get_residue(199),protein.get_residue(200),protein.get_residue(203),nullptr};
    // Ion energy is independent of which donor residues enter the scoring sphere.
    for(int count:{0,1,3})
    {
        AminoAcid* selected[4]={nullptr,nullptr,nullptr,nullptr};
        float residue_energy=0;
        for(int i=0;i<count;++i)
        {
            selected[i]=donors[i];
            residue_energy+=ligand.get_intermol_binding(donors[i],false).summed();
        }
        DockResult result(&protein,&ligand,count,selected);
        assert(std::fabs(result.bytype[mcoord-covalent]-expected_coord)<.01);
        assert(std::fabs(result.kJmol_raw-(expected+residue_energy))<.1);
        int metrics=0;
        for(int i=0;result.metric[i] && result.metric[i][0];++i)
            if(!strcmp(result.metric[i],"Metals"))
            {
                ++metrics;assert(std::fabs(result.mkJmol[i]-expected)<.01);
            }
        assert(metrics==1);
    }
    // A serialized ion still contributes before its donor associations are rebuilt.
    FILE* saved=tmpfile();assert(saved);
    protein.save_pdb(saved);protein.end_pdb(saved);rewind(saved);
    Protein loaded("loaded metal scoring");loaded.load_pdb(saved);fclose(saved);
    assert(loaded.get_metals_count()==1 && !loaded._mcoords[0].mtl);
    std::unique_ptr<Molecule> loaded_ions(loaded.metals_as_molecule());
    for(int i=0;i<_INTER_TYPES_LIMIT;++i)total_binding_by_type[i]=0;
    expected=ligand.get_intermol_binding(loaded_ions.get(),false).summed();
    expected_coord=total_binding_by_type[mcoord-covalent];
    AminoAcid* empty[]={nullptr};
    DockResult result(&loaded,&ligand,0,empty);
    assert(std::fabs(result.kJmol_raw-expected)<.01);
    assert(std::fabs(result.bytype[mcoord-covalent]-expected_coord)<.01);
    std::cout<<"Metal scoring includes each owned ion once, with or without donor associations.\n";
}
