#define REEC_NO_MAIN
#include "../reec.cpp"
#include <cassert>

int main(int argc, char** argv)
{
    using namespace reec;
    if(argc==2)
    {
        Protein protein;
        FILE* file=fopen(argv[1],"r");assert(file);
        protein.load_pdb(file);fclose(file);
        int links=0;
        for(int a=1;a<=protein.get_end_resno();++a)
            for(int b=a+1;b<=protein.get_end_resno();++b)
            {
                Atom* first=protein.get_atom(a,"SG"),*second=protein.get_atom(b,"SG");
                if(!first||!second||first->distance_to(second)>2.3)continue;
                assert(first->get_bond_between(second));
                protein.get_residue(a)->hydrogenate();
                protein.get_residue(b)->hydrogenate();
                assert(!first->is_bonded_to("H")&&!second->is_bonded_to("H"));
                ++links;
            }
        assert(links>0);
        std::cout<<"REEC disulfides reload as covalent bonds.\n";
        return 0;
    }
    Model model;
    model.xyz = {{0,0,0}, {1.5,.1,0}, {.2,1.4,.2}, {.1,.3,1.6}};
    model.anchors = model.xyz;
    model.anchor_weight.assign(4, .3);
    model.geometry.push_back({0,1,1.2,1.2,5,true});
    model.volumes.push_back({0,1,2,3,2.0,3});
    model.torsions.push_back({0,1,2,3,0,.7});
    model.dihedral_restraints.push_back({0,1,2,3,.8,2.4});
    Contact contact{1,2,"CA","CA","contact",2.5,3.0,4,true,{{1,2}}};
    model.contacts.push_back(contact);
    XYZ gradient;
    model.energy(model.xyz, gradient);
    const double epsilon = 1e-6;
    for (size_t i=0; i<model.xyz.size(); ++i)
    {
        for (int axis=0; axis<3; ++axis)
        {
            double* component = axis==0 ? &model.xyz[i].x : axis==1 ? &model.xyz[i].y : &model.xyz[i].z;
            *component += epsilon;
            XYZ scratch;
            double plus = model.energy(model.xyz, scratch);
            *component -= 2*epsilon;
            double minus = model.energy(model.xyz, scratch);
            *component += epsilon;
            double analytical = axis==0 ? gradient[i].x : axis==1 ? gradient[i].y : gradient[i].z;
            assert(std::fabs((plus-minus)/(2*epsilon)-analytical) < 1e-5);
        }
    }
    Template reference, moving;
    XYZ points{{1,2,3},{4,3,2},{-1,3,7},{2,8,1}};
    for (size_t i=0; i<points.size(); ++i)
    {
        V p=points[i];
        reference.residues.push_back({int(i+1),'A',{{"CA",p}}});
        moving.residues.push_back({int(i+1),'A',{{"CA",{-p.y+8,p.x-4,p.z+3}}}});
        reference.mapping.push_back(i);
        moving.mapping.push_back(i);
    }
    superpose(moving, reference);
    for (size_t i=0; i<points.size(); ++i)
        assert((moving.residues[i].atoms.at("CA")-points[i]).norm() < 1e-8);
    Template fragment=pdb("testdata/reec-small.pdb",'A');
    std::string sequence;
    for(size_t i=0;i<fragment.residues.size();++i)
    {
        sequence+=fragment.residues[i].aa;
        fragment.mapping.push_back(i);
    }
    std::vector<Template> templates{fragment};
    Model backbone;
    backbone.build(sequence,templates);
    bool peptide_bond=false;
    for(auto pair:backbone.geometry)
    {
        const Node& first=backbone.nodes[pair.a];
        const Node& second=backbone.nodes[pair.b];
        // A fixed CB-O distance would incorrectly freeze the backbone psi torsion.
        assert(!(first.res==second.res&&((first.name=="CB"&&second.name=="O")
            ||(first.name=="O"&&second.name=="CB"))));
        if(pair.bond&&first.res+1==second.res&&first.name=="C"&&second.name=="N")
            peptide_bond=true;
    }
    assert(peptide_bond);
    std::cout << "REEC analytical gradients and rigid template fit pass.\n";
}
