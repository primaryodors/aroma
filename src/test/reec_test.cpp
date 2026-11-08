#define REEC_NO_MAIN
#include "../reec.cpp"
#include <cassert>

int main(int argc, char** argv)
{
    using namespace reec;
    // Similarity, local context, and coverage all affect the same blend.
    Template close, distant;
    std::string blend_sequence(25,'L');
    for(int r=0;r<25;++r)
    {
        close.mapping.push_back(r);
        distant.mapping.push_back(r);
        close.residues.push_back({r+1,'L', {}});
        distant.residues.push_back({r+1,'D', {}});
    }
    TemplateBlend dominant(blend_sequence,{close,distant});
    assert(dominant.weights[0][12]>.99);
    // Conservative substitutions get intermediate support.
    assert(sequence_credit('L','I')>sequence_credit('L','D'));
    assert(sequence_credit('L','I')<sequence_credit('L','L'));
    for(int r=0;r<25;++r)
    {
        close.residues[r].aa=r<12?'L':'D';
        distant.residues[r].aa=r<12?'D':'L';
    }
    TemplateBlend regional(blend_sequence,{close,distant});
    assert(regional.weights[0][3]>regional.weights[1][3]);
    assert(regional.weights[0][21]<regional.weights[1][21]);
    close.mapping[12]=-1;
    TemplateBlend primary(blend_sequence,{close,distant},0);
    assert(primary.weights[0][12]==0&&primary.weights[1][12]==1);
    assert(primary.weights[0][3]==1&&primary.weights[1][3]==0);
    assert(primary.weights[0][11]==.25&&std::fabs(primary.weights[1][11]-.75)<1e-12);
    assert(primary.weights[0][10]==.5&&primary.weights[0][9]==.75);
    assert(primary.pair_weight(1,3,21)==0);
    for(int r=0;r<25;++r)
        assert(std::fabs(primary.weights[0][r]+primary.weights[1][r]-1)<1e-12);
    TemplateBlend reversed(blend_sequence,{distant,close});
    TemplateBlend forward(blend_sequence,{close,distant});
    for(int r=0;r<25;++r)assert(std::fabs(reversed.weights[1][r]-forward.weights[0][r])<1e-12);
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
    // Equal templates average coordinates, while primary mode preserves its input.
    Template displaced_template=fragment;
    displaced_template.id="donor";
    fragment.id="primary";
    for(auto& residue:displaced_template.residues)
        for(auto& atom:residue.atoms)atom.second+=V{.2,.1,0};
    std::vector<Template> blended_templates{fragment,displaced_template};
    Model blended, preserved;
    blended.build(sequence,blended_templates);
    preserved.primary_template="primary";
    preserved.build(sequence,blended_templates);
    for(size_t r=0;r<fragment.residues.size();++r)
    {
        int ca=blended.id(r+1,"CA");
        V original=fragment.residues[r].atoms.at("CA");
        assert((blended.anchors[ca]-original-V{.1,.05,0}).norm()<1e-10);
        assert((preserved.anchors[ca]-original).norm()<1e-10);
    }
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
    XYZ base=backbone.xyz;
    backbone.seed=17;
    backbone.alternate_start(0);
    XYZ seeded=backbone.xyz;
    backbone.xyz=base;
    backbone.alternate_start(0);
    for(size_t i=0;i<seeded.size();++i)assert((seeded[i]-backbone.xyz[i]).norm()<1e-12);
    backbone.xyz=base;
    backbone.seed=18;
    backbone.alternate_start(0);
    double change=0;
    for(size_t i=0;i<seeded.size();++i)change+=(seeded[i]-backbone.xyz[i]).norm();
    assert(change>1e-6);
    for(size_t r=1;r<=sequence.size();++r)backbone.anneal_residues.insert(r);
    auto rotations=backbone.targeted_rotations();
    assert(!rotations.empty());
    XYZ before_rotation=backbone.xyz;
    backbone.rotate(rotations.front(),.75);
    for(auto pair:backbone.geometry)if(pair.bond)
        assert(std::fabs((before_rotation[pair.a]-before_rotation[pair.b]).norm()
            -(backbone.xyz[pair.a]-backbone.xyz[pair.b]).norm())<1e-9);
    for(size_t i=0;i<backbone.nodes.size();++i)if(backbone.nodes[i].backbone)
        assert((before_rotation[i]-backbone.xyz[i]).norm()<1e-12);
    backbone.xyz=before_rotation;
    // Construct a displaced rotamer with a known reference contact.
    double largest_change=0,reference_distance=0;
    int moving_atom=-1,reference_atom=-1;
    XYZ displaced=before_rotation;
    for(const auto& rotation:rotations)
    {
        backbone.xyz=before_rotation;
        backbone.rotate(rotation,1.7);
        for(int atom:rotation.moving)for(size_t anchor=0;anchor<backbone.nodes.size();++anchor)
            if(backbone.nodes[anchor].backbone)
            {
                double original=(before_rotation[atom]-before_rotation[anchor]).norm();
                double current=(backbone.xyz[atom]-backbone.xyz[anchor]).norm();
                if(std::fabs(current-original)>largest_change)
                {
                    largest_change=std::fabs(current-original);
                    reference_distance=original;
                    moving_atom=atom;
                    reference_atom=anchor;
                    displaced=backbone.xyz;
                }
            }
    }
    assert(largest_change>.1);
    backbone.pair(moving_atom,reference_atom,reference_distance,100000);
    before_rotation=displaced;
    backbone.xyz=before_rotation;
    backbone.anneal(64,1000000,0);
    assert(backbone.sampling.trials==64);
    assert(backbone.sampling.accepted>0&&backbone.sampling.uphill>0);
    assert(backbone.sampling.after<backbone.sampling.before);
    XYZ annealed=backbone.xyz;
    backbone.xyz=before_rotation;
    backbone.anneal(64,1000000,0);
    for(size_t i=0;i<annealed.size();++i)assert((annealed[i]-backbone.xyz[i]).norm()<1e-12);
    for(size_t i=0;i<backbone.nodes.size();++i)if(backbone.nodes[i].backbone)
        assert((before_rotation[i]-backbone.xyz[i]).norm()<1e-12);
    backbone.reset_template_release();
    auto original_weights=backbone.anchor_weight;
    std::vector<double> original_geometry;
    for(const auto& pair:backbone.geometry)original_geometry.push_back(pair.weight);
    backbone.refresh(backbone.xyz);
    assert(!backbone.neighbors.empty());
    auto clash=backbone.neighbors.front();
    backbone.xyz[clash.second]=backbone.xyz[clash.first]+V{.1,0,0};
    backbone.release_clashing_templates();
    assert(backbone.released_template_residues.count(backbone.nodes[clash.first].res));
    assert(backbone.anchor_weight[clash.first]==original_weights[clash.first]*.1);
    for(size_t i=0;i<backbone.geometry.size();++i)if(backbone.geometry[i].validate)
        assert(backbone.geometry[i].weight==original_geometry[i]);
    backbone.reset_template_release();
    assert(backbone.released_template_residues.empty());
    assert(backbone.anchor_weight==original_weights);
    for(size_t i=0;i<backbone.geometry.size();++i)assert(backbone.geometry[i].weight==original_geometry[i]);
    std::cout << "REEC analytical gradients and rigid template fit pass.\n";
}
