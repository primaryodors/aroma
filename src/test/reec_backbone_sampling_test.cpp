#define REEC_NO_MAIN
#include "../reec.cpp"
#include <cassert>

int main()
{
    using namespace reec;
    Template reference=pdb("testdata/reec-small.pdb",'A');
    std::string sequence;
    for(size_t i=0;i<reference.residues.size();++i)
    {
        reference.mapping.push_back(i);
        sequence+=reference.residues[i].aa;
    }
    Model model;
    std::vector<Template> templates{reference};
    model.build(sequence,templates);
    model.hard_scale=1000;
    model.anchor_scale=.1;
    model.minimize(3000);
    assert(Model::acceptable(model.stats()));
    XYZ origin=model.xyz;
    // Restricted gradients cannot move the surrounding backbone.
    model.local_mobile.assign(model.nodes.size(),false);
    for(size_t i=0;i<model.nodes.size();++i)
        model.local_mobile[i]=model.nodes[i].res>2&&model.nodes[i].res<9;
    model.minimize(20);
    for(size_t i=0;i<model.nodes.size();++i)if(!model.local_mobile[i])
        assert((model.xyz[i]-origin[i]).norm()==0);
    model.local_mobile.clear();
    model.xyz=origin;
    model.seed=7;
    model.sample_backbone(2,9,4,2,1500);
    assert(model.backbone_sampling.trials==4);
    assert(model.local_mobile.empty());
    assert(Model::acceptable(model.stats()));
    assert(model.backbone_sampling.max_displacement<=.75);
    assert(model.backbone_sampling.closure_error<=.08);
    XYZ first=model.xyz;
    for(size_t i=0;i<model.nodes.size();++i)
    {
        int r=model.nodes[i].res;
        if(r==2||r==9||(model.nodes[i].backbone&&(r<2||r>9)))
            assert((first[i]-origin[i]).norm()==0);
    }
    model.xyz=origin;
    model.sample_backbone(2,9,4,2,1500);
    for(size_t i=0;i<first.size();++i)assert((first[i]-model.xyz[i]).norm()==0);
    std::cout<<"Bounded backbone closure, frozen surroundings, and reproducibility pass.\n";
}
