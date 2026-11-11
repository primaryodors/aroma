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
        reference.mapping.push_back(i);sequence+=reference.residues[i].aa;
    }
    Model model;
    std::vector<Template> templates{reference};
    model.build(sequence,templates);
    model.local_mobile.assign(model.nodes.size(),true);
    model.local_mobile[0]=false;
    XYZ origin=model.xyz;
    double before=model.stats().score;
    model.conjugate_gradient(100);
    auto stop=model.minimizations.back();
    assert(stop.method=="conjugate_gradient"&&stop.completed<=100);
    assert(std::isfinite(stop.max_gradient)&&stop.after<=stop.before);
    assert(model.stats().score<before);
    assert((model.xyz[0]-origin[0]).norm()==0);
    XYZ result=model.xyz;
    model.xyz=origin;
    model.conjugate_gradient(100);
    for(size_t i=0;i<result.size();++i)assert((model.xyz[i]-result[i]).norm()==0);
    std::cout<<"Conjugate-gradient descent, bounded budget, fixed coordinates, and reproducibility pass.\n";
}
