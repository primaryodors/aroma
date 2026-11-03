#define REEC_NO_MAIN
#include "../reec.cpp"
#include <cassert>

int main()
{
    using namespace reec;
    Model model;
    model.xyz = {{0,0,0}, {1.5,.1,0}, {.2,1.4,.2}, {.1,.3,1.6}};
    model.anchors = model.xyz;
    model.anchor_weight.assign(4, .3);
    model.geometry.push_back({0,1,1.2,1.2,5,true});
    model.volumes.push_back({0,1,2,3,2.0,3});
    model.torsions.push_back({0,1,2,3,0,.7});
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
    std::cout << "REEC analytical gradients and rigid template fit pass.\n";
}
