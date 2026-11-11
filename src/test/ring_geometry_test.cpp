#include "../classes/molecule.h"
#include <cassert>
#include <cmath>
#include <fstream>
#include <iterator>
#include <vector>
#include <string>

int main(int argc,char** argv)
{
    bool probe=argc>1 && std::string(argv[1])=="--probe";
    std::ifstream input("testdata/N-methylpiperidine.sdf");assert(input);
    std::string text((std::istreambuf_iterator<char>(input)),std::istreambuf_iterator<char>());
    Molecule molecule("N-methylpiperidine");
    assert(molecule.from_sdf(text.c_str())==21);
    assert(std::fabs(molecule.get_charge()-1)<.001);
    struct Length{Atom* first;Atom* second;float initial;};
    struct Angle{Atom* center;Atom* first;Atom* second;float initial;};
    struct Stereo{Atom* center;Atom* first;Atom* second;Atom* third;double initial;};
    std::vector<Length> lengths;
    std::vector<Angle> angles;
    std::vector<Stereo> stereocenters;
    auto volume=[](Atom* center,Atom* first,Atom* second,Atom* third)
    {
        Point a=first->loc.subtract(center->loc),b=second->loc.subtract(center->loc),c=third->loc.subtract(center->loc);
        return a.x*(b.y*c.z-b.z*c.y)-a.y*(b.x*c.z-b.z*c.x)+a.z*(b.x*c.y-b.y*c.x);
    };
    for(int i=0;i<molecule.get_atom_count();++i)
    {
        Atom* atom=molecule.get_atom(i);
        std::vector<Atom*> neighbors;
        for(int slot=0;slot<std::abs(atom->get_geometry());++slot)
        {
            Bond* bond=atom->get_bond_by_idx(slot);
            if(bond && bond->atom2 && atom<bond->atom2)
                lengths.push_back({atom,bond->atom2,atom->distance_to(bond->atom2)});
            if(bond && bond->atom2)neighbors.push_back(bond->atom2);
        }
        for(size_t a=0;a<neighbors.size();++a)for(size_t b=a+1;b<neighbors.size();++b)
            angles.push_back({atom,neighbors[a],neighbors[b],find_3d_angle(neighbors[a]->loc,neighbors[b]->loc,atom->loc)});
        if(neighbors.size()==4)
            stereocenters.push_back({atom,neighbors[0],neighbors[1],neighbors[2],volume(atom,neighbors[0],neighbors[1],neighbors[2])});
    }
    auto check=[&](const char* stage)
    {
        float maximum=0;
        float angle_change=0;
        for(auto bond:lengths)maximum=std::max(maximum,std::fabs(bond.first->distance_to(bond.second)-bond.initial));
        for(auto angle:angles)angle_change=std::max(angle_change,std::fabs(find_3d_angle(angle.first->loc,angle.second->loc,angle.center->loc)-angle.initial));
        std::cout<<stage<<": C2-C3 "<<molecule.get_atom("C2")->distance_to(molecule.get_atom("C3"))
                 <<" A, maximum bond change "<<maximum<<" A, angle change "<<angle_change*fiftyseven<<" degrees"<<std::endl;
        if(!probe)
        {
            assert(maximum<.001);
            assert(angle_change<10*fiftyseventh);
            for(auto center:stereocenters)
                assert(center.initial*volume(center.center,center.first,center.second,center.third)>0);
        }
    };
    check("input");
    molecule.optimize();check("preparation");
    molecule.minimize_internal_clashes();check("torsion refinement");
    Ring** rings=molecule.get_atom(0)->get_rings();assert(rings && rings[0]);
    int changed=0;
    for(int round=0;round<10;++round)
    {
        for(int j=0;j<rings[0]->get_atom_count();++j)
        {
            Atom* atom=rings[0]->get_atom(j);Point before=atom->loc;
            rings[0]->flip_atom(atom);
            if(before.get_3d_distance(atom->loc)>.001)++changed;
            if(!probe)check("ring flip");
        }
    }
    check("repeated flips");
    std::cout<<"Nontrivial ring moves: "<<changed<<std::endl;
    if(!probe)assert(changed>0);
}
