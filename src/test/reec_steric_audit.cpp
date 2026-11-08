// Audit serialized coordinates with REEC radii and bonded-pair exclusions.
#define REEC_NO_MAIN
#include "../reec.cpp"

int main(int argc,char** argv)
{
    using namespace reec;
    try
    {
        if(argc!=2)throw std::runtime_error("Usage: reec_steric_audit MODEL.pdb");
        Template reference=pdb(argv[1],'A');
        std::string sequence;
        for(size_t i=0;i<reference.residues.size();++i)
        {
            if(reference.residues[i].no!=int(i+1))
                throw std::runtime_error("Audit expects contiguous protein numbering starting at 1");
            sequence+=reference.residues[i].aa;
            reference.mapping.push_back(i);
        }
        Model model;
        std::vector<Template> templates{reference};
        std::ostringstream construction;
        auto previous=std::cout.rdbuf(construction.rdbuf());
        try
        {
            model.build(sequence,templates);
        }
        catch(...)
        {
            std::cout.rdbuf(previous);
            throw;
        }
        std::cout.rdbuf(previous);
        std::ifstream serialized(argv[1]);
        std::string line;
        bool ligand=false;
        while(std::getline(serialized,line))
            if(line.size()>=54&&line.compare(0,6,"HETATM")==0&&!letter(line.substr(17,3)))
            {
                std::string element=line.size()>=78?trim(line.substr(76,2)):trim(line.substr(12,4)).substr(0,1);
                if(element!="H"&&element!="D")ligand=true;
            }
        if(ligand)model.environment(argv[1]);
        std::vector<bool> present(model.nodes.size(),true);
        int absent=0;
        for(size_t i=0;i<model.nodes.size();++i)
        {
            if(model.fixed_records.count(i))continue;
            const auto& node=model.nodes[i];
            const auto& residue=reference.residues[node.res-1];
            auto position=residue.atoms.find(node.name);
            if(position==residue.atoms.end())
            {
                present[i]=false;
                ++absent;
            }
            else model.xyz[i]=position->second;
        }
        int severe=0;
        double maximum=0;
        std::ostringstream pairs;
        for(size_t i=0;i<model.nodes.size();++i)if(present[i])
            for(size_t j=i+1;j<model.nodes.size();++j)if(present[j])
            {
                if(model.fixed_records.count(i)&&model.fixed_records.count(j))continue;
                if(model.excluded.count({int(i),int(j)}))continue;
                double overlap=model.nodes[i].radius+model.nodes[j].radius-(model.xyz[i]-model.xyz[j]).norm();
                maximum=std::max(maximum,overlap);
                if(overlap<=1.)continue;
                ++severe;
                if(severe>30)continue;
                if(severe>1)pairs<<",";
                pairs<<"{\"res1\":"<<model.nodes[i].res<<",\"atom1\":"<<std::quoted(model.nodes[i].name)
                     <<",\"res2\":"<<model.nodes[j].res<<",\"atom2\":"<<std::quoted(model.nodes[j].name)
                     <<",\"overlap_angstrom\":"<<overlap<<"}";
            }
        std::cout<<"{\"coordinate_source\":\"serialized PDB without minimization\","
                 <<"\"severe_clashes\":"<<severe<<",\"max_overlap_angstrom\":"<<maximum
                 <<",\"unmatched_constructed_atoms\":"<<absent
                 <<",\"fixed_heavy_atoms\":"<<model.fixed_records.size()
                 <<",\"severe_pairs_first_30\":["<<pairs.str()<<"]}\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<"Steric audit: "<<error.what()<<"\n";
        return 1;
    }
}
