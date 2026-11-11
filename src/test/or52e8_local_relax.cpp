// Isolated ligand-present experiment using the native REEC objective.
// This executable does not publish or replace canonical models.
#define REEC_NO_MAIN
#include "../reec.cpp"

int main(int argc, char** argv)
{
    using namespace reec;
    try
    {
        if(argc!=7) throw std::runtime_error("Usage: local_relax BASE.pdb START.pdb REFERENCE.pdb OUT.pdb REPORT.json ITERATIONS");
        Template base=pdb(argv[1],'A'), start=pdb(argv[2],'A'), reference=pdb(argv[3],'A');
        std::string sequence;
        for(size_t i=0;i<base.residues.size();++i)
        {
            if(base.residues[i].no!=int(i+1)) throw std::runtime_error("Contiguous full receptor numbering is required");
            sequence+=base.residues[i].aa;
            base.mapping.push_back(i);
        }
        Model model;
        std::vector<Template> templates{base};
        model.build(sequence,templates);
        std::map<std::pair<int,std::string>,std::string> records;
        std::ifstream input(argv[1]);
        std::string line;
        while(std::getline(input,line))
            if(line.compare(0,6,"ATOM  ")==0 && line.size()>=54)
                records[{integer(trim(line.substr(22,4))),trim(line.substr(12,4))}]=line;
        std::map<std::pair<int,std::string>,V> initial;
        for(const auto& residue:start.residues)
            for(const auto& atom:residue.atoms) initial[{residue.no,atom.first}]=atom.second;
        std::map<std::string,V> reference_asn;
        for(const auto& residue:reference.residues)
            for(const auto& atom:residue.atoms) reference_asn[atom.first]=atom.second;
        int mobile_count=0;
        for(size_t i=0;i<model.nodes.size();++i)
        {
            const auto& node=model.nodes[i];
            const auto& residue=base.residues[node.res-1];
            if(!residue.atoms.count(node.name)) throw std::runtime_error("Absent input protein heavy atom");
            V original=residue.atoms.at(node.name);
            model.xyz[i]=original;
            model.anchors[i]=original;
            model.anchor_weight[i]=0;
            bool mobile=(node.res==161 || node.res==208) && !node.backbone;
            if(mobile)
            {
                ++mobile_count;
                model.xyz[i]=initial.at({node.res,node.name});
                // Identical restraints for both starting conformers.
                model.anchor_weight[i]=.03;
                if(node.res==208 && (node.name=="CG" || node.name=="OD1" || node.name=="ND2"))
                {
                    model.anchors[i]=reference_asn.at(node.name);
                    model.anchor_weight[i]=.3;
                }
            }
            else
            {
                // The native fixed-record mechanism zeros these gradients
                // and omits fixed/fixed nonbonded terms. Coordinates remain exact.
                model.fixed_records[int(i)]=records.at({node.res,node.name});
            }
        }
        model.environment(argv[1],RigidFit(),true);
        // Frozen-only terms are constant and contribute no mobile gradient.
        // Remove them so the reported objective describes the local experiment.
        auto fixed=[&](int atom){return model.fixed_records.count(atom)>0;};
        model.geometry.erase(std::remove_if(model.geometry.begin(),model.geometry.end(),
            [&](const Pair& pair){return fixed(pair.a)&&fixed(pair.b);}),model.geometry.end());
        auto fixed_volume=[&](const Volume& volume)
            {return fixed(volume.a)&&fixed(volume.b)&&fixed(volume.c)&&fixed(volume.d);};
        model.volumes.erase(std::remove_if(model.volumes.begin(),model.volumes.end(),fixed_volume),model.volumes.end());
        model.torsions.erase(std::remove_if(model.torsions.begin(),model.torsions.end(),fixed_volume),model.torsions.end());
        // A soft contact hypothesis, shared by both arms, not an acceptance rule.
        Contact contact{208,999,"OD1","O10","contact",2.65,3.15,5,false,
                        {{model.id(208,"OD1"),model.id(999,"O10")}}};
        model.contacts.push_back(contact);
        auto before=model.stats();
        model.minimize(integer(argv[6]));
        auto after=model.stats();
        model.save(argv[4],after);
        std::ofstream report(argv[5]);
        report<<std::setprecision(10)
              <<"{\"diagnostic_only\":true,\"mobile_heavy_atoms\":"<<mobile_count
              <<",\"initial_objective\":"<<before.score<<",\"final_objective\":"<<after.score
              <<",\"initial_local_severe_pairs\":"<<before.severe<<",\"final_local_severe_pairs\":"<<after.severe
              <<",\"final_local_maximum_overlap\":"<<after.overlap_max
              <<",\"initial_mobile_bond_error\":"<<before.bond_max
              <<",\"maximum_mobile_bond_error\":"<<after.bond_max
              <<",\"chirality_errors\":"<<after.chirality
              <<",\"hydroxyl_od1_distance\":"<<(model.xyz[model.id(208,"OD1")]-model.xyz[model.id(999,"O10")]).norm()
              <<"}\n";
        if(!report) throw std::runtime_error("Report write failed");
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<error.what()<<"\n";
        return 1;
    }
}
