#include <fstream>
#include <iostream>
#include <iterator>
#include <cstdio>
#include <stdexcept>
#include "../classes/protein.h"

// Verify the actual reaction after repair removes the adduct and restores the amine.
int main(int argc,char** argv)
{
    try
    {
        if(argc!=4)throw std::runtime_error("Expected receptor PDB, ligand SDF, and lysine number");
        Protein protein;
        FILE* file=std::fopen(argv[1],"r");
        if(!file)throw std::runtime_error("Cannot open receptor");
        protein.load_pdb(file);
        std::fclose(file);
        AminoAcid* lysine=protein.get_residue(std::stoi(argv[3]));
        if(!lysine)throw std::runtime_error("Missing lysine");
        std::ifstream input(argv[2]);
        std::string sdf((std::istreambuf_iterator<char>(input)),std::istreambuf_iterator<char>());
        Molecule ligand("fresh ligand");
        if(!ligand.from_sdf(sdf.c_str()))throw std::runtime_error("Cannot read ligand");
        Atom *n=nullptr,*h1=nullptr,*h2=nullptr,*c=nullptr,*o=nullptr;
        if(!lysine->identify_Schiff_amine(&n,&h1,&h2))throw std::runtime_error("Repaired lysine is not a Schiff-capable amine");
        if(!ligand.identify_Schiff_carbonyl(&c,&o))throw std::runtime_error("Fresh ligand lost its carbonyl");
        Molecule* water=ligand.create_Schiff_base(lysine);
        if(!water || n->is_bonded_to(c)!=2 || n->distance_to(c)>1.4)
            throw std::runtime_error("Fresh Schiff formation failed");
        delete water;
        std::cout<<"Fresh ligand formed a Schiff bond with the repaired lysine.\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<error.what()<<"\n";
        return 1;
    }
}
