#include "../classes/aminoacid.h"
#include <cassert>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <sys/stat.h>
#include <limits.h>

int main()
{
    char original[PATH_MAX];assert(getcwd(original,sizeof(original)));
    char work[]="/tmp/aroma-amino-cache-XXXXXX";assert(mkdtemp(work));
    std::string root=work;
    assert(!symlink((std::string(original)+"/data").c_str(),(root+"/data").c_str()));
    assert(!mkdir((root+"/sdf").c_str(),0700));
    assert(!chdir(work));
    FILE* cache=fopen("sdf/alanine.sdf","w");assert(cache);
    fputs("truncated cache\n",cache);fclose(cache);
    AminoAcid residue('A');
    assert(residue.get_CA() && residue.get_atom("N") && residue.get_atom("C"));
    std::ifstream input("sdf/alanine.sdf");
    std::ostringstream text;text<<input.rdbuf();input.close();
    assert(text.str().find("M  END")!=std::string::npos);
    Molecule restored;
    assert(restored.from_sdf(text.str().c_str())>0);
    const std::string letters="ACDEFGHIKLMNPQRSTVWY";
    for (char letter : letters)
    {
        printf("Constructing amino acid %c from an isolated cache.\n",letter);fflush(stdout);
        AminoAcid amino(letter,nullptr,false);
        assert(amino.get_CA() && amino.get_atom("N") && amino.get_atom("C"));
    }
    assert(!chdir(original));
    for (char letter : letters) unlink((root+"/sdf/"+aa_defs[int(letter)].name+".sdf").c_str());
    rmdir((root+"/sdf").c_str());unlink((root+"/data").c_str());rmdir(work);
    puts("Corrupt amino-acid caches rebuild; all 20 standard amino acids construct successfully.");
}
