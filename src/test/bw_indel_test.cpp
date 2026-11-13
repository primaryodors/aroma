#include "../classes/protein.h"
#include <cassert>
#include <iostream>
#include <sstream>
#include <memory>

int main(int argc,char** argv)
{
    if(argc==2 && std::string(argv[1])=="--probe")
    {
        std::unique_ptr<Protein> p(new Protein());
        std::string line;
        while(std::getline(std::cin,line))
        {
            std::istringstream input(line);
            char command;int helix=0,position=0,value=0;
            input>>command>>helix>>position>>value;
            if(command=='R')p.reset(new Protein());
            else if(command=='A')p->set_bw50(helix,position);
            else if(command=='I'||command=='D')p->set_bw_indel(helix,position,value,command=='I');
            else if(command=='Q')std::cout<<p->resolve_bw_resno(helix,position,value)<<'\n';
        }
        return 0;
    }
    Protein p;
    assert(!p.has_bwnos());
    p.set_bw50(5,215);
    assert(p.has_bwnos());
    std::string sequence(317,'A');sequence[207]='N';sequence[208]='I';sequence[214]='D';
    for(size_t i=0;i<sequence.size();++i)assert(p.add_residue(i+1,sequence[i]));
    p.set_bw50(1,48);p.set_bw50(5,215);p.set_bw50(7,294);
    p.set_region("TMR5",198,234);p.set_region("TMR7",274,299);
    p.set_bw_indel(5,49,1,false);
    assert(p.get_residue_bw("N5.42")->get_residue_no()==208);
    assert(p.get_residue_bw("STNQ5.43")==nullptr);
    assert(p.get_residue_bw(5,43)->get_residue_no()==209);
    assert(p.get_residue_bw(5,49)==nullptr);
    p.set_bw_indel(5,51,1,true);p.set_bw_indel(5,56,1,false);
    assert(p.resolve_bw_resno(5,51)==216);
    assert(p.resolve_bw_resno(5,51,1)==217);
    assert(p.resolve_bw_resno(5,52)==218);
    assert(p.resolve_bw_resno(5,56)==0);
    assert(p.resolve_bw_resno(5,57)==222);
    auto reverse=p.get_bw_from_resno(208);
    assert(reverse.helix_no==5 && reverse.member_no==42);
    ResiduePlaceholder placeholder;placeholder.set("N5.42");placeholder.resolve_resno(&p);
    assert(placeholder.resno==208);
    placeholder.set("5.43");placeholder.resolve_resno(&p);assert(placeholder.resno==209);
    placeholder.set("5.s");placeholder.resolve_resno(&p);assert(placeholder.resno==198);
    placeholder.set("5.e");placeholder.resolve_resno(&p);assert(placeholder.resno==234);
    p.set_bw_indel(7,36,1,false);p.set_bw_indel(7,36,2,true);
    assert(p.resolve_bw_resno(7,36)==0);
    assert(p.resolve_bw_resno(7,36,1)==279);
    assert(p.resolve_bw_resno(7,36,2)==280);
    assert(p.resolve_bw_resno(7,37)==281);
    assert(p.resolve_bw_resno(7,36,3)==0);
    reverse=p.get_bw_from_resno(280);
    assert(reverse.to_string()=="7.36+2");
    assert(p.get_residue(BallesterosWeinstein("7.36+2"))->get_residue_no()==280);
    FILE* file=tmpfile();assert(file);p.save_pdb(file);p.end_pdb(file);rewind(file);
    Protein restored;restored.load_pdb(file);fclose(file);
    assert(restored.resolve_bw_resno(5,42)==208);
    assert(restored.resolve_bw_resno(7,36,2)==280);
    Protein copied;copied.set_bw50(5,215);copied.copy_bw_indels_from(&restored);
    assert(copied.resolve_bw_resno(5,42)==208);
    copied.copy_bw_indels_from(&restored);
    assert(copied.resolve_bw_resno(5,42)==208);
    copied.clear_bw_indels();
    assert(copied.resolve_bw_resno(5,42)==207);
    copied.copy_bw_indels_from(&restored);
    assert(copied.resolve_bw_resno(5,42)==208);
    file=tmpfile();assert(file);
    fputs("REMARK 800 SITE BW 5.50 215\nEND\n",file);rewind(file);
    restored.load_pdb(file);fclose(file);
    assert(restored.resolve_bw_resno(5,42)==207);
    std::cout<<"BW indel lookup, reverse lookup, placeholders, persistence, and reload checks pass.\n";
}
