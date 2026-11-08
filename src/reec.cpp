// REEC: Restrained Ensemble Energy Conformer.
// Geometry comes from the project's amino acid library; optimization is local.
#include "classes/protein.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <vector>
#include <cstdio>
#include <unistd.h>
namespace reec
{
struct V
{
    double x=0,y=0,z=0;
    V() = default;
    V(double a,double b,double c):x(a),y(b),z(c)
    {
    }
    V(const Point& p):x(p.x),y(p.y),z(p.z)
    {
    }
    V operator+(V b) const
    {
        return
        {
            x+b.x,y+b.y,z+b.z
        };
    }
    V operator-(V b) const
    {
        return
        {
            x-b.x,y-b.y,z-b.z
        };
    }
    V operator*(double b) const
    {
        return
        {
            x*b,y*b,z*b
        };
    }
    V& operator+=(V b)
    {
        x+=b.x;
        y+=b.y;
        z+=b.z;
        return *this;
    }
    double dot(V b) const
    {
        return x*b.x+y*b.y+z*b.z;
    }
    V cross(V b) const
    {
        return
        {
            y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x
        };
    }
    double norm() const
    {
        return std::sqrt(dot(*this));
    }
    V unit() const
    {
        double d=norm();
        if(d<1e-8) throw std::runtime_error("Degenerate backbone frame");
        return *this*(1/d);
    }
};
using XYZ=std::vector<V>;
double dot(const XYZ& a,const XYZ& b)
{
    double r=0;
    for(size_t i=0; i<a.size(); ++i) r+=a[i].dot(b[i]);
    return r;
}
std::string trim(std::string s)
{
    auto a=s.find_first_not_of(" \r\n\t");
    return a==s.npos?"":s.substr(a,s.find_last_not_of(" \r\n\t")-a+1);
}
std::vector<std::string> split(std::string s,char c)
{
    std::vector<std::string> v;
    std::stringstream f(s);
    while(std::getline(f,s,c))v.push_back(trim(s));
    return v;
}
int integer(const std::string& s)
{
    size_t p=0;
    int n=std::stoi(s,&p);
    if(p!=s.size())throw std::runtime_error("Invalid integer: "+s);
    return n;
}
double number(const std::string& s)
{
    size_t p=0;
    double n=std::stod(s,&p);
    if(p!=s.size()||!std::isfinite(n))throw std::runtime_error("Invalid number: "+s);
    return n;
}
const std::string codes="ACDEFGHIKLMNPQRSTVWY";
const std::vector<std::string> names=
{
    "ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL",
    "TRP","TYR"
};
char letter(std::string s)
{
    for(size_t i=0; i<names.size(); ++i)if(s==names[i])return codes[i];
    if(s=="MSE")return 'M';
    return 0;
}
struct Entry
{
    std::string id,desc,seq;
    bool structure=false;
};
std::vector<Entry> alignment(const std::string& path)
{
    std::ifstream f(path);
    if(!f)throw std::runtime_error("Cannot read alignment: "+path);
    std::vector<Entry> a;
    std::string s;
    bool desc=false,done=true;
    while(std::getline(f,s))
    {
        s=trim(s);
        if(s.empty())continue;
        if(s[0]=='>')
        {
            if(!done)throw std::runtime_error("Unterminated PIR entry");
            auto p=s.find(';');
            if(p==s.npos)throw std::runtime_error("Expected PIR header >P1;ID");
            a.push_back(
            {
                s.substr(p+1),"","",false
            });
            desc=true;
            done=false;
            continue;
        }
        if(done&&!a.empty()&&s.find_first_not_of("-* \t")==s.npos)continue;
        if(a.empty()||done)throw std::runtime_error("Unexpected PIR content");
        if(desc)
        {
            a.back().desc=s;
            a.back().structure=s.compare(0,9,"structure")==0;
            if(!a.back().structure&&s.compare(0,8,"sequence")!=0)throw std::runtime_error("Expected sequence or structure descriptor");
            desc=false;
            continue;
        }
        for(char c:s)
        {
            if(std::isspace(static_cast<unsigned char>(c)))continue;
            if(c=='*')
            {
                done=true;
                break;
            }
            if(c!='-'&&c!='/'&&c!='.'&&codes.find(c)==codes.npos)throw std::runtime_error("Unsupported residue in PIR");
            a.back().seq+=c;
        }
    }
    if(!done||a.empty())throw std::runtime_error("Incomplete PIR alignment");
    std::set<std::string> ids;
    for(auto& e:a)if(e.seq.empty()||!ids.insert(e.id).second)throw std::runtime_error("Empty or duplicate PIR entry");
    return a;
}
void protein_row(const Entry& entry)
{
    if(entry.seq.find_first_of("/.")!=std::string::npos)
        throw std::runtime_error("Only a single protein chain is supported in PIR; omit ligand columns: "+entry.id);
}
struct Res
{
    int no;
    char aa;
    std::map<std::string,V> atoms;
};
struct Template
{
    std::string path,id;
    std::vector<Res> residues;
    std::set<std::pair<int,int>> links;
    std::vector<int> mapping;
};
Template pdb(const std::string& arg,char chain)
{
    Template t;
    auto eq=arg.find('=');
    t.path=eq==arg.npos?arg:arg.substr(eq+1);
    t.id=eq==arg.npos?"":arg.substr(0,eq);
    if(eq==arg.npos&&arg.size()==4&&std::all_of(arg.begin(),arg.end(),[](unsigned char c){return std::isalnum(c);}))
    {
        std::ifstream local(arg);
        if(!local)
        {
            t.id=arg;
            t.path="hm/tpl/"+arg+".pdb";
        }
    }
    std::ifstream f(t.path);
    if(!f)throw std::runtime_error("Cannot read template: "+t.path);
    std::string s;
    std::map<int,int> serial;
    std::vector<std::pair<int,int>> conect;
    std::map<int,int> index;
    bool ended=false;
    while(std::getline(f,s))
    {
        if(s.compare(0,6,"ENDMDL")==0)
        {
            ended=true;
            continue;
        }
        if(!ended&&s.size()>=54&&(s.compare(0,6,"ATOM  ")==0||s.compare(0,6,"HETATM")==0)&&s[21]==chain)
        {
            if(s[16]!=' '&&s[16]!='A')continue;
            char aa=letter(s.substr(17,3));
            if(!aa)continue;
            if(s[26]!=' ')throw std::runtime_error("Insertion codes require a renumbered template: "+t.path);
            int n=integer(trim(s.substr(22,4)));
            if(!index.count(n))
            {
                index[n]=t.residues.size();
                t.residues.push_back(
                {
                    n,aa, {
                    }
                });
            }
            Res& r=t.residues[index[n]];
            std::string atom=trim(s.substr(12,4));
            if(s.substr(17,3)=="MSE"&&atom=="SE")atom="SD";
            V p
            {
                number(trim(s.substr(30,8))),number(trim(s.substr(38,8))),number(trim(s.substr(46,8)))
            };
            if(!r.atoms.count(atom)||s[16]==' ')r.atoms[atom]=p;
            if(atom=="SG"&&aa=='C')serial[integer(trim(s.substr(6,5)))]=n;
        }
        if(s.compare(0,6,"SSBOND")==0&&s.size()>=35&&s[15]==chain&&s[29]==chain)
        {
            int a=integer(trim(s.substr(17,4))),b=integer(trim(s.substr(31,4)));
            t.links.insert(std::minmax(a,b));
        }
        if(s.compare(0,6,"CONECT")==0&&s.size()>=11)
        {
            int a=integer(trim(s.substr(6,5)));
            for(size_t i=11; i<s.size(); i+=5)
            {
                auto v=trim(s.substr(i,5));
                if(!v.empty())conect.push_back(
                {
                    a,integer(v)
                });
            }
        }
    }
    for(auto p:conect)if(serial.count(p.first)&&serial.count(p.second)&&p.first!=p.second)t.links.insert(std::minmax(serial[p.first],
                    serial[p.second]));
    if(t.residues.empty())throw std::runtime_error("No protein atoms in requested template chain");
    for(auto& r:t.residues)for(auto a:
                                   {
                                       "N","CA","C"
                                   })if(!r.atoms.count(a))throw std::runtime_error("Template backbone incomplete at "+std::to_string(r.no));
    return t;
}
// Global sequence alignment for a target-only PIR file. Explicit PIR columns take priority.
std::vector<int> align(const std::string& target,const Template& t)
{
    int n=target.size(),m=t.residues.size();
    std::vector<int> d((n+1)*(m+1));
    std::vector<char> trace(d.size());
    for(int i=1; i<=n; ++i)
    {
        d[i*(m+1)]=-3*i;
        trace[i*(m+1)]=1;
    }
    for(int j=1; j<=m; ++j)
    {
        d[j]=-3*j;
        trace[j]=2;
    }
    for(int i=1; i<=n; ++i)for(int j=1; j<=m; ++j)
        {
            int k=i*(m+1)+j;
            int v=d[k-m-2]+(target[i-1]==t.residues[j-1].aa?4:-1);
            char c=0;
            if(d[k-m-1]-3>v)
            {
                v=d[k-m-1]-3;
                c=1;
            }
            if(d[k-1]-3>v)
            {
                v=d[k-1]-3;
                c=2;
            }
            d[k]=v;
            trace[k]=c;
        }
    std::vector<int> map(n,-1);
    int i=n,j=m;
    while(i||j)
    {
        char c=trace[i*(m+1)+j];
        if(!c)
        {
            map[--i]=--j;
        }
        else if(c==1)--i;
        else --j;
    }
    return map;
}
void map_template(Template& t,const std::vector<Entry>& entries,const Entry& target,const std::string& seq,bool single_template=false)
{
    const Entry* e=nullptr;
    for(auto& a:entries)if(a.structure&&(a.id==t.id||t.path==a.id+".pdb"||t.path.substr(t.path.find_last_of('/')+1)==a.id+".pdb"))
        {
            if(e)throw std::runtime_error("Ambiguous template alignment");
            e=&a;
        }
    if(!e&&t.id.empty()&&single_template)
    {
        int count=0;
        for(auto& a:entries)if(a.structure)
            {
                e=&a;
                ++count;
            }
        if(count!=1)e=nullptr;
    }
    if(!e)
    {
        if(!t.id.empty())throw std::runtime_error("Template ID absent from PIR: "+t.id);
        t.mapping=align(seq,t);
        return;
    }
    protein_row(*e);
    auto fields=split(e->desc,':');
    int start=0;
    if(fields.size()>2&&fields[2]!="FIRST"&&!fields[2].empty())
    {
        int n=integer(fields[2]);
        while(start<int(t.residues.size())&&t.residues[start].no<n)++start;
    }
    int j=start;
    for(size_t c=0; c<std::max(target.seq.size(),e->seq.size()); ++c)
    {
        char template_aa=c<e->seq.size()?e->seq[c]:'-';
        char target_aa=c<target.seq.size()?target.seq[c]:'-';
        int k=-1;
        if(template_aa!='-')
        {
            if(j>=int(t.residues.size())||template_aa!=t.residues[j].aa)throw std::runtime_error("PIR/template sequence mismatch in "+t.path+" at column "+std::to_string(c+1));
            k=j++;
        }
        if(target_aa!='-')t.mapping.push_back(k);
    }
    if(j!=int(t.residues.size()))throw std::runtime_error("Template contains residues outside PIR row");
}
struct Frame
{
    V o,x,y,z;
    Frame(V n,V ca,V c):o(ca)
    {
        x=(c-ca).unit();
        z=x.cross(n-ca).unit();
        y=z.cross(x);
    }
    V local(V p)const
    {
        p=p-o;
        return
        {
            p.dot(x),p.dot(y),p.dot(z)
        };
    }
    V world(V p)const
    {
        return o+x*p.x+y*p.y+z*p.z;
    }
};
struct RigidFit
{
    V from,to,axis;
    double scalar=1;
    V apply(V p) const
    {
        p=p-from;
        return p+axis.cross(p)*(2*scalar)+axis.cross(axis.cross(p))*2+to;
    }
};
// Horn's quaternion fit, solved by Jacobi diagonalization, without an external library.
RigidFit superpose(Template& moving,const Template& reference)
{
    XYZ a,b;
    for(size_t i=0; i<moving.mapping.size(); ++i)if(moving.mapping[i]>=0&&reference.mapping[i]>=0)
        {
            a.push_back(moving.residues[moving.mapping[i]].atoms.at("CA"));
            b.push_back(reference.residues[reference.mapping[i]].atoms.at("CA"));
        }
    if(a.size()<3)throw std::runtime_error("Templates have fewer than three common CA anchors");
    V ca,cb;
    for(size_t i=0; i<a.size(); ++i)
    {
        ca+=a[i];
        cb+=b[i];
    }
    ca=ca*(1./a.size());
    cb=cb*(1./a.size());
    double s[3][3]=
    {
    };
    for(size_t i=0; i<a.size(); ++i)
    {
        V u=a[i]-ca,v=b[i]-cb;
        double x[]=
        {
            u.x,u.y,u.z
        },y[]=
        {
            v.x,v.y,v.z
        };
        for(int j=0; j<3; ++j)for(int k=0; k<3; ++k)s[j][k]+=x[j]*y[k];
    }
    double q[4][4]=
    {
        {
            s[0][0]+s[1][1]+s[2][2],s[1][2]-s[2][1],s[2][0]-s[0][2],s[0][1]-s[1][0]
        }, {
            s[1][2]-s[2][1],s[0][0]-s[1][1]-s[2][2],s[0][1]+s[1][0],s[0][2]+s[2][0]
        }, {
            s[2][0]-s[0][2],s[0][1]+s[1][0],-s[0][0]+s[1][1]-s[2][2],s[1][2]+s[2][1]
        }, {
            s[0][1]-s[1][0],s[0][2]+s[2][0],s[1][2]+s[2][1],-s[0][0]-s[1][1]+s[2][2]
        }
    };
    double ev[4][4]=
    {
    };
    for(int i=0; i<4; ++i)ev[i][i]=1;
    for(int it=0; it<60; ++it)
    {
        int p=0,r=1;
        for(int j=0; j<4; ++j)for(int k=j+1; k<4; ++k)if(std::fabs(q[j][k])>std::fabs(q[p][r]))
                {
                    p=j;
                    r=k;
                }
        if(std::fabs(q[p][r])<1e-12)break;
        double angle=.5*std::atan2(2*q[p][r],q[r][r]-q[p][p]),c=std::cos(angle),sine=std::sin(angle);
        double pp=q[p][p],rr=q[r][r],pr=q[p][r];
        q[p][p]=c*c*pp-2*c*sine*pr+sine*sine*rr;
        q[r][r]=sine*sine*pp+2*c*sine*pr+c*c*rr;
        q[p][r]=q[r][p]=0;
        for(int k=0; k<4; ++k)
        {
            if(k!=p&&k!=r)
            {
                double kp=q[k][p],kr=q[k][r];
                q[k][p]=q[p][k]=c*kp-sine*kr;
                q[k][r]=q[r][k]=sine*kp+c*kr;
            }
            double kp=ev[k][p],kr=ev[k][r];
            ev[k][p]=c*kp-sine*kr;
            ev[k][r]=sine*kp+c*kr;
        }
    }
    int best=0;
    for(int i=1; i<4; ++i)if(q[i][i]>q[best][best])best=i;
    double w=ev[0][best];
    V v
    {
        ev[1][best],ev[2][best],ev[3][best]
    };
    RigidFit fit{ca,cb,v,w};
    for(auto& r:moving.residues)for(auto& atom:r.atoms) atom.second=fit.apply(atom.second);
    return fit;
}
struct Node
{
    Atom* atom;
    int res;
    std::string name;
    double radius;
    bool backbone;
};
struct Pair
{
    int a,b;
    double lo,hi,weight;
    bool bond=false;
    bool validate=true;
};
struct Volume
{
    int a,b,c,d;
    double value,weight;
};
struct Contact
{
    int r1,r2;
    std::string a1,a2,kind;
    double lo,hi,weight;
    bool hard;
    std::vector<std::pair<int,int>> alternatives;
};
struct Stats
{
    double score=0,bond_max=0,angle_max=0,overlap_max=0,ca_rmsd=0,planarity_max=0;
    int severe=0,chirality=0,hard_fail=0;
};
struct Model
{
    struct Rotation
    {
        int root,tip;
        std::vector<int> moving;
    };
    Protein protein;
    std::vector<Node> nodes;
    XYZ xyz,ideal,anchors;
    std::vector<double> anchor_weight;
    std::vector<Pair> geometry;
    std::vector<Volume> volumes;
    std::vector<Volume> torsions;
    std::vector<Volume> dihedral_restraints;
    std::vector<std::pair<int,int>> mobile_ranges;
    std::vector<double> original_anchor_weight;
    std::vector<Contact> contacts;
    std::map<std::pair<int,std::string>,int> ids;
    std::vector<std::vector<int>> adj;
    std::set<std::pair<int,int>> excluded;
    std::vector<std::pair<int,int>> neighbors;
    std::vector<std::pair<int,int>> disulfides;
    XYZ neighbor_origin;
    std::vector<std::unique_ptr<Atom>> external;
    std::map<int,std::string> fixed_records;
    std::vector<std::pair<int,int>> fixed_links;
    double hard_scale=1,anchor_scale=1;
    unsigned seed=0;
    double clash_scale=1;
    int candidate_index=0;
    std::string start_signature;
    std::set<int> released_template_residues;
    std::vector<double> template_anchor_weights,template_geometry_weights;
    std::set<int> anneal_residues;
    struct Annealing
    {
        int trials=0,accepted=0,uphill=0,incremental=0,bonds=0;
        double before=0,after=0;
    } sampling;
    int id(int r,const std::string& a)const
    {
        auto it=ids.find(
        {
            r,a
        });
        if(it==ids.end())throw std::runtime_error("Unknown atom "+std::to_string(r)+":"+a);
        return it->second;
    }
    void pair(int a,int b,double target,double weight,bool bond=false)
    {
        geometry.push_back(
        {
            a,b,target,target,weight,bond
        });
    }
    void pair(int r,const char* a,int s,const char* b,double distance,double w,bool bond=false)
    {
        pair(id(r,a),id(s,b),distance,w,bond);
    }
    void edge(int a,int b)
    {
        adj[a].push_back(b);
        adj[b].push_back(a);
    }
    void build(const std::string& seq,std::vector<Template>& templates)
    {
        // Use Protein/Atom topology directly; no separate residue dictionary.
        for(size_t r=0; r<seq.size(); ++r)
        {
            if(!protein.add_residue(r+1,seq[r]))throw std::runtime_error("Cannot create residue");
            AminoAcid* aa=protein.get_residue(r+1);
            // The shared residue builder omits the terminal carboxyl oxygen.
            if(r+1==seq.size() && !aa->get_atom("OXT"))
            {
                Atom* carbon=aa->get_atom("C");
                V center(carbon->loc);
                V u=(V(aa->get_atom("CA")->loc)-center).unit();
                V v=(V(aa->get_atom("O")->loc)-center).unit();
                V terminal=center+(u+v).unit()*(-1.25);
                Point point(terminal.x,terminal.y,terminal.z);
                Atom* oxygen=aa->add_atom("O","OXT",&point,carbon,1,-1);
                if(!oxygen) throw std::runtime_error("Cannot create terminal oxygen");
                oxygen->is_backbone=true;
            }
            for(int j=0; j<aa->get_atom_count(); ++j)
            {
                Atom* a=aa->get_atom(j);
                if(a->Z<=1||(std::string(a->name)=="OXT"&&r+1<seq.size()))continue;
                int i=nodes.size();
                ids[
                {
                    int(r+1),a->name
                }]=i;
                nodes.push_back(
                {
                    a,int(r+1),a->name,a->vdW_radius,a->is_backbone
                });
                ideal.emplace_back(a->loc);
            }
        }
        // Normalize peptide references while retaining the library's sidechain topology.
        // Some cached residue conformers have nonstandard carbonyl lengths or handedness.
        const double radians=std::acos(-1.)/180;
        for(size_t r=1; r<=seq.size(); ++r)
        {
            int ni=id(r,"N"),ca=id(r,"CA"),ci=id(r,"C"),oi=id(r,"O");
            Frame frame(ideal[ni],ideal[ca],ideal[ci]);
            V old_oxygen=ideal[oi]-ideal[ci];
            V transverse=old_oxygen-frame.x*old_oxygen.dot(frame.x);
            if(transverse.norm()<1e-6) transverse=frame.y;
            ideal[ni]=ideal[ca]+frame.x*(1.458*std::cos(111.2*radians))+frame.y*(1.458*std::sin(111.2*radians));
            ideal[ci]=ideal[ca]+frame.x*1.525;
            ideal[oi]=ideal[ci]+frame.x*(1.231*std::cos(59.2*radians))+transverse.unit()*(1.231*std::sin(59.2*radians));
            if(ids.count({r,"OXT"}))
            {
                V u=(ideal[ca]-ideal[ci]).unit(),v=(ideal[oi]-ideal[ci]).unit();
                ideal[id(r,"OXT")]=ideal[ci]+(u+v).unit()*(-1.25);
            }
        }
        xyz=ideal;
        adj.resize(nodes.size());
        for(size_t i=0; i<nodes.size(); ++i)for(size_t j=i+1; j<nodes.size(); ++j)if(nodes[i].res==nodes[j].res&&nodes[i].atom->is_bonded_to(nodes[j].atom))
                {
                    edge(i,j);
                    double d=(ideal[i]-ideal[j]).norm();
                    if(d<.9||d>2.2)throw std::runtime_error("Invalid amino acid library geometry");
                    pair(i,j,d,3000,true);
                }
        // Distances over two bonds preserve bond angles; rigid central bonds preserve planarity.
        for(size_t b=0; b<adj.size(); ++b)for(size_t j=0; j<adj[b].size(); ++j)for(size_t k=j+1; k<adj[b].size(); ++k)
                {
                    int a=adj[b][j],c=adj[b][k];
                    pair(a,c,(ideal[a]-ideal[c]).norm(),1200);
                }
        std::set<std::pair<int,int>> rigid;
        for(size_t b=0; b<adj.size(); ++b)for(int c:adj[b])if(int(b)<c)
                {
                    Bond* bond=nodes[b].atom->get_bond_between(nodes[c].atom);
                    // Library movement locks are not chemical rigidity: phi and psi remain free.
                    const std::string& first=nodes[b].name;
                    const std::string& second=nodes[c].name;
                    bool backbone_rotation=(first=="CA"&&(second=="N"||second=="C"))
                        ||(second=="CA"&&(first=="N"||first=="C"));
                    if(backbone_rotation)continue;
                    if(bond&&!bond->can_rotate)for(int a:adj[b])for(int d:adj[c])if(a!=c&&d!=int(b)&&a!=d)
                                {
                                    auto p=std::minmax(a,d);
                                    if(rigid.insert(p).second)pair(a,d,(ideal[a]-ideal[d]).norm(),600);
                                }
                }
        for(size_t b=0; b<adj.size(); ++b) for(int c:adj[b]) if(int(b)<c)
                {
                    Bond* bond=nodes[b].atom->get_bond_between(nodes[c].atom);
                    if(!bond || !bond->can_rotate || nodes[b].atom->Z!=6 || nodes[c].atom->Z!=6) continue;
                    if(nodes[b].atom->get_count_pi_bonds() || nodes[c].atom->get_count_pi_bonds()) continue;
                    if(nodes[b].name=="C" || nodes[c].name=="C") continue;
                    if(protein.get_residue(nodes[b].res)->get_letter()=='P') continue;
                    for(int a:adj[b]) for(int d:adj[c]) if(a!=c && d!=int(b) && a!=d)
                                torsions.push_back({a,int(b),c,d,0,.3});
                }
        for(size_t a=0; a<adj.size(); ++a)if(adj[a].size()>=3)
            {
                if(nodes[a].name=="CA") continue;
                if(nodes[a].name=="CB" && (seq[nodes[a].res-1]=='I' || seq[nodes[a].res-1]=='T')) continue;
                int b=adj[a][0],c=adj[a][1],d=adj[a][2];
                double v=(ideal[b]-ideal[a]).dot((ideal[c]-ideal[a]).cross(ideal[d]-ideal[a]));
                volumes.push_back(
                {
                    int(a),b,c,d,v,100
                });
            }
        // Named signed volumes enforce standard L alpha centers and I/T beta centers.
        for(size_t r=1; r<=seq.size(); ++r)
        {
            if(ids.count({r,"CB"}))
            {
                int ca=id(r,"CA"),nn=id(r,"N"),cc=id(r,"C"),cb=id(r,"CB");
                double v=(ideal[nn]-ideal[ca]).dot((ideal[cc]-ideal[ca]).cross(ideal[cb]-ideal[ca]));
                volumes.push_back({ca,nn,cc,cb,std::fabs(v),100});
            }
            if(seq[r-1]=='I' || seq[r-1]=='T')
            {
                int cb=id(r,"CB"),ca=id(r,"CA"),branch=id(r,seq[r-1]=='I'?"CG1":"OG1"),other=id(r,"CG2");
                double v=(ideal[ca]-ideal[cb]).dot((ideal[branch]-ideal[cb]).cross(ideal[other]-ideal[cb]));
                volumes.push_back({cb,ca,branch,other,std::fabs(v),100});
            }
        }
        const int n=seq.size();
        // Standard trans peptide geometry: 1.329 A C-N, 116.2/121.7 degree angles.
        const double pi=std::acos(-1.),a=116.2*pi/180,b=121.7*pi/180,o=123.0*pi/180;
        V C
        {
            0,0,0
        },N
        {
            1.329,0,0
        },CA
        {
            1.525*std::cos(a),1.525*std::sin(a),0
        };
        V next=N+V
        {
            -1.458*std::cos(b),-1.458*std::sin(b),0
            },O
        {
            1.231*std::cos(o),-1.231*std::sin(o),0
        };
        for(int r=1; r<n; ++r)
        {
            int c=id(r,"C"),nn=id(r+1,"N");
            edge(c,nn);
            pair(c,nn,1.329,3000,true);
            pair(r,"CA",r+1,"N",(CA-N).norm(),1200);
            pair(r,"O",r+1,"N",(O-N).norm(),1200);
            pair(r,"C",r+1,"CA",next.norm(),1200);
            pair(r,"CA",r+1,"CA",(CA-next).norm(),800);
            pair(r,"O",r+1,"CA",(O-next).norm(),800);
            volumes.push_back(
            {
                c,id(r,"CA"),id(r,"O"),nn,0,100
            });
            volumes.push_back({c,id(r,"CA"),nn,id(r+1,"CA"),0,100});
        }
        for(size_t aidx=0; aidx<adj.size(); ++aidx)
        {
            excluded.insert(
            {
                aidx,aidx
            });
            for(int bidx:adj[aidx])
            {
                excluded.insert(std::minmax(int(aidx),bidx));
                for(int cidx:adj[bidx])
                {
                    excluded.insert(std::minmax(int(aidx),cidx));
                    for(int didx:adj[cidx])excluded.insert(std::minmax(int(aidx),didx));
                }
            }
        }
        std::vector<int> source(n,-1),mapped(n,-1);
        for(int r=0; r<n; ++r)
        {
            double best=-1;
            for(size_t t=0; t<templates.size(); ++t)
            {
                int k=templates[t].mapping[r];
                if(k<0)continue;
                double quality=templates[t].residues[k].aa==seq[r]?2:1;
                if(quality>best)
                {
                    best=quality;
                    source[r]=t;
                    mapped[r]=k;
                }
            }
        }
        int count=std::count_if(source.begin(),source.end(),[](int s)
        {
            return s>=0;
        });
        if(!count)throw std::runtime_error("No aligned template residues");
        std::vector<V> cas(n);
        for(int r=0; r<n; ++r)if(source[r]>=0)cas[r]=templates[source[r]].residues[mapped[r]].atoms.at("CA");
        // Place untemplated loops along an arched bridge. Optimization closes all peptide bonds.
        for(int first=0; first<n;)
        {
            if(source[first]>=0)
            {
                ++first;
                continue;
            }
            int last=first;
            while(last<n&&source[last]<0)++last;
            V left,right,normal;
            if(first>0&&last<n)
            {
                left=cas[first-1];
                right=cas[last];
                V tangent=(right-left).unit();
                normal=tangent.cross(
                {
                    0,0,1
                });
                if(normal.norm()<.1)normal=tangent.cross(
                {
                    0,1,0
                });
                normal=normal.unit();
                double span=(right-left).norm(),length=3.8*(last-first+1),height=std::min(8.,std::sqrt(std::max(0.,length*length-span*span))*.3);
                for(int r=first; r<last; ++r)
                {
                    double u=double(r-first+1)/(last-first+1);
                    cas[r]=left*(1-u)+right*u+normal*(height*std::sin(pi*u));
                }
            }
            else
            {
                int anchor=first>0?first-1:last;
                V dir
                {
                    1,0,0
                };
                if(first>1)dir=(cas[first-1]-cas[first-2]).unit();
                else if(last+1<n)dir=(cas[last]-cas[last+1]).unit();
                for(int r=first; r<last; ++r)cas[r]=cas[anchor]+dir*(3.8*std::abs(r-anchor));
            }
            first=last;
        }
        anchor_weight.resize(nodes.size());
        for(int r=0; r<n; ++r)
        {
            Frame original(ideal[id(r+1,"N")],ideal[id(r+1,"CA")],ideal[id(r+1,"C")]);
            V nn,cc;
            if(source[r]>=0)
            {
                auto& rr=templates[source[r]].residues[mapped[r]];
                nn=rr.atoms.at("N");
                cc=rr.atoms.at("C");
            }
            else
            {
                V tangent=(r+1<n?cas[r+1]-cas[r]:cas[r]-cas[r-1]).unit();
                V norm=tangent.cross(
                {
                    0,0,1
                });
                if(norm.norm()<.1)norm=tangent.cross(
                {
                    0,1,0
                });
                norm=norm.unit();
                nn=cas[r]-tangent*.9+norm*1.15;
                cc=cas[r]+tangent*1.3+norm*.8;
            }
            Frame frame(nn,cas[r],cc);
            for(size_t i=0; i<nodes.size(); ++i)if(nodes[i].res==r+1)
                {
                    xyz[i]=frame.world(original.local(ideal[i]));
                    if(source[r]>=0)
                    {
                        auto& rr=templates[source[r]].residues[mapped[r]];
                        auto p=rr.atoms.find(nodes[i].name);
                        if(p!=rr.atoms.end()&&(nodes[i].backbone||rr.aa==seq[r]))xyz[i]=p->second;
                        anchor_weight[i]=nodes[i].backbone?1.5:.03;
                    }
                }
        }
        // Repair initial sidechain handedness without distorting covalent distances.
        for(int r=1; r<=n; ++r)
        {
            if(!ids.count({r,"CB"})) continue;
            int ca=id(r,"CA"),nn=id(r,"N"),cc=id(r,"C"),cb=id(r,"CB");
            V normal=(xyz[nn]-xyz[ca]).cross(xyz[cc]-xyz[ca]).unit();
            if((xyz[cb]-xyz[ca]).dot(normal)<0 && seq[r-1]!='P')
            {
                for(size_t i=0; i<nodes.size(); ++i) if(nodes[i].res==r && !nodes[i].backbone)
                    {
                        V v=xyz[i]-xyz[ca];
                        xyz[i]=xyz[i]-normal*(2*v.dot(normal));
                    }
            }
            if(seq[r-1]=='I' || seq[r-1]=='T')
            {
                int branch=id(r,seq[r-1]=='I'?"CG1":"OG1"),other=id(r,"CG2");
                normal=(xyz[ca]-xyz[cb]).cross(xyz[branch]-xyz[cb]).unit();
                if((xyz[other]-xyz[cb]).dot(normal)<0)
                {
                    for(size_t i=0; i<nodes.size(); ++i) if(nodes[i].res==r && !nodes[i].backbone && int(i)!=cb)
                        {
                            V v=xyz[i]-xyz[cb];
                            xyz[i]=xyz[i]-normal*(2*v.dot(normal));
                        }
                }
            }
        }
        anchors=xyz;
        // Each template contributes a soft, rotation-invariant local fold restraint.
        for(auto& t:templates) for(int r=0; r<n; ++r) if(t.mapping[r]>=0)
                    for(int q=r+3; q<n; ++q) if(t.mapping[q]>=0)
                        {
                            double d=(t.residues[t.mapping[r]].atoms.at("CA")-t.residues[t.mapping[q]].atoms.at("CA")).norm();
                            if(d<9) geometry.push_back({id(r+1,"CA"),id(q+1,"CA"),std::max(.5,d-.75),d+.75,.15/templates.size(),false,false});
                        }
        for(auto& t:templates)for(auto link:t.links)
            {
                int ares=0,bres=0;
                for(int r=0; r<n; ++r)if(t.mapping[r]>=0)
                    {
                        int no=t.residues[t.mapping[r]].no;
                        if(no==link.first)ares=r+1;
                        if(no==link.second)bres=r+1;
                    }
                if(ares&&bres&&seq[ares-1]=='C'&&seq[bres-1]=='C')add_disulfide(ares,bres);
            }
        std::cout<<"REEC built "<<n<<" residues, "<<nodes.size()<<" heavy atoms; template coverage "<<count<<"/"<<n
                 <<"\n";
    }
    void resolve(Contact& c)
    {
        auto aa=split(c.a1,'|'),bb=split(c.a2,'|');
        for(auto& a:aa)for(auto& b:bb)c.alternatives.push_back(
            {
                id(c.r1,a),id(c.r2,b)
            });
        for(auto p:c.alternatives)if(p.first==p.second)throw std::runtime_error("Self contact");
    }
    void add_disulfide(int a,int b)
    {
        std::pair<int,int> p=std::minmax(a,b);
        for(auto link:disulfides)if(link==p)return;
        for(auto link:disulfides)if(link.first==a||link.second==a||link.first==b||link.second==b)throw std::runtime_error("Conflicting disulfide topology");
        Contact c
        {
            a,b,"SG","SG","disulfide",1.95,2.15,100,true, {
            }
        };
        resolve(c);
        contacts.push_back(c);
        disulfides.push_back(p);
        int x=id(a,"SG"),y=id(b,"SG");
        edge(x,y);
        pair(x,y,2.05,3000,true);
        excluded.insert(std::minmax(x,y));
        for(int v:adj[x])excluded.insert(std::minmax(v,y));
        for(int v:adj[y])excluded.insert(std::minmax(v,x));
        // Approximately 104 degree C-S-S angles and either handed 90 degree torsion.
        pair(id(a,"CB"),y,3.05,500);
        pair(id(b,"CB"),x,3.05,500);
        pair(id(a,"CB"),id(b,"CB"),3.84,150);
    }
    void read_contacts(const std::string& path)
    {
        std::ifstream f(path);
        if(!f)throw std::runtime_error("Cannot read contact file: "+path);
        std::string s;
        int line=0;
        while(std::getline(f,s))
        {
            ++line;
            s=s.substr(0,s.find('#'));
            s=trim(s);
            if(s.empty())continue;
            try
            {
                std::stringstream ss(s);
                std::string kind;
                ss>>kind;
                if(kind=="disulfide")
                {
                    std::string a,b,extra;
                    if(!(ss>>a>>b)||(ss>>extra))throw std::runtime_error("Expected disulfide RES RES");
                    add_disulfide(integer(a),integer(b));
                    continue;
                }
                if(kind=="dihedral")
                {
                    std::string r[4], atom[4], mean, deviation, extra;
                    for(int i=0;i<4;++i)
                        if(!(ss>>r[i]>>atom[i]))throw std::runtime_error("Expected four dihedral atoms");
                    if(!(ss>>mean>>deviation)||(ss>>extra))throw std::runtime_error("Expected dihedral mean and deviation in radians");
                    double sd=number(deviation);
                    if(sd<=0)throw std::runtime_error("Invalid dihedral deviation");
                    dihedral_restraints.push_back({id(integer(r[0]),atom[0]),id(integer(r[1]),atom[1]),
                        id(integer(r[2]),atom[2]),id(integer(r[3]),atom[3]),number(mean),1/(sd*sd)});
                    continue;
                }
                if(kind=="mobile")
                {
                    std::string first,last,extra;
                    if(!(ss>>first>>last)||(ss>>extra))throw std::runtime_error("Expected mobile FIRST LAST");
                    int start=integer(first),end=integer(last);
                    if(start<1||end<start||end>protein.get_end_resno())throw std::runtime_error("Invalid mobile residue range");
                    mobile_ranges.push_back({start,end});
                    if(original_anchor_weight.empty())original_anchor_weight=anchor_weight;
                    for(size_t i=0;i<nodes.size();++i)
                    {
                        bool mobile=false;
                        for(auto range:mobile_ranges)
                            if(nodes[i].res>=range.first&&nodes[i].res<=range.second)mobile=true;
                        anchor_weight[i]=original_anchor_weight[i]*(mobile?.25:1);
                    }
                    continue;
                }
                Contact c;
                std::string mode,r1,r2,lo,hi,w,extra;
                if(!(ss>>mode>>r1>>c.a1>>r2>>c.a2>>lo>>hi>>w)||(ss>>extra))throw std::runtime_error("Expected contact hard|soft RES ATOM RES ATOM MIN MAX WEIGHT");
                if(kind!="contact"||(mode!="hard"&&mode!="soft"))throw std::runtime_error("Invalid contact type");
                c.kind="contact";
                c.hard=mode=="hard";
                c.r1=integer(r1);
                c.r2=integer(r2);
                c.lo=number(lo);
                c.hi=number(hi);
                c.weight=number(w);
                if(c.lo<=0||c.hi<c.lo||c.weight<=0)throw std::runtime_error("Invalid contact bounds or weight");
                resolve(c);
                contacts.push_back(c);
            }
            catch(const std::exception& e)
            {
                throw std::runtime_error(path+":"+std::to_string(line)+": "+e.what());
            }
        }
    }
    // Ligands remain rigid. Their coordinates still contribute to steric and contact terms.
    void environment(const std::string& path, const RigidFit& fit=RigidFit(), bool include_hydrogens=false)
    {
        std::ifstream f(path);
        if (!f) throw std::runtime_error("Cannot read fixed environment: " + path);
        std::string s;
        std::map<int,int> serial_to_index;
        std::vector<std::pair<int,int>> input_links;
        bool ended=false;
        while (std::getline(f,s))
        {
            if(s.compare(0,6,"ENDMDL")==0)ended=true;
            if(s.compare(0,6,"CONECT")==0&&s.size()>=11)
            {
                int a=integer(trim(s.substr(6,5)));
                for(size_t k=11; k<s.size(); k+=5)
                {
                    auto field=trim(s.substr(k,5));
                    if(!field.empty())input_links.push_back(
                    {
                        a,integer(field)
                    });
                }
            }
            if (ended || s.compare(0,6,"HETATM")!=0 || s.size()<54) continue;
            if (s[16]!=' ' && s[16]!='A') continue;
            if (letter(s.substr(17,3))) continue;
            int r=integer(trim(s.substr(22,4)));
            std::string name=trim(s.substr(12,4));
            if (ids.count(
        {
            r,name
        })) throw std::runtime_error("Fixed environment atom conflicts with target atom");
            std::string element=s.size()>=78?trim(s.substr(76,2)):"";
            if(element.empty())
            {
                for(char c:name) if(std::isalpha(static_cast<unsigned char>(c)))
                    {
                        element=std::string(1,c);
                        break;
                    }
            }
            if(element.empty())throw std::runtime_error("Unknown environment element");
            V p
            {
                number(trim(s.substr(30,8))),number(trim(s.substr(38,8))),number(trim(s.substr(46,8)))
            };
            p=fit.apply(p);
            Point point(p.x,p.y,p.z);
            std::unique_ptr<Atom> atom(new Atom(element.c_str(),&point));
            if(atom->Z<=0)throw std::runtime_error("Unknown environment element");
            if(atom->Z==1&&!include_hydrogens)continue;
            int i=nodes.size();
            ids[
            {
                r,name
            }]=i;
            nodes.push_back(
            {
                atom.get(),r,name,atom->vdW_radius,false
            });
            char coordinates[25];
            snprintf(coordinates,sizeof(coordinates),"%8.3f%8.3f%8.3f",p.x,p.y,p.z);
            s.replace(30,24,coordinates);
            fixed_records[i]=s;
            serial_to_index[integer(trim(s.substr(6,5)))]=i;
            external.push_back(std::move(atom));
            xyz.push_back(p);
            ideal.push_back(p);
            anchors.push_back(p);
            anchor_weight.push_back(0);
        }
        for(auto p:input_links)if(serial_to_index.count(p.first)&&serial_to_index.count(p.second))fixed_links.push_back(
            {
                serial_to_index[p.first],serial_to_index[p.second]
            });
        if (fixed_records.empty()) throw std::runtime_error("Fixed environment has no ligand heavy atoms");
    }
    void refresh(const XYZ& x)
    {
        neighbors.clear();
        std::map<std::array<int,3>,std::vector<int>> cells;
        double cell=6.;
        for(size_t i=0; i<x.size(); ++i)
        {
            std::array<int,3> key
            {
                {
                    int(std::floor(x[i].x/cell)),int(std::floor(x[i].y/cell)),int(std::floor(x[i].z/cell))
                }
            };
            for(int a=-1; a<=1; ++a)for(int b=-1; b<=1; ++b)for(int c=-1; c<=1; ++c)
                    {
                        auto it=cells.find(
                        {
                            {
                                key[0]+a,key[1]+b,key[2]+c
                            }
                        });
                        if(it!=cells.end())for(int j:it->second)if(!(fixed_records.count(j)&&fixed_records.count(i))&&!excluded.count(
                            {
                                j,int(i)
                                })&&(x[i]-x[j]).norm()<6.)neighbors.push_back(
                        {
                            j,int(i)
                        });
                    }
            cells[key].push_back(i);
        }
        neighbor_origin=x;
    }
    bool stale(const XYZ& x)const
    {
        if(neighbor_origin.size()!=x.size())return true;
        for(size_t i=0; i<x.size(); ++i)if((x[i]-neighbor_origin[i]).norm()>.6)return true;
        return false;
    }
    static double distance(const XYZ& x,XYZ& g,int a,int b,double lo,double hi,double w)
    {
        V v=x[a]-x[b];
        double d=std::max(1e-10,v.norm()),err=d<lo?d-lo:(d>hi?d-hi:0);
        V force=v*(2*w*err/d);
        if(d<1e-8)force=
        {
            2*w*err,0,0
        };
        g[a]+=force;
        g[b]+=force*(-1);
        return w*err*err;
    }
    // Threefold alkyl torsion potential; bonds and chiral volumes are independent terms.
    static double torsion_energy(const XYZ& x, XYZ& g, const Volume& t)
    {
        V u=x[t.b]-x[t.a], v=x[t.c]-x[t.b], w=x[t.d]-x[t.c];
        V n=u.cross(v), m=v.cross(w);
        double n2=n.dot(n), m2=m.dot(m);
        if(n2<1e-12 || m2<1e-12) return 0;
        double denominator=std::sqrt(n2*m2);
        double cosine=std::max(-1.,std::min(1.,n.dot(m)/denominator));
        double derivative=t.weight*(12*cosine*cosine-3);
        V gn=(m*(1/denominator)-n*(cosine/n2))*derivative;
        V gm=(n*(1/denominator)-m*(cosine/m2))*derivative;
        V gu=v.cross(gn), gv=gn.cross(u)+w.cross(gm), gw=gm.cross(v);
        g[t.a]+=gu*(-1);
        g[t.b]+=gu-gv;
        g[t.c]+=gv-gw;
        g[t.d]+=gw;
        return t.weight*(1+4*cosine*cosine*cosine-3*cosine);
    }
    static double dihedral_value(const XYZ& x,const Volume& t)
    {
        V u=x[t.b]-x[t.a],v=x[t.c]-x[t.b],w=x[t.d]-x[t.c];
        V n=u.cross(v),m=v.cross(w);
        return std::atan2(n.cross(m).dot(v)/std::max(1e-12,v.norm()),n.dot(m));
    }
    static double restrained_dihedral(const XYZ& x,XYZ& g,const Volume& t)
    {
        const int atoms[]={t.a,t.b,t.c,t.d};
        const Volume local{0,1,2,3,t.value,t.weight};
        auto energy=[&](const XYZ& points)
        {
            double error=std::remainder(dihedral_value(points,local)-t.value,2*std::acos(-1.));
            return t.weight*error*error;
        };
        XYZ shifted={x[t.a],x[t.b],x[t.c],x[t.d]};
        const double step=1e-6;
        for(int local_atom=0;local_atom<4;++local_atom)for(int axis=0;axis<3;++axis)
        {
            int atom=atoms[local_atom];
            double* position=axis==0?&shifted[local_atom].x:axis==1?&shifted[local_atom].y:&shifted[local_atom].z;
            *position+=step;
            double plus=energy(shifted);
            *position-=2*step;
            double minus=energy(shifted);
            *position+=step;
            double derivative=(plus-minus)/(2*step);
            if(axis==0)g[atom].x+=derivative;
            else if(axis==1)g[atom].y+=derivative;
            else g[atom].z+=derivative;
        }
        return energy(shifted);
    }
    void alternate_start(int candidate)
    {
        auto mixed=[](uint32_t key)
        {
            key^=key>>16;
            key*=0x7feb352du;
            key^=key>>15;
            key*=0x846ca68bu;
            return key^(key>>16);
        };
        if(seed||candidate)for(size_t i=0;i<xyz.size();++i)if(!fixed_records.count(i))
        {
            uint32_t key=mixed(uint32_t(seed)^uint32_t(candidate)*0x9e3779b9u^uint32_t(i)*0x85ebca6bu);
            xyz[i].x+=.02*(double(key&65535)/65535.-.5);
            key=mixed(key);
            xyz[i].y+=.02*(double(key&65535)/65535.-.5);
            key=mixed(key);
            xyz[i].z+=.02*(double(key&65535)/65535.-.5);
        }
        refresh(xyz);
        std::set<int> flexible;
        for(auto p:neighbors)
        {
            double overlap=nodes[p.first].radius+nodes[p.second].radius-(xyz[p.first]-xyz[p.second]).norm();
            if(overlap>.8)
            {
                if(!fixed_records.count(p.first)) flexible.insert(nodes[p.first].res);
                if(!fixed_records.count(p.second)) flexible.insert(nodes[p.second].res);
            }
        }
        for(auto& c:contacts) if(c.hard)
            {
                flexible.insert(c.r1);
                flexible.insert(c.r2);
            }
        for(int r:flexible)
        {
            if(!protein.get_residue(r) || !ids.count({r,"CB"}) || !ids.count({r,"CA"})) continue;
            if(protein.get_residue(r)->get_letter()=='P') continue;
            V origin=xyz[id(r,"CA")], axis=(xyz[id(r,"CB")]-origin).unit();
            uint32_t key=mixed(uint32_t(seed)^uint32_t(candidate)*0x9e3779b9u^uint32_t(r)*0x85ebca6bu);
            double angle=2*std::acos(-1.)/3*(key%3)+.2*(double(key>>8)/16777215.-.5);
            double c=std::cos(angle),s=std::sin(angle);
            for(size_t i=0; i<nodes.size(); ++i) if(nodes[i].res==r&&!nodes[i].backbone&&!fixed_records.count(i))
                {
                    V p=xyz[i]-origin;
                    xyz[i]=origin+p*c+axis.cross(p)*s+axis*(axis.dot(p)*(1-c));
                }
        }
    }
    std::vector<Rotation> targeted_rotations()
    {
        std::set<int> selected=anneal_residues;
        refresh(xyz);
        for(auto p:neighbors)
        {
            double distance=(xyz[p.first]-xyz[p.second]).norm();
            bool ligand=(fixed_records.count(p.first)||fixed_records.count(p.second));
            if(ligand||nodes[p.first].radius+nodes[p.second].radius-distance>.5)
            {
                selected.insert(nodes[p.first].res);
                selected.insert(nodes[p.second].res);
            }
        }
        for(auto range:mobile_ranges)for(int r=range.first;r<=range.second;++r)selected.insert(r);
        std::vector<Rotation> result;
        for(size_t a=0;a<adj.size();++a)for(int b:adj[a])if(int(a)<b&&nodes[a].res==nodes[b].res
            &&selected.count(nodes[a].res))
        {
            Bond* bond=nodes[a].atom->get_bond_between(nodes[b].atom);
            if(!bond||!bond->can_rotate)continue;
            for(int orientation=0;orientation<2;++orientation)
            {
                int root=orientation?b:int(a),tip=orientation?int(a):b;
                if(nodes[tip].backbone)continue;
                std::set<int> visited{tip};
                std::vector<int> moving{tip};
                bool blocked=false;
                for(size_t i=0;i<moving.size()&&!blocked;++i)
                {
                    int atom=moving[i];
                    if(nodes[atom].backbone||fixed_records.count(atom)){blocked=true;break;}
                    for(int next:adj[atom])
                    {
                        if(atom==tip&&next==root)continue;
                        if(next==root){blocked=true;break;}
                        if(visited.insert(next).second)moving.push_back(next);
                    }
                }
                if(!blocked&&moving.size()>1)
                {
                    result.push_back({root,tip,moving});
                    break;
                }
            }
        }
        return result;
    }
    void rotate(const Rotation& bond,double angle)
    {
        V origin=xyz[bond.root],axis=(xyz[bond.tip]-origin).unit();
        double cosine=std::cos(angle),sine=std::sin(angle);
        for(int atom:bond.moving)
        {
            V point=xyz[atom]-origin;
            xyz[atom]=origin+point*cosine+axis.cross(point)*sine+axis*(axis.dot(point)*(1-cosine));
        }
    }
    void anneal(int steps,double search_heat,int candidate)
    {
        sampling.trials=sampling.accepted=sampling.uphill=sampling.incremental=sampling.bonds=0;
        XYZ gradient;
        refresh(xyz);
        double score=energy(xyz,gradient);
        sampling.before=sampling.after=score;
        if(!steps)return;
        auto rotations=targeted_rotations();
        sampling.bonds=rotations.size();
        if(rotations.empty())return;
        std::mt19937 random(seed^unsigned(candidate+1)*0x9e3779b9u);
        auto uniform=[&](){return (double(random())+.5)/4294967296.;};
        XYZ best=xyz;
        double best_score=score;
        const double pi=std::acos(-1.);
        auto evaluate=[&](){refresh(xyz);return energy(xyz,gradient);};
        for(int step=0;step<steps;++step)
        {
            double fraction=steps>1?double(step)/(steps-1):1.;
            double heat=search_heat*std::pow(.01,fraction);
            const Rotation& bond=rotations[random()%rotations.size()];
            XYZ prior=xyz;
            double angle=(2*uniform()-1)*pi;
            rotate(bond,angle);
            double trial=evaluate(),delta=trial-score;
            ++sampling.trials;
            if(delta<=0||uniform()<std::exp(-delta/heat))
            {
                ++sampling.accepted;
                if(delta>0)++sampling.uphill;
                score=trial;
            }
            else xyz=prior;
            if(score<best_score){best=xyz;best_score=score;}

            // Match docking's incremental search with both rotation directions.
            prior=xyz;
            XYZ local_best=prior;
            double local_score=score;
            double increment=(2+10*(1-fraction))*pi/180.;
            for(int direction:{-1,1})
            {
                xyz=prior;
                rotate(bond,direction*increment);
                trial=evaluate();
                if(trial<local_score){local_score=trial;local_best=xyz;}
            }
            if(local_score<score)++sampling.incremental;
            xyz=local_best;
            score=local_score;
            if(score<best_score){best=xyz;best_score=score;}
        }
        xyz=best;
        refresh(xyz);
        sampling.after=best_score;
        std::cout<<"REEC annealing: "<<sampling.bonds<<" torsions, "<<sampling.accepted<<"/"<<sampling.trials
                 <<" Monte Carlo moves accepted, "<<sampling.incremental<<" incremental improvements; score "
                 <<sampling.before<<" -> "<<sampling.after<<"\n";
    }
    double energy(const XYZ& x,XYZ& g)
    {
        g.assign(x.size(),
        {
        });
        double e=0;
        for(auto& p:geometry)e+=distance(x,g,p.a,p.b,p.lo,p.hi,p.weight);
        for(auto& v:volumes)
        {
            V b=x[v.b]-x[v.a],c=x[v.c]-x[v.a],d=x[v.d]-x[v.a];
            double err=b.dot(c.cross(d))-v.value,w=v.weight;
            if(std::fabs(v.value)<.1)w*=3;
            e+=w*err*err;
            V gb=c.cross(d)*(2*w*err),gc=d.cross(b)*(2*w*err),gd=b.cross(c)*(2*w*err);
            g[v.b]+=gb;
            g[v.c]+=gc;
            g[v.d]+=gd;
            g[v.a]+=(gb+gc+gd)*(-1);
        }
        for(auto& t:torsions) e+=torsion_energy(x,g,t);
        for(auto& t:dihedral_restraints)e+=restrained_dihedral(x,g,t);
        for(size_t i=0; i<x.size(); ++i)
        {
            V v=x[i]-anchors[i];
            double w=anchor_weight[i]*anchor_scale;
            e+=w*v.dot(v);
            g[i]+=v*(2*w);
        }
        for(auto& c:contacts)
        {
            double best=1e100;
            std::pair<int,int> chosen;
            for(auto p:c.alternatives)
            {
                double d=(x[p.first]-x[p.second]).norm(),err=d<c.lo?c.lo-d:(d>c.hi?d-c.hi:0);
                if(err<best)
                {
                    best=err;
                    chosen=p;
                }
            }
            e+=distance(x,g,chosen.first,chosen.second,c.lo,c.hi,c.weight*(c.hard?hard_scale:1));
        }
        for(auto p:neighbors)
        {
            double limit=nodes[p.first].radius+nodes[p.second].radius-.45;
            if((nodes[p.first].atom->Z==7||nodes[p.first].atom->Z==8)&&(nodes[p.second].atom->Z==7||nodes[p.second].atom->Z==8))limit=std::min(limit,
                        2.65);
            e+=distance(x,g,p.first,p.second,limit,1e10,30*clash_scale);
        }
        for(auto& fixed:fixed_records)g[fixed.first]=
        {
        };
        if(!std::isfinite(e))throw std::runtime_error("Nonfinite optimization energy");
        return e;
    }
    void minimize(int iterations)
    {
        refresh(xyz);
        XYZ g;
        double e=energy(xyz,g);
        std::vector<XYZ> sh,yh;
        std::vector<double> rho;
        for(int it=0; it<iterations; ++it)
        {
            if(stale(xyz))
            {
                refresh(xyz);
                e=energy(xyz,g);
                sh.clear();
                yh.clear();
                rho.clear();
            }
            XYZ direction=g;
            std::vector<double> alpha(sh.size());
            for(int h=int(sh.size())-1; h>=0; --h)
            {
                alpha[h]=rho[h]*dot(sh[h],direction);
                for(size_t i=0; i<g.size(); ++i)direction[i]+=yh[h][i]*(-alpha[h]);
            }
            double scale=sh.empty()?1./6000:dot(sh.back(),yh.back())/std::max(1e-20,dot(yh.back(),yh.back()));
            for(auto& v:direction)v=v*scale;
            for(size_t h=0; h<sh.size(); ++h)
            {
                double beta=rho[h]*dot(yh[h],direction);
                for(size_t i=0; i<g.size(); ++i)direction[i]+=sh[h][i]*(alpha[h]-beta);
            }
            for(auto& v:direction)v=v*(-1);
            double slope=dot(g,direction);
            if(slope>=0)
            {
                sh.clear();
                yh.clear();
                rho.clear();
                direction=g;
                for(auto& v:direction)v=v*(-1./6000);
                slope=dot(g,direction);
            }
            double maxstep=0;
            for(auto& v:direction)maxstep=std::max(maxstep,v.norm());
            if(maxstep<1e-9)break;
            double step=std::min(1.,.15/maxstep);
            XYZ trial,ng;
            double en=e;
            bool accepted=false;
            for(int ls=0; ls<25; ++ls)
            {
                trial=xyz;
                for(size_t i=0; i<g.size(); ++i)trial[i]+=direction[i]*step;
                en=energy(trial,ng);
                if(en<=e+1e-4*step*slope)
                {
                    accepted=true;
                    break;
                }
                step*=.5;
            }
            if(!accepted)break;
            XYZ s(g.size()),y(g.size());
            for(size_t i=0; i<g.size(); ++i)
            {
                s[i]=trial[i]-xyz[i];
                y[i]=ng[i]-g[i];
            }
            double sy=dot(s,y);
            if(sy>1e-10)
            {
                if(sh.size()==10)
                {
                    sh.erase(sh.begin());
                    yh.erase(yh.begin());
                    rho.erase(rho.begin());
                }
                sh.push_back(s);
                yh.push_back(y);
                rho.push_back(1/sy);
            }
            xyz.swap(trial);
            g.swap(ng);
            e=en;
            if(it%250==0)std::cout<<"  iteration "<<it<<" energy "<<e<<"\n";
        }
    }
    Stats stats()
    {
        refresh(xyz);
        XYZ g;
        Stats s;
        s.score=energy(xyz,g);
        int ca_count=0;
        for(size_t i=0; i<nodes.size(); ++i) if(nodes[i].name=="CA" && anchor_weight[i]>0)
            {
                V delta=xyz[i]-anchors[i];
                s.ca_rmsd+=delta.dot(delta);
                ++ca_count;
            }
        if(ca_count) s.ca_rmsd=std::sqrt(s.ca_rmsd/ca_count);
        for(auto p:geometry)
        {
            if(!p.validate) continue;
            double err=std::fabs((xyz[p.a]-xyz[p.b]).norm()-p.lo);
            if(p.bond)s.bond_max=std::max(s.bond_max,err);
            else s.angle_max=std::max(s.angle_max,err);
        }
        for(auto v:volumes) if(std::fabs(v.value)<.1)
            {
                double actual=(xyz[v.b]-xyz[v.a]).dot((xyz[v.c]-xyz[v.a]).cross(xyz[v.d]-xyz[v.a]));
                s.planarity_max=std::max(s.planarity_max,std::fabs(actual-v.value));
            }
        for(auto v:volumes)if(std::fabs(v.value)>.5)
            {
                double actual=(xyz[v.b]-xyz[v.a]).dot((xyz[v.c]-xyz[v.a]).cross(xyz[v.d]-xyz[v.a]));
                if(actual*v.value<=0)++s.chirality;
            }
        for(auto p:neighbors)
        {
            double overlap=nodes[p.first].radius+nodes[p.second].radius-(xyz[p.first]-xyz[p.second]).norm();
            s.overlap_max=std::max(s.overlap_max,overlap);
            if(overlap>1.0)++s.severe;
        }
        for(auto& c:contacts)
        {
            double err=1e100;
            for(auto p:c.alternatives)
            {
                double d=(xyz[p.first]-xyz[p.second]).norm();
                err=std::min(err,std::max(
                {
                    0.,c.lo-d,d-c.hi
                }));
            }
            if(c.hard&&err>.05)++s.hard_fail;
        }
        return s;
    }
    void apply_template_release()
    {
        for(size_t i=0;i<anchor_weight.size();++i)
            anchor_weight[i]=template_anchor_weights[i]*(released_template_residues.count(nodes[i].res)?.1:1.);
        for(size_t i=0;i<geometry.size();++i)
        {
            const Pair& pair=geometry[i];
            bool release=!pair.validate&&(released_template_residues.count(nodes[pair.a].res)
                ||released_template_residues.count(nodes[pair.b].res));
            geometry[i].weight=template_geometry_weights[i]*(release?.1:1.);
        }
    }
    void reset_template_release()
    {
        if(template_anchor_weights.empty())
        {
            template_anchor_weights=anchor_weight;
            for(const Pair& pair:geometry)template_geometry_weights.push_back(pair.weight);
        }
        released_template_residues.clear();
        apply_template_release();
    }
    void release_clashing_templates()
    {
        refresh(xyz);
        size_t previous=released_template_residues.size();
        for(auto pair:neighbors)
        {
            double overlap=nodes[pair.first].radius+nodes[pair.second].radius-(xyz[pair.first]-xyz[pair.second]).norm();
            if(overlap<=1.)continue;
            for(int atom:{pair.first,pair.second})if(!fixed_records.count(atom))
                for(int residue=nodes[atom].res-1;residue<=nodes[atom].res+1;++residue)
                    if(residue>=1&&residue<=protein.get_end_resno())released_template_residues.insert(residue);
        }
        apply_template_release();
        if(released_template_residues.size()>previous)
            std::cout<<"REEC softened template anchors around "<<released_template_residues.size()
                     <<" residues affected by clashes\n";
    }
    std::string coordinate_signature() const
    {
        uint64_t hash=14695981039346656037ull;
        for(const V& point:xyz)for(double coordinate:{point.x,point.y,point.z})
        {
            uint64_t value=uint64_t(std::llround(coordinate*1000000.));
            for(int byte=0;byte<8;++byte){hash^=(value>>(byte*8))&255;hash*=1099511628211ull;}
        }
        std::ostringstream out;
        out<<std::hex<<hash;
        return out.str();
    }
    void report(const std::string& path,const Stats& before,const Stats& after,bool valid)
    {
        std::ofstream f(path);
        if(!f) throw std::runtime_error("Cannot write report");
        f << std::setprecision(8)
          << "{\n  \"engine\": \"REEC\",\n  \"valid\": " << (valid?"true":"false")
          << ",\n  \"seed\": " << seed
          << ",\n  \"candidate_index\": " << candidate_index
          << ",\n  \"start_signature\": " << std::quoted(start_signature)
          << ",\n  \"clash_scale\": " << clash_scale
          << ",\n  \"annealing\": {\"torsions\": " << sampling.bonds << ", \"trials\": " << sampling.trials
          << ", \"accepted\": " << sampling.accepted << ", \"uphill\": " << sampling.uphill
          << ", \"incremental_improvements\": " << sampling.incremental
          << ", \"score_before\": " << sampling.before << ", \"score_after\": " << sampling.after << "}"
          << ",\n  \"heavy_atoms\": " << std::count_if(nodes.begin(),nodes.end(),[](const Node& n)
        {
            return n.atom->Z>1;
        })
                << ",\n  \"fixed_atoms\": " << fixed_records.size()
                << ",\n  \"dihedral_restraints\": " << dihedral_restraints.size()
                << ",\n  \"mobile_ranges\": " << mobile_ranges.size()
                << ",\n  \"initial_energy\": " << before.score
                << ",\n  \"final_energy\": " << after.score
                << ",\n  \"initial_severe_clashes\": " << before.severe
                << ",\n  \"severe_clashes\": " << after.severe
                << ",\n  \"max_bond_error\": " << after.bond_max
                << ",\n  \"max_geometry_distance_error\": " << after.angle_max
                << ",\n  \"template_ca_rmsd\": " << after.ca_rmsd
                << ",\n  \"max_planarity_volume_error\": " << after.planarity_max
                << ",\n  \"chirality_errors\": " << after.chirality
                << ",\n  \"mandatory_violations\": " << after.hard_fail
                << ",\n  \"geometry_outliers\": [";
        bool comma=false;
        for(auto p:geometry)if(p.validate)
        {
            double actual=(xyz[p.a]-xyz[p.b]).norm(),error=std::fabs(actual-p.lo);
            if(error<=(p.bond?.08:.18))continue;
            if(comma)f<<",";
            comma=true;
            f<<"\n    {\"res1\": "<<nodes[p.a].res<<", \"atom1\": "<<std::quoted(nodes[p.a].name)
             <<", \"res2\": "<<nodes[p.b].res<<", \"atom2\": "<<std::quoted(nodes[p.b].name)
             <<", \"distance\": "<<actual<<", \"ideal\": "<<p.lo<<", \"error\": "<<error<<"}";
        }
        f << "\n  ],\n  \"planarity_outliers\": [";
        comma=false;
        for(auto v:volumes)if(std::fabs(v.value)<.1)
        {
            double actual=(xyz[v.b]-xyz[v.a]).dot((xyz[v.c]-xyz[v.a]).cross(xyz[v.d]-xyz[v.a]));
            double error=std::fabs(actual-v.value);
            if(error<=.25)continue;
            if(comma)f<<",";
            comma=true;
            f<<"\n    {\"residue\": "<<nodes[v.a].res<<", \"atom\": "<<std::quoted(nodes[v.a].name)
             <<", \"error\": "<<error<<"}";
        }
        f << "\n  ]"
                << ",\n  \"released_template_residues\": [";
        comma=false;
        for(int residue:released_template_residues)
        {
            if(comma)f<<",";
            comma=true;
            f<<residue;
        }
        f << "],\n  \"contacts\": [";
        for(size_t i=0; i<contacts.size(); ++i)
        {
            auto& c=contacts[i];
            double err=1e100,distance=0;
            std::pair<int,int> selected;
            for(auto p:c.alternatives)
            {
                double d=(xyz[p.first]-xyz[p.second]).norm();
                double v=std::max({0.,c.lo-d,d-c.hi});
                if(v<err)
                {
                    err=v;
                    distance=d;
                    selected=p;
                }
            }
            if(i) f << ",";
            f << "\n    {\"kind\": " << std::quoted(c.kind)
              << ", \"res1\": " << c.r1 << ", \"atom1\": " << std::quoted(nodes[selected.first].name)
              << ", \"res2\": " << c.r2 << ", \"atom2\": " << std::quoted(nodes[selected.second].name)
              << ", \"distance\": " << distance << ", \"min\": " << c.lo << ", \"max\": " << c.hi
              << ", \"hard\": " << (c.hard?"true":"false") << ", \"violation\": " << err << "}";
        }
        f << "\n  ]\n}\n";
        if(!f) throw std::runtime_error("Report write failed");
    }
    void save(const std::string& path,const Stats& s)
    {
        if(nodes.size()>99999) throw std::runtime_error("Atom count exceeds PDB capacity");
        for(auto p:xyz) if(p.x < -999.999 || p.y < -999.999 || p.z < -999.999 ||
                               p.x > 9999.999 || p.y > 9999.999 || p.z > 9999.999)
                throw std::runtime_error("Coordinates exceed PDB capacity");
        std::string temp=path+".reec-"+std::to_string(getpid());
        FILE* f=std::fopen(temp.c_str(),"w");
        if(!f)throw std::runtime_error("Cannot write PDB: "+path);
        fprintf(f,"REMARK 250 REEC Restrained Ensemble Energy Conformer\nREMARK 250 ENERGY %.6f MAX_BOND_ERROR %.5f HARD_VIOLATIONS %d\n",
                s.score,s.bond_max,s.hard_fail);
        int serial=0;
        std::vector<int> serials;
        for(size_t i=0; i<nodes.size(); ++i)
        {
            auto& n=nodes[i];
            if(fixed_records.count(i))
            {
                serials.push_back(++serial);
                auto line=fixed_records.at(i);
                char number[6];
                snprintf(number,sizeof(number),"%5d",serial);
                line.replace(6,5,number);
                fprintf(f,"%s\n",line.c_str());
                continue;
            }
            n.atom->move(Point(xyz[i].x,xyz[i].y,xyz[i].z));
            serials.push_back(++serial);
            std::string atom=n.name.size()<4?" "+n.name:n.name;
            auto aa=protein.get_residue(n.res);
            fprintf(f,"ATOM  %5d %-4s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n",serial,atom.c_str(),names[codes.find(aa->get_letter())].c_str(),
                    n.res,xyz[i].x,xyz[i].y,xyz[i].z,Atom::esym_from_Z(n.atom->Z));
        }
        for(auto p:disulfides)fprintf(f,"CONECT%5d%5d\n",serials[id(p.first,"SG")],serials[id(p.second,"SG")]);
        for(auto p:fixed_links)fprintf(f,"CONECT%5d%5d\n",serials[p.first],serials[p.second]);
        fprintf(f,"TER\n");
        fprintf(f,"END\n");
        bool failed=std::ferror(f);
        if(std::fclose(f))failed=true;
        if(failed||std::rename(temp.c_str(),path.c_str()))
        {
            std::remove(temp.c_str());
            throw std::runtime_error("PDB publication failed");
        }
    }
};
void help()
{
    std::cout<<"REEC - Restrained Ensemble Energy Conformer\nUsage: bin/reec --alignment FILE.ali --template [ID=]FILE.pdb [--template ...]\n       --output FILE.pdb [--target ID] [--contacts FILE] [--chain A]\n       [--iterations 3000] [--models 1] [--report FILE.json] [--environment FILE.pdb] [--fit-environment]\n       [--environment-hydrogens] [--candidate-prefix PATH] [--seed 0]\n       [--anneal-steps 0] [--anneal-temperature 12] [--anneal-residues 103,107]\nAnnealing is off by default; --anneal-steps 64 enables it.\nRun from the Aroma repository root. Exit 2 means validation failed; no PDB is published.\nContacts: disulfide RES RES\n          contact hard|soft RES ATOM RES ATOM MIN MAX WEIGHT\n          dihedral RES ATOM RES ATOM RES ATOM RES ATOM MEAN STDEV (radians)\n          mobile FIRST LAST\nAtom alternatives use |, for example OD1|OD2. Distances are Angstroms.\n";
}
}
#ifndef REEC_NO_MAIN
int main(int argc,char** argv)
{
    using namespace reec;
    try
    {
        std::cout.setf(std::ios::unitbuf);
        std::string ali,out,target,report,environment,candidate_prefix;
        bool fit_environment=false,environment_hydrogens=false;
        std::vector<std::string> templates,contact_files;
        char chain='A';
        int iterations=3000,models=1,seed=0,anneal_steps=0;
        double anneal_temperature=12,clash_scale=1;
        std::set<int> anneal_residues;
        for(int i=1; i<argc; ++i)
        {
            std::string a=argv[i];
            if(a=="--environment-hydrogens")
            {
                environment_hydrogens=true;
                continue;
            }
            if(a=="--fit-environment")
            {
                fit_environment=true;
                continue;
            }
            if(a=="--help"||a=="-h")
            {
                help();
                return 0;
            }
            if(i+1>=argc)throw std::runtime_error("Missing option value: "+a);
            std::string v=argv[++i];
            if(a=="--alignment")ali=v;
            else if(a=="--output")out=v;
            else if(a=="--template")templates.push_back(v);
            else if(a=="--target")target=v;
            else if(a=="--contacts")contact_files.push_back(v);
            else if(a=="--report")report=v;
            else if(a=="--candidate-prefix")candidate_prefix=v;
            else if(a=="--environment")environment=v;
            else if(a=="--chain")
            {
                if(v.size()!=1)throw std::runtime_error("Chain must be one character");
                chain=v[0];
            }
            else if(a=="--iterations")iterations=integer(v);
            else if(a=="--models")models=integer(v);
            else if(a=="--seed")seed=integer(v);
            else if(a=="--anneal-steps")anneal_steps=integer(v);
            else if(a=="--anneal-temperature")anneal_temperature=number(v);
            else if(a=="--clash-scale")clash_scale=number(v);
            else if(a=="--anneal-residues")for(auto residue:split(v,','))
            {
                int r=integer(residue);
                if(r<1)throw std::runtime_error("Annealing residues must be positive");
                anneal_residues.insert(r);
            }
            else throw std::runtime_error("Unknown option: "+a);
        }
        if(ali.empty()||out.empty()||templates.empty())
        {
            help();
            return 1;
        }
        if(iterations<0||iterations>100000||models<1||models>32)throw std::runtime_error("Invalid iteration or model count");
        if(seed<0)throw std::runtime_error("Seed must be nonnegative");
        if(anneal_steps<0||anneal_steps>100000||anneal_temperature<=0)throw std::runtime_error("Invalid annealing settings");
        if(clash_scale<=0)throw std::runtime_error("Clash scale must be positive");
        if(report.empty())report=out+".json";
        if(report==out)throw std::runtime_error("Report and PDB paths must differ");
        auto entries=alignment(ali);
        Entry* chosen=nullptr;
        for(auto& e:entries)if(!e.structure&&(target.empty()||e.id==target))
            {
                if(chosen)throw std::runtime_error("Select one target with --target");
                chosen=&e;
            }
        if(!chosen)throw std::runtime_error("Target absent from PIR");
        protein_row(*chosen);
        std::string seq;
        for(char c:chosen->seq)if(c!='-')seq+=c;
        if(seq.size()<2||seq.size()>5000)throw std::runtime_error("Target length outside 2..5000");
        std::vector<Template> t;
        for(auto& path:templates)
        {
            t.push_back(pdb(path,chain));
            map_template(t.back(),entries,*chosen,seq,templates.size()==1);
        }
        for(size_t i=1; i<t.size(); ++i)superpose(t[i],t[0]);
        Model m;
        m.seed=unsigned(seed);
        m.clash_scale=clash_scale;
        for(int r:anneal_residues)if(r>int(seq.size()))throw std::runtime_error("Annealing residue exceeds target length");
        m.anneal_residues=anneal_residues;
        m.build(seq,t);
        if(!environment.empty())
        {
            RigidFit fit;
            if(fit_environment)
            {
                Template env=pdb(environment,chain);
                env.mapping=align(seq,env);
                fit=superpose(env,t[0]);
            }
            m.environment(environment,fit,environment_hydrogens);
        }
        for(auto& path:contact_files)m.read_contacts(path);
        Stats initial=m.stats();
        XYZ start=m.xyz,best;
        double best_score=1e100;
        Stats final;
        Model::Annealing best_sampling;
        int best_candidate=0;
        std::string best_signature;
        std::set<int> best_released;
        bool valid=false;
        for(int candidate=0; candidate<models; ++candidate)
        {
            m.xyz=start;
            m.reset_template_release();
            if(candidate||seed) m.alternate_start(candidate);
            m.candidate_index=candidate;
            m.start_signature=m.coordinate_signature();
            std::cout<<"REEC candidate "<<candidate+1<<"/"<<models<<" seed "<<seed<<" start "<<m.start_signature<<"\n";
            for(int stage=0; stage<3; ++stage)
            {
                if(stage)m.release_clashing_templates();
                m.hard_scale=stage==0?10:(stage==1?100:1000);
                m.anchor_scale=stage==0?1:(stage==1?.3:.1);
                m.minimize(iterations/3+(stage==2?iterations%3:0));
                if(stage==0)m.anneal(anneal_steps,anneal_temperature,candidate);
            }
            Stats s=m.stats();
            bool ok=s.hard_fail==0&&s.chirality==0&&s.bond_max<=.08&&s.angle_max<=.18&&s.planarity_max<=.25&&s.severe==0;
            if(!candidate_prefix.empty())
            {
                std::string path=candidate_prefix+"."+std::to_string(candidate)+".pdb";
                m.report(path+".json",initial,s,ok);
                if(ok)m.save(path,s);
            }
            std::cout<<"Candidate energy "<<s.score<<" bond error "<<s.bond_max<<" geometry error "<<s.angle_max<<" severe clashes "
                     <<s.severe<<" mandatory violations "<<s.hard_fail<<"\n";
            if((ok&&!valid)||(ok==valid&&s.score<best_score))
            {
                best=m.xyz;
                best_score=s.score;
                final=s;
                best_sampling=m.sampling;
                best_candidate=candidate;
                best_signature=m.start_signature;
                best_released=m.released_template_residues;
                valid=ok;
            }
        }
        m.xyz=best;
        m.sampling=best_sampling;
        m.candidate_index=best_candidate;
        m.start_signature=best_signature;
        m.released_template_residues=best_released;
        m.apply_template_release();
        m.report(report,initial,final,valid);
        if(!valid)
        {
            std::cerr<<"REEC validation failed; inspect "<<report<<". Existing output remains intact.\n";
            return 2;
        }
        m.save(out,final);
        std::cout<<"REEC wrote "<<out<<" and "<<report<<"\n";
        return 0;
    }
    catch(const std::exception& e)
    {
        std::cerr<<"REEC: "<<e.what()<<"\n";
        return 1;
    }
    catch(...)
    {
        std::cerr<<"REEC: molecular library error\n";
        return 1;
    }
}
#endif
