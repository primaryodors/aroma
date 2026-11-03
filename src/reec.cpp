// REEC: Restrained Ensemble Energy Conformer.
// Geometry comes from the project's amino acid library; optimization is local.
#include "classes/protein.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <vector>
#include <cstdio>
#include <unistd.h>

namespace reec {
struct V {
    double x=0,y=0,z=0;
    V() = default;
    V(double a,double b,double c):x(a),y(b),z(c) {}
    V(const Point& p):x(p.x),y(p.y),z(p.z) {}
    V operator+(V b) const { return {x+b.x,y+b.y,z+b.z}; }
    V operator-(V b) const { return {x-b.x,y-b.y,z-b.z}; }
    V operator*(double b) const { return {x*b,y*b,z*b}; }
    V& operator+=(V b) { x+=b.x;y+=b.y;z+=b.z;return *this; }
    double dot(V b) const { return x*b.x+y*b.y+z*b.z; }
    V cross(V b) const { return {y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x}; }
    double norm() const { return std::sqrt(dot(*this)); }
    V unit() const { double d=norm(); if(d<1e-8) throw std::runtime_error("Degenerate backbone frame"); return *this*(1/d); }
};
using XYZ=std::vector<V>;
double dot(const XYZ& a,const XYZ& b) { double r=0;for(size_t i=0;i<a.size();++i) r+=a[i].dot(b[i]);return r; }
std::string trim(std::string s) { auto a=s.find_first_not_of(" \r\n\t");return a==s.npos?"":s.substr(a,s.find_last_not_of(" \r\n\t")-a+1); }
std::vector<std::string> split(std::string s,char c) { std::vector<std::string> v;std::stringstream f(s);while(std::getline(f,s,c))v.push_back(trim(s));return v; }
int integer(const std::string& s) { size_t p=0;int n=std::stoi(s,&p);if(p!=s.size())throw std::runtime_error("Invalid integer: "+s);return n; }
double number(const std::string& s) { size_t p=0;double n=std::stod(s,&p);if(p!=s.size()||!std::isfinite(n))throw std::runtime_error("Invalid number: "+s);return n; }
const std::string codes="ACDEFGHIKLMNPQRSTVWY";
const std::vector<std::string> names={"ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR"};
char letter(std::string s) { for(size_t i=0;i<names.size();++i)if(s==names[i])return codes[i];if(s=="MSE")return 'M';return 0; }
struct Entry { std::string id,desc,seq; bool structure=false; };
std::vector<Entry> alignment(const std::string& path) {
    std::ifstream f(path);if(!f)throw std::runtime_error("Cannot read alignment: "+path);
    std::vector<Entry> a;std::string s;bool desc=false,done=true;
    while(std::getline(f,s)) {
        s=trim(s);if(s.empty())continue;
        if(s[0]=='>') { if(!done)throw std::runtime_error("Unterminated PIR entry");auto p=s.find(';');if(p==s.npos)throw std::runtime_error("Expected PIR header >P1;ID");a.push_back({s.substr(p+1),"","",false});desc=true;done=false;continue; }
        if(a.empty()||done)throw std::runtime_error("Unexpected PIR content");
        if(desc){a.back().desc=s;a.back().structure=s.compare(0,9,"structure")==0;if(!a.back().structure&&s.compare(0,8,"sequence")!=0)throw std::runtime_error("Expected sequence or structure descriptor");desc=false;continue;}
        for(char c:s){if(std::isspace(static_cast<unsigned char>(c)))continue;if(c=='*'){done=true;break;}if(c=='/'||c=='.')throw std::runtime_error("Only a single protein chain is supported in PIR; omit ligand columns");if(c!='-'&&codes.find(c)==codes.npos)throw std::runtime_error("Unsupported residue in PIR");a.back().seq+=c;}
    }
    if(!done||a.empty())throw std::runtime_error("Incomplete PIR alignment");
    std::set<std::string> ids;for(auto& e:a)if(e.seq.empty()||!ids.insert(e.id).second)throw std::runtime_error("Empty or duplicate PIR entry");return a;
}
struct Res { int no;char aa;std::map<std::string,V> atoms; };
struct Template { std::string path,id;std::vector<Res> residues;std::set<std::pair<int,int>> links;std::vector<int> mapping; };
Template pdb(const std::string& arg,char chain) {
    Template t;auto eq=arg.find('=');t.path=eq==arg.npos?arg:arg.substr(eq+1);t.id=eq==arg.npos?"":arg.substr(0,eq);
    std::ifstream f(t.path);if(!f)throw std::runtime_error("Cannot read template: "+t.path);
    std::string s;std::map<int,int> serial;std::vector<std::pair<int,int>> conect;std::map<int,int> index;bool ended=false;
    while(std::getline(f,s)) {
        if(s.compare(0,6,"ENDMDL")==0){ended=true;continue;}
        if(!ended&&s.size()>=54&&(s.compare(0,6,"ATOM  ")==0||s.compare(0,6,"HETATM")==0)&&s[21]==chain) {
            if(s[16]!=' '&&s[16]!='A')continue;
            char aa=letter(s.substr(17,3));if(!aa)continue;
            if(s[26]!=' ')throw std::runtime_error("Insertion codes require a renumbered template: "+t.path);
            int n=integer(trim(s.substr(22,4)));
            if(!index.count(n)){index[n]=t.residues.size();t.residues.push_back({n,aa,{}});}
            Res& r=t.residues[index[n]];std::string atom=trim(s.substr(12,4));if(s.substr(17,3)=="MSE"&&atom=="SE")atom="SD";
            V p{number(trim(s.substr(30,8))),number(trim(s.substr(38,8))),number(trim(s.substr(46,8)))};
            if(!r.atoms.count(atom)||s[16]==' ')r.atoms[atom]=p;
            if(atom=="SG"&&aa=='C')serial[integer(trim(s.substr(6,5)))]=n;
        }
        if(s.compare(0,6,"SSBOND")==0&&s.size()>=35&&s[15]==chain&&s[29]==chain){int a=integer(trim(s.substr(17,4))),b=integer(trim(s.substr(31,4)));t.links.insert(std::minmax(a,b));}
        if(s.compare(0,6,"CONECT")==0&&s.size()>=11){int a=integer(trim(s.substr(6,5)));for(size_t i=11;i<s.size();i+=5){auto v=trim(s.substr(i,5));if(!v.empty())conect.push_back({a,integer(v)});}}
    }
    for(auto p:conect)if(serial.count(p.first)&&serial.count(p.second)&&p.first!=p.second)t.links.insert(std::minmax(serial[p.first],serial[p.second]));
    if(t.residues.empty())throw std::runtime_error("No protein atoms in requested template chain");
    for(auto& r:t.residues)for(auto a:{"N","CA","C"})if(!r.atoms.count(a))throw std::runtime_error("Template backbone incomplete at "+std::to_string(r.no));
    return t;
}
// Global sequence alignment for a target-only PIR file. Explicit PIR columns take priority.
std::vector<int> align(const std::string& target,const Template& t) {
    int n=target.size(),m=t.residues.size();std::vector<int> d((n+1)*(m+1));std::vector<char> trace(d.size());
    for(int i=1;i<=n;++i){d[i*(m+1)]=-3*i;trace[i*(m+1)]=1;}for(int j=1;j<=m;++j){d[j]=-3*j;trace[j]=2;}
    for(int i=1;i<=n;++i)for(int j=1;j<=m;++j){int k=i*(m+1)+j;int v=d[k-m-2]+(target[i-1]==t.residues[j-1].aa?4:-1);char c=0;if(d[k-m-1]-3>v){v=d[k-m-1]-3;c=1;}if(d[k-1]-3>v){v=d[k-1]-3;c=2;}d[k]=v;trace[k]=c;}
    std::vector<int> map(n,-1);int i=n,j=m;while(i||j){char c=trace[i*(m+1)+j];if(!c){map[--i]=--j;}else if(c==1)--i;else --j;}return map;
}
void map_template(Template& t,const std::vector<Entry>& entries,const Entry& target,const std::string& seq) {
    const Entry* e=nullptr;
    for(auto& a:entries)if(a.structure&&(a.id==t.id||t.path==a.id+".pdb"||t.path.substr(t.path.find_last_of('/')+1)==a.id+".pdb")){if(e)throw std::runtime_error("Ambiguous template alignment");e=&a;}
    if(!e&&t.id.empty()){int count=0;for(auto& a:entries)if(a.structure){e=&a;++count;}if(count!=1)e=nullptr;}
    if(!e){if(!t.id.empty())throw std::runtime_error("Template ID absent from PIR: "+t.id);t.mapping=align(seq,t);return;}
    if(e->seq.size()!=target.seq.size())throw std::runtime_error("PIR rows have different column counts");
    auto fields=split(e->desc,':');int start=0;
    if(fields.size()>2&&fields[2]!="FIRST"&&!fields[2].empty()){int n=integer(fields[2]);while(start<int(t.residues.size())&&t.residues[start].no<n)++start;}
    int j=start;for(size_t c=0;c<target.seq.size();++c){int k=-1;if(e->seq[c]!='-'){if(j>=int(t.residues.size())||e->seq[c]!=t.residues[j].aa)throw std::runtime_error("PIR/template sequence mismatch in "+t.path+" at column "+std::to_string(c+1));k=j++;}if(target.seq[c]!='-')t.mapping.push_back(k);}
    if(j!=int(t.residues.size()))throw std::runtime_error("Template contains residues outside PIR row");
}
struct Frame {
    V o,x,y,z;
    Frame(V n,V ca,V c):o(ca){x=(c-ca).unit();z=x.cross(n-ca).unit();y=z.cross(x);}
    V local(V p)const{p=p-o;return {p.dot(x),p.dot(y),p.dot(z)};}
    V world(V p)const{return o+x*p.x+y*p.y+z*p.z;}
};
// Horn's quaternion fit, solved by Jacobi diagonalization, without an external library.
void superpose(Template& moving,const Template& reference) {
    XYZ a,b;for(size_t i=0;i<moving.mapping.size();++i)if(moving.mapping[i]>=0&&reference.mapping[i]>=0){a.push_back(moving.residues[moving.mapping[i]].atoms.at("CA"));b.push_back(reference.residues[reference.mapping[i]].atoms.at("CA"));}
    if(a.size()<3)throw std::runtime_error("Templates have fewer than three common CA anchors");
    V ca,cb;for(size_t i=0;i<a.size();++i){ca+=a[i];cb+=b[i];}ca=ca*(1./a.size());cb=cb*(1./a.size());
    double s[3][3]={};for(size_t i=0;i<a.size();++i){V u=a[i]-ca,v=b[i]-cb;double x[]={u.x,u.y,u.z},y[]={v.x,v.y,v.z};for(int j=0;j<3;++j)for(int k=0;k<3;++k)s[j][k]+=x[j]*y[k];}
    double q[4][4]={{s[0][0]+s[1][1]+s[2][2],s[1][2]-s[2][1],s[2][0]-s[0][2],s[0][1]-s[1][0]},
        {s[1][2]-s[2][1],s[0][0]-s[1][1]-s[2][2],s[0][1]+s[1][0],s[0][2]+s[2][0]},
        {s[2][0]-s[0][2],s[0][1]+s[1][0],-s[0][0]+s[1][1]-s[2][2],s[1][2]+s[2][1]},
        {s[0][1]-s[1][0],s[0][2]+s[2][0],s[1][2]+s[2][1],-s[0][0]-s[1][1]+s[2][2]}};
    double ev[4][4]={};for(int i=0;i<4;++i)ev[i][i]=1;
    for(int it=0;it<60;++it){int p=0,r=1;for(int j=0;j<4;++j)for(int k=j+1;k<4;++k)if(std::fabs(q[j][k])>std::fabs(q[p][r])){p=j;r=k;}if(std::fabs(q[p][r])<1e-12)break;double angle=.5*std::atan2(2*q[p][r],q[r][r]-q[p][p]),c=std::cos(angle),sine=std::sin(angle);double pp=q[p][p],rr=q[r][r],pr=q[p][r];q[p][p]=c*c*pp-2*c*sine*pr+sine*sine*rr;q[r][r]=sine*sine*pp+2*c*sine*pr+c*c*rr;q[p][r]=q[r][p]=0;for(int k=0;k<4;++k){if(k!=p&&k!=r){double kp=q[k][p],kr=q[k][r];q[k][p]=q[p][k]=c*kp-sine*kr;q[k][r]=q[r][k]=sine*kp+c*kr;}double kp=ev[k][p],kr=ev[k][r];ev[k][p]=c*kp-sine*kr;ev[k][r]=sine*kp+c*kr;}}
    int best=0;for(int i=1;i<4;++i)if(q[i][i]>q[best][best])best=i;double w=ev[0][best];V v{ev[1][best],ev[2][best],ev[3][best]};
    for(auto& r:moving.residues)for(auto& atom:r.atoms){V p=atom.second-ca;atom.second=p+v.cross(p)*(2*w)+v.cross(v.cross(p))*2+cb;}
}
struct Node { Atom* atom;int res;std::string name;double radius;bool backbone; };
struct Pair { int a,b;double lo,hi,weight;bool bond=false; };
struct Volume { int a,b,c,d;double value,weight; };
struct Contact { int r1,r2;std::string a1,a2,kind;double lo,hi,weight;bool hard;std::vector<std::pair<int,int>> alternatives; };
struct Stats { double score=0,bond_max=0,angle_max=0,overlap_max=0;int severe=0,chirality=0,hard_fail=0; };
struct Model {
    Protein protein;
    std::vector<Node> nodes;
    XYZ xyz,ideal,anchors;
    std::vector<double> anchor_weight;
    std::vector<Pair> geometry;
    std::vector<Volume> volumes;
    std::vector<Contact> contacts;
    std::map<std::pair<int,std::string>,int> ids;
    std::vector<std::vector<int>> adj;
    std::set<std::pair<int,int>> excluded;
    std::vector<std::pair<int,int>> neighbors;
    std::vector<std::pair<int,int>> disulfides;
    XYZ neighbor_origin;
    double hard_scale=1,anchor_scale=1;
    int id(int r,const std::string& a)const{auto it=ids.find({r,a});if(it==ids.end())throw std::runtime_error("Unknown atom "+std::to_string(r)+":"+a);return it->second;}
    void pair(int a,int b,double target,double weight,bool bond=false){geometry.push_back({a,b,target,target,weight,bond});}
    void pair(int r,const char* a,int s,const char* b,double distance,double w,bool bond=false){pair(id(r,a),id(s,b),distance,w,bond);}
    void edge(int a,int b){adj[a].push_back(b);adj[b].push_back(a);}
    void build(const std::string& seq,std::vector<Template>& templates) {
        // Use Protein/Atom topology directly; no separate residue dictionary.
        for(size_t r=0;r<seq.size();++r){if(!protein.add_residue(r+1,seq[r]))throw std::runtime_error("Cannot create residue");AminoAcid* aa=protein.get_residue(r+1);for(int j=0;j<aa->get_atom_count();++j){Atom* a=aa->get_atom(j);if(a->Z<=1||std::string(a->name)=="OXT")continue;int i=nodes.size();ids[{int(r+1),a->name}]=i;nodes.push_back({a,int(r+1),a->name,a->vdW_radius,a->is_backbone});ideal.emplace_back(a->loc);}}
        xyz=ideal;adj.resize(nodes.size());
        for(size_t i=0;i<nodes.size();++i)for(size_t j=i+1;j<nodes.size();++j)if(nodes[i].res==nodes[j].res&&nodes[i].atom->is_bonded_to(nodes[j].atom)){edge(i,j);double d=(ideal[i]-ideal[j]).norm();if(d<.9||d>2.2)throw std::runtime_error("Invalid amino acid library geometry");pair(i,j,d,3000,true);}
        // Distances over two bonds preserve bond angles; rigid central bonds preserve planarity.
        for(size_t b=0;b<adj.size();++b)for(size_t j=0;j<adj[b].size();++j)for(size_t k=j+1;k<adj[b].size();++k){int a=adj[b][j],c=adj[b][k];pair(a,c,(ideal[a]-ideal[c]).norm(),1200);}
        std::set<std::pair<int,int>> rigid;
        for(size_t b=0;b<adj.size();++b)for(int c:adj[b])if(int(b)<c){Bond* bond=nodes[b].atom->get_bond_between(nodes[c].atom);if(bond&&!bond->can_rotate)for(int a:adj[b])for(int d:adj[c])if(a!=c&&d!=int(b)&&a!=d){auto p=std::minmax(a,d);if(rigid.insert(p).second)pair(a,d,(ideal[a]-ideal[d]).norm(),600);}}
        for(size_t a=0;a<adj.size();++a)if(adj[a].size()>=3){int b=adj[a][0],c=adj[a][1],d=adj[a][2];double v=(ideal[b]-ideal[a]).dot((ideal[c]-ideal[a]).cross(ideal[d]-ideal[a]));volumes.push_back({int(a),b,c,d,v,100});}
        const int n=seq.size();
        // Standard trans peptide geometry: 1.329 A C-N, 116.2/121.7 degree angles.
        const double pi=std::acos(-1.),a=116.2*pi/180,b=121.7*pi/180,o=123.0*pi/180;
        V C{0,0,0},N{1.329,0,0},CA{1.525*std::cos(a),1.525*std::sin(a),0};
        V next=N+V{-1.458*std::cos(b),-1.458*std::sin(b),0},O{1.231*std::cos(o),-1.231*std::sin(o),0};
        for(int r=1;r<n;++r){int c=id(r,"C"),nn=id(r+1,"N");edge(c,nn);pair(c,nn,1.329,3000,true);pair(r,"CA",r+1,"N",(CA-N).norm(),1200);pair(r,"O",r+1,"N",(O-N).norm(),1200);pair(r,"C",r+1,"CA",next.norm(),1200);pair(r,"CA",r+1,"CA",(CA-next).norm(),800);pair(r,"O",r+1,"CA",(O-next).norm(),800);volumes.push_back({c,id(r,"CA"),id(r,"O"),nn,0,100});}
        for(size_t aidx=0;aidx<adj.size();++aidx){excluded.insert({aidx,aidx});for(int bidx:adj[aidx]){excluded.insert(std::minmax(int(aidx),bidx));for(int cidx:adj[bidx]){excluded.insert(std::minmax(int(aidx),cidx));for(int didx:adj[cidx])excluded.insert(std::minmax(int(aidx),didx));}}}
        std::vector<int> source(n,-1),mapped(n,-1);
        for(int r=0;r<n;++r){double best=-1;for(size_t t=0;t<templates.size();++t){int k=templates[t].mapping[r];if(k<0)continue;double quality=templates[t].residues[k].aa==seq[r]?2:1;if(quality>best){best=quality;source[r]=t;mapped[r]=k;}}}
        int count=std::count_if(source.begin(),source.end(),[](int s){return s>=0;});if(!count)throw std::runtime_error("No aligned template residues");
        std::vector<V> cas(n);for(int r=0;r<n;++r)if(source[r]>=0)cas[r]=templates[source[r]].residues[mapped[r]].atoms.at("CA");
        // Place untemplated loops along an arched bridge. Optimization closes all peptide bonds.
        for(int first=0;first<n;){if(source[first]>=0){++first;continue;}int last=first;while(last<n&&source[last]<0)++last;V left,right,normal;
            if(first>0&&last<n){left=cas[first-1];right=cas[last];V tangent=(right-left).unit();normal=tangent.cross({0,0,1});if(normal.norm()<.1)normal=tangent.cross({0,1,0});normal=normal.unit();double span=(right-left).norm(),length=3.8*(last-first+1),height=std::min(8.,std::sqrt(std::max(0.,length*length-span*span))*.3);for(int r=first;r<last;++r){double u=double(r-first+1)/(last-first+1);cas[r]=left*(1-u)+right*u+normal*(height*std::sin(pi*u));}}
            else {int anchor=first>0?first-1:last;V dir{1,0,0};if(first>1)dir=(cas[first-1]-cas[first-2]).unit();else if(last+1<n)dir=(cas[last]-cas[last+1]).unit();for(int r=first;r<last;++r)cas[r]=cas[anchor]+dir*(3.8*std::abs(r-anchor));}first=last;}
        anchor_weight.resize(nodes.size());
        for(int r=0;r<n;++r){Frame original(ideal[id(r+1,"N")],ideal[id(r+1,"CA")],ideal[id(r+1,"C")]);V nn,cc;
            if(source[r]>=0){auto& rr=templates[source[r]].residues[mapped[r]];nn=rr.atoms.at("N");cc=rr.atoms.at("C");}
            else {V tangent=(r+1<n?cas[r+1]-cas[r]:cas[r]-cas[r-1]).unit();V norm=tangent.cross({0,0,1});if(norm.norm()<.1)norm=tangent.cross({0,1,0});norm=norm.unit();nn=cas[r]-tangent*.9+norm*1.15;cc=cas[r]+tangent*1.3+norm*.8;}
            Frame frame(nn,cas[r],cc);
            for(size_t i=0;i<nodes.size();++i)if(nodes[i].res==r+1){xyz[i]=frame.world(original.local(ideal[i]));if(source[r]>=0){auto& rr=templates[source[r]].residues[mapped[r]];auto p=rr.atoms.find(nodes[i].name);if(p!=rr.atoms.end()&&(nodes[i].backbone||rr.aa==seq[r]))xyz[i]=p->second;anchor_weight[i]=nodes[i].backbone?1.5:.03;}}
        }
        anchors=xyz;
        for(auto& t:templates)for(auto link:t.links){int ares=0,bres=0;for(int r=0;r<n;++r)if(t.mapping[r]>=0){int no=t.residues[t.mapping[r]].no;if(no==link.first)ares=r+1;if(no==link.second)bres=r+1;}if(ares&&bres&&seq[ares-1]=='C'&&seq[bres-1]=='C')add_disulfide(ares,bres);}
        std::cout<<"REEC built "<<n<<" residues, "<<nodes.size()<<" heavy atoms; template coverage "<<count<<"/"<<n<<"\n";
    }
    void resolve(Contact& c){auto aa=split(c.a1,'|'),bb=split(c.a2,'|');for(auto& a:aa)for(auto& b:bb)c.alternatives.push_back({id(c.r1,a),id(c.r2,b)});for(auto p:c.alternatives)if(p.first==p.second)throw std::runtime_error("Self contact");}
    void add_disulfide(int a,int b){std::pair<int,int> p=std::minmax(a,b);for(auto link:disulfides)if(link==p)return;for(auto link:disulfides)if(link.first==a||link.second==a||link.first==b||link.second==b)throw std::runtime_error("Conflicting disulfide topology");Contact c{a,b,"SG","SG","disulfide",1.95,2.15,100,true,{}};resolve(c);contacts.push_back(c);disulfides.push_back(p);int x=id(a,"SG"),y=id(b,"SG");edge(x,y);excluded.insert(std::minmax(x,y));for(int v:adj[x])excluded.insert(std::minmax(v,y));for(int v:adj[y])excluded.insert(std::minmax(v,x));pair(id(a,"CB"),y,2.80,500);pair(id(b,"CB"),x,2.80,500);}
    void read_contacts(const std::string& path){std::ifstream f(path);if(!f)throw std::runtime_error("Cannot read contact file: "+path);std::string s;int line=0;while(std::getline(f,s)){++line;s=s.substr(0,s.find('#'));s=trim(s);if(s.empty())continue;try {std::stringstream ss(s);std::string kind;ss>>kind;if(kind=="disulfide"){std::string a,b,extra;if(!(ss>>a>>b)||(ss>>extra))throw std::runtime_error("Expected disulfide RES RES");add_disulfide(integer(a),integer(b));continue;}Contact c;std::string mode,r1,r2,lo,hi,w,extra;if(!(ss>>mode>>r1>>c.a1>>r2>>c.a2>>lo>>hi>>w)||(ss>>extra))throw std::runtime_error("Expected contact hard|soft RES ATOM RES ATOM MIN MAX WEIGHT");if(kind!="contact"||(mode!="hard"&&mode!="soft"))throw std::runtime_error("Invalid contact type");c.kind="contact";c.hard=mode=="hard";c.r1=integer(r1);c.r2=integer(r2);c.lo=number(lo);c.hi=number(hi);c.weight=number(w);if(c.lo<=0||c.hi<c.lo||c.weight<=0)throw std::runtime_error("Invalid contact bounds or weight");resolve(c);contacts.push_back(c);}catch(const std::exception& e){throw std::runtime_error(path+":"+std::to_string(line)+": "+e.what());}}}
    void refresh(const XYZ& x){neighbors.clear();std::map<std::array<int,3>,std::vector<int>> cells;double cell=6.;for(size_t i=0;i<x.size();++i){std::array<int,3> key{{int(std::floor(x[i].x/cell)),int(std::floor(x[i].y/cell)),int(std::floor(x[i].z/cell))}};for(int a=-1;a<=1;++a)for(int b=-1;b<=1;++b)for(int c=-1;c<=1;++c){auto it=cells.find({{key[0]+a,key[1]+b,key[2]+c}});if(it!=cells.end())for(int j:it->second)if(!excluded.count({j,int(i)})&&(x[i]-x[j]).norm()<6.)neighbors.push_back({j,int(i)});}cells[key].push_back(i);}neighbor_origin=x;}
    bool stale(const XYZ& x)const {if(neighbor_origin.size()!=x.size())return true;for(size_t i=0;i<x.size();++i)if((x[i]-neighbor_origin[i]).norm()>.6)return true;return false;}
    static double distance(const XYZ& x,XYZ& g,int a,int b,double lo,double hi,double w){V v=x[a]-x[b];double d=std::max(1e-10,v.norm()),err=d<lo?d-lo:(d>hi?d-hi:0);V force=v*(2*w*err/d);if(d<1e-8)force={2*w*err,0,0};g[a]+=force;g[b]+=force*(-1);return w*err*err;}
    double energy(const XYZ& x,XYZ& g){g.assign(x.size(),{});double e=0;for(auto& p:geometry)e+=distance(x,g,p.a,p.b,p.lo,p.hi,p.weight);
        for(auto& v:volumes){V b=x[v.b]-x[v.a],c=x[v.c]-x[v.a],d=x[v.d]-x[v.a];double err=b.dot(c.cross(d))-v.value,w=v.weight;e+=w*err*err;V gb=c.cross(d)*(2*w*err),gc=d.cross(b)*(2*w*err),gd=b.cross(c)*(2*w*err);g[v.b]+=gb;g[v.c]+=gc;g[v.d]+=gd;g[v.a]+=(gb+gc+gd)*(-1);}
        for(size_t i=0;i<x.size();++i){V v=x[i]-anchors[i];double w=anchor_weight[i]*anchor_scale;e+=w*v.dot(v);g[i]+=v*(2*w);}
        for(auto& c:contacts){double best=1e100;std::pair<int,int> chosen;for(auto p:c.alternatives){double d=(x[p.first]-x[p.second]).norm(),err=d<c.lo?c.lo-d:(d>c.hi?d-c.hi:0);if(err<best){best=err;chosen=p;}}e+=distance(x,g,chosen.first,chosen.second,c.lo,c.hi,c.weight*(c.hard?hard_scale:1));}
        for(auto p:neighbors){double limit=nodes[p.first].radius+nodes[p.second].radius-.45;if((nodes[p.first].atom->Z==7||nodes[p.first].atom->Z==8)&&(nodes[p.second].atom->Z==7||nodes[p.second].atom->Z==8))limit=std::min(limit,2.65);e+=distance(x,g,p.first,p.second,limit,1e10,30);}
        if(!std::isfinite(e))throw std::runtime_error("Nonfinite optimization energy");return e;
    }
    void minimize(int iterations){refresh(xyz);XYZ g;double e=energy(xyz,g);std::vector<XYZ> sh,yh;std::vector<double> rho;
        for(int it=0;it<iterations;++it){if(stale(xyz)){refresh(xyz);e=energy(xyz,g);sh.clear();yh.clear();rho.clear();}XYZ direction=g;std::vector<double> alpha(sh.size());for(int h=int(sh.size())-1;h>=0;--h){alpha[h]=rho[h]*dot(sh[h],direction);for(size_t i=0;i<g.size();++i)direction[i]+=yh[h][i]*(-alpha[h]);}double scale=sh.empty()?1./6000:dot(sh.back(),yh.back())/std::max(1e-20,dot(yh.back(),yh.back()));for(auto& v:direction)v=v*scale;for(size_t h=0;h<sh.size();++h){double beta=rho[h]*dot(yh[h],direction);for(size_t i=0;i<g.size();++i)direction[i]+=sh[h][i]*(alpha[h]-beta);}for(auto& v:direction)v=v*(-1);
            double slope=dot(g,direction);if(slope>=0){sh.clear();yh.clear();rho.clear();direction=g;for(auto& v:direction)v=v*(-1./6000);slope=dot(g,direction);}double maxstep=0;for(auto& v:direction)maxstep=std::max(maxstep,v.norm());if(maxstep<1e-9)break;double step=std::min(1.,.15/maxstep);XYZ trial,ng;double en=e;bool accepted=false;for(int ls=0;ls<25;++ls){trial=xyz;for(size_t i=0;i<g.size();++i)trial[i]+=direction[i]*step;en=energy(trial,ng);if(en<=e+1e-4*step*slope){accepted=true;break;}step*=.5;}if(!accepted)break;
            XYZ s(g.size()),y(g.size());for(size_t i=0;i<g.size();++i){s[i]=trial[i]-xyz[i];y[i]=ng[i]-g[i];}double sy=dot(s,y);if(sy>1e-10){if(sh.size()==10){sh.erase(sh.begin());yh.erase(yh.begin());rho.erase(rho.begin());}sh.push_back(s);yh.push_back(y);rho.push_back(1/sy);}xyz.swap(trial);g.swap(ng);e=en;
            if(it%250==0)std::cout<<"  iteration "<<it<<" energy "<<e<<"\n";
        }
    }
    Stats stats(){refresh(xyz);XYZ g;Stats s;s.score=energy(xyz,g);for(auto p:geometry){double err=std::fabs((xyz[p.a]-xyz[p.b]).norm()-p.lo);if(p.bond)s.bond_max=std::max(s.bond_max,err);else s.angle_max=std::max(s.angle_max,err);}for(auto v:volumes)if(std::fabs(v.value)>.5){double actual=(xyz[v.b]-xyz[v.a]).dot((xyz[v.c]-xyz[v.a]).cross(xyz[v.d]-xyz[v.a]));if(actual*v.value<=0)++s.chirality;}for(auto p:neighbors){double overlap=nodes[p.first].radius+nodes[p.second].radius-(xyz[p.first]-xyz[p.second]).norm();s.overlap_max=std::max(s.overlap_max,overlap);if(overlap>1.0)++s.severe;}for(auto& c:contacts){double err=1e100;for(auto p:c.alternatives){double d=(xyz[p.first]-xyz[p.second]).norm();err=std::min(err,std::max({0.,c.lo-d,d-c.hi}));}if(c.hard&&err>.05)++s.hard_fail;}return s;}
    void report(const std::string& path,const Stats& before,const Stats& after,bool valid){std::ofstream f(path);if(!f)throw std::runtime_error("Cannot write report");f<<std::setprecision(8)<<"{\n  \"engine\": \"REEC\",\n  \"valid\": "<<(valid?"true":"false")<<",\n  \"heavy_atoms\": "<<nodes.size()<<",\n  \"initial_energy\": "<<before.score<<",\n  \"final_energy\": "<<after.score<<",\n  \"initial_severe_clashes\": "<<before.severe<<",\n  \"severe_clashes\": "<<after.severe<<",\n  \"max_bond_error\": "<<after.bond_max<<",\n  \"max_geometry_distance_error\": "<<after.angle_max<<",\n  \"chirality_errors\": "<<after.chirality<<",\n  \"mandatory_violations\": "<<after.hard_fail<<",\n  \"contacts\": [";for(size_t i=0;i<contacts.size();++i){auto& c=contacts[i];double err=1e100,distance=0;for(auto p:c.alternatives){double d=(xyz[p.first]-xyz[p.second]).norm(),v=std::max({0.,c.lo-d,d-c.hi});if(v<err){err=v;distance=d;}}if(i)f<<",";f<<"\n    {\"res1\": "<<c.r1<<", \"res2\": "<<c.r2<<", \"distance\": "<<distance<<", \"min\": "<<c.lo<<", \"max\": "<<c.hi<<", \"hard\": "<<(c.hard?"true":"false")<<", \"violation\": "<<err<<"}";}f<<"\n  ]\n}\n";if(!f)throw std::runtime_error("Report write failed");}
    void save(const std::string& path,const Stats& s){std::string temp=path+".reec-"+std::to_string(getpid());FILE* f=std::fopen(temp.c_str(),"w");if(!f)throw std::runtime_error("Cannot write PDB: "+path);fprintf(f,"REMARK 250 REEC Restrained Ensemble Energy Conformer\nREMARK 250 ENERGY %.6f MAX_BOND_ERROR %.5f HARD_VIOLATIONS %d\n",s.score,s.bond_max,s.hard_fail);int serial=0;std::vector<int> serials;for(size_t i=0;i<nodes.size();++i){auto& n=nodes[i];n.atom->move(Point(xyz[i].x,xyz[i].y,xyz[i].z));serials.push_back(++serial);std::string atom=n.name.size()<4?" "+n.name:n.name;auto aa=protein.get_residue(n.res);fprintf(f,"ATOM  %5d %-4s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n",serial,atom.c_str(),names[codes.find(aa->get_letter())].c_str(),n.res,xyz[i].x,xyz[i].y,xyz[i].z,Atom::esym_from_Z(n.atom->Z));}fprintf(f,"TER\n");for(auto p:disulfides)fprintf(f,"CONECT%5d%5d\n",serials[id(p.first,"SG")],serials[id(p.second,"SG")]);fprintf(f,"END\n");bool failed=std::ferror(f);if(std::fclose(f))failed=true;if(failed||std::rename(temp.c_str(),path.c_str())){std::remove(temp.c_str());throw std::runtime_error("PDB publication failed");}}
};
void help(){std::cout<<"REEC - Restrained Ensemble Energy Conformer\nUsage: bin/reec --alignment FILE.ali --template [ID=]FILE.pdb [--template ...]\n       --output FILE.pdb [--target ID] [--contacts FILE] [--chain A]\n       [--iterations 3000] [--models 1] [--report FILE.json]\nRun from the Aroma repository root. Exit 2 means validation failed; no PDB is published.\nContacts: disulfide RES RES\n          contact hard|soft RES ATOM RES ATOM MIN MAX WEIGHT\nAtom alternatives use |, for example OD1|OD2. Distances are Angstroms.\n";}
}
#ifndef REEC_NO_MAIN
int main(int argc,char** argv){using namespace reec;try{std::string ali,out,target,report;std::vector<std::string> templates,contact_files;char chain='A';int iterations=3000,models=1;for(int i=1;i<argc;++i){std::string a=argv[i];if(a=="--help"||a=="-h"){help();return 0;}if(i+1>=argc)throw std::runtime_error("Missing option value: "+a);std::string v=argv[++i];if(a=="--alignment")ali=v;else if(a=="--output")out=v;else if(a=="--template")templates.push_back(v);else if(a=="--target")target=v;else if(a=="--contacts")contact_files.push_back(v);else if(a=="--report")report=v;else if(a=="--chain"){if(v.size()!=1)throw std::runtime_error("Chain must be one character");chain=v[0];}else if(a=="--iterations")iterations=integer(v);else if(a=="--models")models=integer(v);else throw std::runtime_error("Unknown option: "+a);}
    if(ali.empty()||out.empty()||templates.empty()){help();return 1;}if(iterations<0||iterations>100000||models<1||models>32)throw std::runtime_error("Invalid iteration or model count");if(report.empty())report=out+".json";if(report==out)throw std::runtime_error("Report and PDB paths must differ");auto entries=alignment(ali);Entry* chosen=nullptr;for(auto& e:entries)if(!e.structure&&(target.empty()||e.id==target)){if(chosen)throw std::runtime_error("Select one target with --target");chosen=&e;}if(!chosen)throw std::runtime_error("Target absent from PIR");std::string seq;for(char c:chosen->seq)if(c!='-')seq+=c;if(seq.size()<2||seq.size()>5000)throw std::runtime_error("Target length outside 2..5000");std::vector<Template> t;for(auto& path:templates){t.push_back(pdb(path,chain));map_template(t.back(),entries,*chosen,seq);}for(size_t i=1;i<t.size();++i)superpose(t[i],t[0]);
    Model m;m.build(seq,t);for(auto& path:contact_files)m.read_contacts(path);Stats initial=m.stats();XYZ start=m.xyz,best;double best_score=1e100;Stats final;bool valid=false;
    for(int candidate=0;candidate<models;++candidate){m.xyz=start;if(candidate)for(size_t i=0;i<m.xyz.size();++i)if(!m.nodes[i].backbone)m.xyz[i]+=V{std::sin(i*1.7+candidate),std::cos(i*.9+candidate),std::sin(i*.4-candidate)}*.15;std::cout<<"REEC candidate "<<candidate+1<<"/"<<models<<"\n";
        for(int stage=0;stage<3;++stage){m.hard_scale=stage==0?10:(stage==1?100:1000);m.anchor_scale=stage==0?1:(stage==1?.3:.1);m.minimize(iterations/3+(stage==2?iterations%3:0));}Stats s=m.stats();bool ok=s.hard_fail==0&&s.chirality==0&&s.bond_max<=.08&&s.angle_max<=.18&&s.severe==0;std::cout<<"Candidate energy "<<s.score<<" bond error "<<s.bond_max<<" geometry error "<<s.angle_max<<" severe clashes "<<s.severe<<" mandatory violations "<<s.hard_fail<<"\n";if((ok&&!valid)||(ok==valid&&s.score<best_score)){best=m.xyz;best_score=s.score;final=s;valid=ok;}}
    m.xyz=best;m.report(report,initial,final,valid);if(!valid){std::cerr<<"REEC validation failed; inspect "<<report<<". Existing output remains intact.\n";return 2;}m.save(out,final);std::cout<<"REEC wrote "<<out<<" and "<<report<<"\n";return 0;
}catch(const std::exception& e){std::cerr<<"REEC: "<<e.what()<<"\n";return 1;}catch(...){std::cerr<<"REEC: molecular library error\n";return 1;}}

#endif
