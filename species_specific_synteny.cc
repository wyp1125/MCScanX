#include "species_specific_synteny.h"

void read_gff(char* path)
{
if(!ifstream(path))
{
cout<<"Cannot read gff_file!"<<endl;
exit(1);
}
ifstream in(path);
string line, word;
Gene_feat gf;
while(!in.eof())
    {
    getline(in,line);
    if(line=="")
    break;
    istringstream test(line);
    getline(test,gf.mol,'\t');
    getline(test,gf.name,'\t');
    getline(test,word,'\t');
    istringstream int_iss(word);
    int_iss>>gf.mid;
    gf.in_blocks=0;
    gf.cr_blocks=0;
    gene_map[gf.name] = gf;
    }
Gene_feat *gf1;
map<string, Gene_feat>::iterator it;
for (it=gene_map.begin(); it!=gene_map.end(); it++)
{
gf1 = &(it->second);
allg.insert(gf1);
}
}
void add_block(Gene_feat* s, Gene_feat* t, int type, string sp)
{
geneSet::iterator it,it1,it11;
it=allg.find(s);
it1=allg.find(t);
it11=allg.end();
it11--;
for(;(*it)->mid<=(*it1)->mid&&(*it)->mol==(*it1)->mol;it++)
{
if(type==0)
{
(*it)->in_blocks++;
}
else
{
(*it)->cr_blocks++;
(*it)->sp.insert(sp);
}
if(it==it11)
{
break;
}
}
}
void read_synteny(char* path)
{
if(!ifstream(path))
{
cout<<"Cannot read synteny_file!"<<endl;
exit(1);
}
int i=0;
int j,strand;
ifstream in(path);
string line, word;
string gene1,gene2,s1,t1,s2,t2;
Gene_feat* gfs1,*gft1,*gfs2,*gft2;
map<string, Gene_feat>::iterator its1,itt1,its2,itt2;
size_t found;
while(!in.eof())
    {
    getline(in,line);
    if(line=="")
    continue;
    found=line.find("#");
    if (found!=string::npos)
    {
    found=line.find("Alignment");
    if (found!=string::npos)
    {
    j=0;
    strand=1;
    found=line.find("minus");
    if (found!=string::npos)
    {
    strand=0;
    }
//////////////////////////////////////////////////////////////////////////////
    if(i>0)
    {
    its1 = gene_map.find(s1);
    itt1 = gene_map.find(t1);
    its2 = gene_map.find(s2);
    itt2 = gene_map.find(t2);
    gfs1 = &(its1->second);
    gft1 = &(itt1->second);
    gfs2 = &(its2->second);
    gft2 = &(itt2->second);
    int type=0;
    string sp1,sp2;
    sp1=gfs1->mol.substr(0,2);
    sp2=gfs2->mol.substr(0,2);
    if(sp1!=sp2)
    type=1;
    add_block(gfs1,gft1,type,sp2);
    add_block(gfs2,gft2,type,sp1);
    //cout<<i-1<<" "<<strand<<" "<<s1<<" "<<s2<<" "<<t1<<" "<<t2<<endl;
    }
//////////////////////////////////////////////////////////////////////////////
    i++;
    }
    continue;
    }
    istringstream test(line);
    getline(test,word,'\t');
    getline(test,gene1,'\t');
    getline(test,gene2,'\t');
    if(j==0)
    {
    s1=gene1;
    if(strand==1)
    s2=gene2;
    else
    t2=gene2;
    }
    else
    {
    t1=gene1;
    if(strand==1)
    t2=gene2;
    else
    s2=gene2;
    }
    j++;
    }
//////////////////////////////////////////////////////////////////////////////
    //cout<<i-1<<" "<<strand<<" "<<s1<<" "<<s2<<" "<<t1<<" "<<t2<<endl;
    its1 = gene_map.find(s1);
    itt1 = gene_map.find(t1);
    its2 = gene_map.find(s2);
    itt2 = gene_map.find(t2);
    gfs1 = &(its1->second);
    gft1 = &(itt1->second);
    gfs2 = &(its2->second);
    gft2 = &(itt2->second);
    int type=0;
    string sp1,sp2;
    sp1=gfs1->mol.substr(0,2);
    sp2=gfs2->mol.substr(0,2);
    if(sp1!=sp2)
    type=1;
    add_block(gfs1,gft1,type,sp2);
    add_block(gfs2,gft2,type,sp1);
//////////////////////////////////////////////////////////////////////////////
}

void print_file(char* path)
{
ofstream result;
result.open(path,ios::out);
Gene_feat *n;
geneSet::iterator it6;
for(it6=allg.begin();it6!=allg.end();it6++)
{
n=(*it6);
result<<n->mol<<"\t"<<n->name<<"\t"<<n->in_blocks<<"\t"<<n->cr_blocks<<"\t"<<n->sp.size()<<"\t";
if(n->sp.size()>0)
{
set<string>::iterator it7;
for(it7=n->sp.begin();it7!=n->sp.end();it7++)
result<<*it7<<",";
}
result<<endl;
}
}

int main(int argc, char *argv[])
{
//cout<<argc<<endl;
if (argc < 7)
{
cout<<"Usage: ./species_specific -g gff_file -s synteny_file -o output_file"<<endl;
exit(1);
}
char gpath[200], spath[200], opath[200];
int c,e_g,e_s,e_o;
e_g=e_s=e_o=1;
while((c = getopt(argc, argv, "g:s:o:")) != -1)
{
switch(c)
{
case 'g':
sprintf(gpath,"%s",optarg);
e_g=0;
break;
case 's':
sprintf(spath,"%s",optarg);
e_s=0;
break;
case 'o':
sprintf(opath,"%s",optarg);
e_o=0;
break;
case '?':
if (optopt!='g' || optopt!='s' || optopt!='o')
{
cout<<"Usage: ./species_specific -g gff_file -s synteny_file -o output_file"<<endl;
}
break;
}
}
if(e_g+e_s+e_o>0)
{
cout<<"Usage: ./species_specific -g gff_file -s synteny_file -o output_file"<<endl;
exit(1);
}
read_gff(gpath);
read_synteny(spath);
print_file(opath);
return 1;
}
