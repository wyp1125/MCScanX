#include "detect_syntenic_tandem_arrays.h"

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
    gene_map[gf.name] = gf;
    }
Gene_feat *gf1;
map<string, Gene_feat>::iterator it;
for (it=gene_map.begin(); it!=gene_map.end(); it++)
{
gf1 = &(it->second);
allg.insert(gf1);
}
int i=0;
geneSet::iterator it6;
for(it6=allg.begin();it6!=allg.end();it6++)
{
(*it6)->geneid=i;
i++;
}
}

void read_blast(char* path)
{
    ifstream in(path);
    int i;
    int total_num=0;
    string line,word,geneids,gene1,gene2;
    double evalue;
    map<string, double>blast_map;
    map<string, double>::iterator it;
    cout<<"Reading BLAST file and pre-processing"<<endl;
    while(!in.eof())
    {
    getline(in,line);
    if(line=="")
    break;
    istringstream test(line);
    getline(test,gene1,'\t');
    getline(test,gene2,'\t');
    getline(test,word,'\t');getline(test,word,'\t');getline(test,word,'\t');getline(test,word,'\t');getline(test,word,'\t');getline(test,word,'\t');getline(test,word,'\t');getline(test,word,'\t');
    getline(test,word,'\t');
    istringstream double_iss(word);
    double_iss>>evalue;
    i=gene1.compare(gene2);
    if(i==0)
    { continue; }
    else if(i<0)
    { geneids=gene1+"&"+gene2; }
    else
    { geneids=gene2+"&"+gene1; }
    it = blast_map.find(geneids);
    if(it==blast_map.end())
    { blast_map[geneids]=evalue; }
    else
    {
    if(evalue<it->second)
    { 
    it->second=evalue;
    }
    }
    total_num++;
    }
    in.close();
    double score;
    Blast_record br;
    map<string, Gene_feat>::iterator it1, it2;
    Gene_feat *gf1, *gf2;
    cout<<"Generating BLAST list"<<endl;
    for(it=blast_map.begin();it!=blast_map.end();it++)
    {
    istringstream test(it->first);
    getline(test,gene1,'&');
    getline(test,gene2,'&');
    it1 = gene_map.find(gene1);
    it2 = gene_map.find(gene2);
    if (it1==gene_map.end() || it2==gene_map.end()) continue;
    gf1 = &(it1->second), gf2 = &(it2->second);
    if (gf1->mol.empty() || gf2->mol.empty()) continue;
//////////////////////////////////////////////////////////////////////////////////
if(gf1->geneid<gf2->geneid)
{
    br.gene1=gene1;
    br.gene2=gene2;
    br.mol1=gf1->mol;
    br.mol2=gf2->mol;
    br.id1=gf1->geneid;
    br.id2=gf2->geneid;
}
else
{
    br.gene1=gene2;
    br.gene2=gene1;
    br.mol1=gf2->mol;
    br.mol2=gf1->mol;
    br.id1=gf2->geneid;
    br.id2=gf1->geneid;
}
//////////////////////////////////////////////////////////////////////////////////
    br.score = it->second;
    match_list.push_back(br);
    }

   int selected_num = match_list.size();
   cout<<selected_num<<" matches imported ("<<total_num-selected_num<<" discarded)"<<endl;
}

void compute_tandem()
{
cout<<"Detecting tandem arrays..."<<endl;
int i,j,k;
tandem_pair temp;
for(i=0;i<match_list.size();i++)
{
if(match_list[i].mol1==match_list[i].mol2&&(match_list[i].id1-match_list[i].id2)==-1)
{
temp.gene1=match_list[i].gene1;
temp.gene2=match_list[i].gene2;
temp.mol=match_list[i].mol1;
temp.sp=match_list[i].mol1.substr(0,2);
temp.id=match_list[i].id1;
alltandempair.push_back(temp);
}
}
sort(alltandempair.begin(),alltandempair.end());
tandem_array temp1;
temp1.gene.push_back(alltandempair[0].gene1);
temp1.gene.push_back(alltandempair[0].gene2);
temp1.sp=alltandempair[0].sp;
temp1.mol=alltandempair[0].mol;
j=alltandempair[0].id;
for(i=1;i<alltandempair.size();i++)
{
if(alltandempair[i].id==j+1)
{
temp1.gene.push_back(alltandempair[i].gene2);
j++;
}
else
{
alltandemarray.push_back(temp1);
temp1.gene.clear();
temp1.gene.push_back(alltandempair[i].gene1);
temp1.gene.push_back(alltandempair[i].gene2);
temp1.sp=alltandempair[i].sp;
temp1.mol=alltandempair[i].mol;
j=alltandempair[i].id;
}
}
for(i=0;i<alltandemarray.size();i++)
{
for(j=0;j<alltandemarray[i].gene.size();j++)
{
tandemgeneid[alltandemarray[i].gene[j]]=i;
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
string gene1,gene2;
ifstream in(path);
string line, word;
size_t found;
tandem_cluster temp;
map<string,int>::iterator it21,it22;
while(!in.eof())
    {
    getline(in,line);
    if(line=="")
    continue;
    found=line.find("#");
    if (found!=string::npos)
    {
    continue;
    }
    istringstream test(line);
    getline(test,word,'\t');
    getline(test,gene1,'\t');
    getline(test,gene2,'\t');
    it21=tandemgeneid.find(gene1);
    it22=tandemgeneid.find(gene2);
    if(it21!=tandemgeneid.end()||it22!=tandemgeneid.end())
    {
    temp.anchor1=gene1;
    temp.anchor2=gene2;
    temp.array_id1=-1;
    temp.array_id2=-1;
    if(it21!=tandemgeneid.end())
    temp.array_id1=it21->second;
    if(it22!=tandemgeneid.end())
    temp.array_id2=it22->second;
    alltandemcluster.push_back(temp);
    }
    }
}

void print_file(char* path)
{
cout<<"Outputting results..."<<endl;
int i,j;
ofstream result;
result.open(path,ios::out);
result<<"Anchor1\tTandem_array1\tAnchor2\tTandem_array2\n";
for(i=0;i<alltandemcluster.size();i++)
{
result<<alltandemcluster[i].anchor1<<"\t";
if(alltandemcluster[i].array_id1<0)
{
result<<alltandemcluster[i].anchor1<<"\t";
}
else
{
for(j=0;j<alltandemarray[alltandemcluster[i].array_id1].gene.size()-1;j++)
{
result<<alltandemarray[alltandemcluster[i].array_id1].gene[j]<<",";
}
result<<alltandemarray[alltandemcluster[i].array_id1].gene[j]<<"\t";
}
result<<alltandemcluster[i].anchor2<<"\t";
if(alltandemcluster[i].array_id2<0)
{
result<<alltandemcluster[i].anchor2<<endl;
}
else
{
for(j=0;j<alltandemarray[alltandemcluster[i].array_id2].gene.size()-1;j++)
{
result<<alltandemarray[alltandemcluster[i].array_id2].gene[j]<<",";
}
result<<alltandemarray[alltandemcluster[i].array_id2].gene[j]<<endl;
}
}
cout<<"Done!"<<endl;
}

int main(int argc, char *argv[])
{
if (argc < 9)
{
cout<<"Usage: ./detect_syntenic_tandem_arrays -g gff_file -b blast_file -s synteny_file -o output_file"<<endl;
exit(1);
}
char gpath[200], bpath[200], spath[200], opath[200];
int c,e_g,e_b,e_s,e_o;
e_g=e_b=e_s=e_o=1;
while((c = getopt(argc, argv, "g:b:s:o:")) != -1)
{
switch(c)
{
case 'g':
sprintf(gpath,"%s",optarg);
e_g=0;
break;
case 'b':
sprintf(bpath,"%s",optarg);
e_b=0;
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
if (optopt!='g' || optopt!='b' || optopt!='s' ||optopt!='o')
{
cout<<"Usage: ./detect_syntenic_tandem_arrays -g gff_file -b blast_file -s synteny_file -o output_file"<<endl;
}
break;
}
}
if(e_g+e_b+e_s+e_o>0)
{
cout<<"Usage: ./detect_syntenic_tandem_arrays -g gff_file -b blast_file -s synteny_file -o output_file"<<endl;
exit(1);
}
read_gff(gpath);
read_blast(bpath);
compute_tandem();
read_synteny(spath);
print_file(opath);
return 1;
}
