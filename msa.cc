
#include "msa.h"

static vector <New_endpoint> endpoints;
int max_level;

void fill_allg()
{
Gene_feat *gf1;
map<string, Gene_feat>::iterator it;
for (it=gene_map.begin(); it!=gene_map.end(); it++)
{
gf1 = &(it->second);
allg.insert(gf1);
}
}

void get_endpoints()
{
int n=seg_list.size();
int i,j;
Seg_feat *s;
for(i=0;i<n;i++)
{
s = &seg_list[i];
s->index=i;
New_endpoint ep;

ep.n=s->s1;
ep.index=2*i;
ep.start=true;
ep.e=s->t1;
endpoints.push_back(ep);
ep.n=s->t1;
ep.index=2*i;
ep.start=false;
ep.e=s->s1;
endpoints.push_back(ep);

ep.n=s->s2;
ep.index=2*i+1;
ep.start=true;
ep.e=s->t2;
endpoints.push_back(ep);
ep.n=s->t2;
ep.index=2*i+1;
ep.start=false;
ep.e=s->s2;
endpoints.push_back(ep);
}
sort(endpoints.begin(), endpoints.end());
}

void add_block(Gene_feat* s, Gene_feat* t, int level)
{
int i,j;
i=0;
geneSet::iterator it,it1,it11;
it=allg.find(s);
it1=allg.find(t);
it11=allg.end();
it11--;
for(;(*it)->mid<=(*it1)->mid&&(*it)->mol==(*it1)->mol;it++)
{
i=(*it)->cursor.size();
if(i<level)
{
for(j=i+1;j<level;j++)
{
(*it)->cursor.push_back(0);
}
(*it)->cursor.push_back(1);
}
else
{
(*it)->cursor[level-1]=1;
}
if(it==it11)
{
break;
}
}
}

void add_matchpoints(int seg_index,int level)
{
Seg_feat *s;
int i,j,k;
map<string, Gene_feat>::iterator it5;
i=(int)(seg_index/2);
s=&seg_list[i];
if(seg_index%2==0)
{
for(j=0;j<s->pids.size();j++)
{
it5=gene_map.find(match_list[s->pids[j]].gene1);
if(it5->second.cursor.size()>=level)
it5->second.cursor[level-1]=s->pids[j]+2;
}
}
else
{
for(j=0;j<s->pids.size();j++)
{
it5=gene_map.find(match_list[s->pids[j]].gene2);
if(it5->second.cursor.size()>=level)
it5->second.cursor[level-1]=-(s->pids[j]+2);
}
}
}

void traverse()
{
int i,j,k,lev;
Gene_feat gf4;
char* temp;
map<string, Gene_feat>::iterator it4;
for(i=0;i<endpoints.size();i++)
{
if(endpoints[i].start==1)
{
it4=gene_map.find(endpoints[i].n->name);
gf4=it4->second;
k=gf4.cursor.size();
if(k==0)
{
add_block(endpoints[i].n,endpoints[i].e,1);
add_matchpoints(endpoints[i].index,1);
}

else
{
for(j=0;j<k;j++)
{
if(gf4.cursor[j]==0)
{
lev=j+1;
break;
}
}
if(j==k)
lev=j+1;
add_block(endpoints[i].n,endpoints[i].e,lev);
add_matchpoints(endpoints[i].index,lev);
if(lev>max_level)
max_level=lev;
}
}
else
{
;
}
}
}

void mark_tandem()
{
int i,j,k;
i=0;
geneSet::const_iterator it7=allg.begin();
more_feat mf;
    for (; it7!=allg.end(); it7++)
    {
        (*it7)->gene_id=i;
        mf.depth=0;
        mf.tandem=0;
        //mf.break_point=false;
        for(j=0;j<(*it7)->cursor.size();j++)
        {
        if((*it7)->cursor[j]!=0)
        {
        mf.depth++;
        }
        }
        gene_more.push_back(mf);       
        i++;
    }
    map<string, Gene_feat>::iterator it8,it9;
    for(i=0;i<match_list.size();i++)
    {
    it8=gene_map.find(match_list[i].gene1);
    it9=gene_map.find(match_list[i].gene2);
    if(fabs(it8->second.gene_id-it9->second.gene_id)==1&&it8->second.mol==it9->second.mol)
    {
    gene_more[it8->second.gene_id].tandem=1;
    gene_more[it9->second.gene_id].tandem=1;
    }
    }
}

void print_html()
{
int i,j,k;
string color;
ofstream result;
string prev_mol="";
Gene_feat *n;
char result_dir[200];
geneSet::iterator it6;
i=0;
for(it6=allg.begin();it6!=allg.end();it6++)
{
n=(*it6);
if(n->mol!=prev_mol)
{
if(i>0)
{
result<<"</table></html>";
result.close();
}
//sprintf(result_dir,"%s.html/%s.html",prefix_fn,n->mol.c_str());
sprintf(result_dir,"%s.html",n->mol.c_str());
cout<<result_dir<<endl;
result.open(result_dir,ios::out);
result<<"<html><table cellspacing='0' cellpadding='0' align='left'>";
result<<"<tr align='center'><td># of blocks</td><td>PIVOT</td><td colspan='"<<2*max_level<<"'>Syntenic blocks</td></tr>"<<endl;
prev_mol=n->mol;
//cout<<prev_mol<<endl;
i++;
}
color="'#dddddd'";
if(gene_more[n->gene_id].tandem)
color="'#ee0000'";
//if(gene_more[n->gene_id].break_point)
//color="'#0000ee'";
result<<"<tr align='center'><td>"<<gene_more[n->gene_id].depth<<"</td><td bgcolor="<<color<<">"<<n->name<<"</td>";
for(j=0;j<n->cursor.size();j++)
{
result<<"<td>&nbsp;&nbsp;</td>";
if(n->cursor[j]==0)
{
result<<"<td>&nbsp;</td>";
}
else if(n->cursor[j]==1)
{
result<<"<td>|&nbsp;|</td>";
}
else if(n->cursor[j]>1)
{
result<<"<td bgcolor='#ffff99'>"<<match_list[n->cursor[j]-2].gene2<<"</td>";
}
else
{
result<<"<td bgcolor='#ffff99'>"<<match_list[-n->cursor[j]-2].gene1<<"</td>";
}
}
for(j=n->cursor.size();j<max_level;j++)
{
result<<"<td>&nbsp;</td>";
}
result<<"</tr>"<<endl;
}
result<<"<html><table>"<<endl;
result.close();
}

void print_test()
{
    int j=0;
    geneSet::const_iterator i=allg.begin();
    for (; i!=allg.end(); i++)
    {
        cout<<(*i)->name.c_str()<<" "<<(*i)->cursor.size()<<endl;
        j++;
    if(j>1000)
    break;
    }

}

void msa_main()
{
max_level=1;
fill_allg();
get_endpoints();
traverse();
mark_tandem();
print_html();
}
