#include "cls.h"
/*
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
*/
void cmpt_duptype()
{
    int i,j;
    i=0;
    geneSet::const_iterator it7=allg.begin();
    more_feat mf;
    for (; it7!=allg.end(); it7++)
    {
        (*it7)->gene_id=i;
        mf.tandem=0;
        gene_more.push_back(mf);
        i++;
    }

    map<string, Gene_feat>::iterator it8,it9;
    for (i=0;i<match_list.size();i++)
    {
        it8=gene_map.find(match_list[i].gene1);
        it9=gene_map.find(match_list[i].gene2);
        gene_more[it8->second.gene_id].tandem=1;
        gene_more[it9->second.gene_id].tandem=1;
        if (fabs(it8->second.gene_id-it9->second.gene_id)<N_PROXIMAL&&it8->second.mol==it9->second.mol)
        {
            gene_more[it8->second.gene_id].tandem=2;
            gene_more[it9->second.gene_id].tandem=2;
            if (fabs(it8->second.gene_id-it9->second.gene_id)==1)
            {
                gene_more[it8->second.gene_id].tandem=3;
                gene_more[it9->second.gene_id].tandem=3;
            }
        }
    }

    int n=seg_list.size();
    Seg_feat *s;
    for (i=0;i<n;i++)
    {
        s = &seg_list[i];
        for (j=0;j<s->pids.size();j++)
        {
            it8=gene_map.find(match_list[s->pids[j]].gene1);
            it9=gene_map.find(match_list[s->pids[j]].gene2);
            gene_more[it8->second.gene_id].tandem=4;
            gene_more[it9->second.gene_id].tandem=4;
        }
    }
}

void print_cls(char* prefix_fn)
{
    int num[5]={0,0,0,0,0};
    char stat_fn[LABEL_LEN];
    sprintf(stat_fn, "%s.gene_type", prefix_fn);
    ofstream result;
    result.open(stat_fn,ios::out);
    geneSet::const_iterator it6=allg.begin();
    Gene_feat *n;
    for (it6=allg.begin();it6!=allg.end();it6++)
    {
        n=(*it6);
        result<<n->name<<"\t"<<gene_more[n->gene_id].tandem<<endl;
        num[gene_more[n->gene_id].tandem]++;
    }
    result.close();
    cout<<"Type of dup\tCode\tNumber"<<endl;
    cout<<"Singleton\t0\t"<<num[0]<<endl;
    cout<<"Dispersed\t1\t"<<num[1]<<endl;
    cout<<"Proximal\t2\t"<<num[2]<<endl;
    cout<<"Tandem\t3\t"<<num[3]<<endl;
    cout<<"WGD\t4\t"<<num[4]<<endl;
}

void cls_main(char* prefix_fn)
{
 //   fill_allg();
    cmpt_duptype();
    print_cls(prefix_fn);
}
