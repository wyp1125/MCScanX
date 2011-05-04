#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <set>
#include <map>
#include <vector>
#include <algorithm>
using namespace std;

struct Blast_record
{
    string gene1,gene2;
    string mol1,mol2;
    int id1,id2;
//    int pair_id;
    double score;
};

vector<Blast_record> match_list;

struct Gene_feat
{
    string name;
    string mol;
    int mid;
    int geneid;
    bool operator < (const Gene_feat &g) const
    {
        return (mol == g.mol && mid < g.mid) || mol < g.mol;
    }
};
struct geneCmp
{
    bool operator() (const Gene_feat *a, const Gene_feat *b) const
    {
        return (a->mol == b->mol && a->mid < b->mid) || a->mol < b->mol;
    }
};
struct tandem_pair
{
    string gene1,gene2;
    string sp;
    string mol;
    int id;
    bool operator <  (const tandem_pair &g) const
    {
        return (mol == g.mol && id < g.id) || mol < g.mol;
    }
};

struct tandem_array
{
    vector<string> gene;
    string mol;
    string sp;
};

struct tandem_cluster
{
    string anchor1;
    string anchor2;
    int array_id1;
    int array_id2;
};

vector <tandem_pair> alltandempair;
vector <tandem_array> alltandemarray;
vector <tandem_cluster> alltandemcluster;
vector <int> is_show;
map <string,int> tandemgeneid;

typedef set<Gene_feat *, geneCmp> geneSet;

map<string, Gene_feat> gene_map;
geneSet allg;

