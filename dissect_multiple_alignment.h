#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

struct Gene_feat
{
    string name;
    string mol;
    int in_blocks;
    int cr_blocks;
    set<string> sp;
    int mid;
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

typedef set<Gene_feat *, geneCmp> geneSet;

map<string, Gene_feat> gene_map;
geneSet allg;

map<string, vector<int> > stat1;
map<string, vector<int> > stat2;
