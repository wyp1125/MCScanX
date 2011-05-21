#include <iostream>
#include <fstream>
#include <sys/stat.h>

#include "struct.h"

struct New_endpoint
{
    Gene_feat *n;
    int index;
    bool start;
    Gene_feat *e;
    bool operator <  (const New_endpoint &g) const
    {
        return (n->mol == g.n->mol && n->mid < g.n->mid) || n->mol < g.n->mol;
    }
};
/*
struct more_feat
{
int gene_id;
int tandem;
int break_point;
int depth;
};
*/
vector<more_feat> gene_more;
//geneSet allg;
