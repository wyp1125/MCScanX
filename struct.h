#ifndef __STRUCT_H
#define __STRUCT_H

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <cstdio>
#include <cstring>
#include <cstdarg>
#include <cerrno>
#include <climits>
#include <cmath>
#include <ctime>
#include <sys/time.h>
#include <getopt.h>
#include <string>
#include <vector>
#include <list>
#include <queue>
#include <set>
#include <map>
#include <algorithm>

#ifndef M_LN2
#define M_LN2 0.69314718055994530942
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

#define sameString(a, b) (strcmp((a), (b))==0)
#define MAX(a, b) ((a)>(b)?(a):(b))
#define MIN(a, b) ((a)<(b)?(a):(b))
#define LABEL_LEN 200

#ifndef __GNUC__
#define __attribute__(x)
#endif

struct Blast_record
{
    string gene1, gene2;
    string mol_pair;
    int pair_id;
    int node;
    double score;
};

struct more_feat
{
//int gene_id;
    int tandem;
//bool break_point;
    int depth;
};

struct Gene_feat
{
    vector<int> cursor;
    string name;
    string mol;
    int mid;
    int gene_id;
//   more_feat *m;
    bool operator < (const Gene_feat &g) const
    {
        return (mol == g.mol && mid < g.mid) || mol < g.mol || (mol == g.mol && mid == g.mid && name.compare(g.name)<0);
    }
};

struct geneCmp
{
    bool operator() (const Gene_feat *a, const Gene_feat *b) const
    {
        return (a->mol == b->mol && a->mid < b->mid) || a->mol < b->mol || (a->mol == b->mol && a->mid == b->mid && a->name.compare(b->name)<0);
    }
};

typedef set<Gene_feat *, geneCmp> geneSet;

struct Seg_feat
{
    vector<int> pids;
    int index;
    Gene_feat *s1, *t1, *s2, *t2;
    double score, e_value;
    string mol_pair;
    bool sameStrand;
};

struct Cell_t
{
    float raw;
int score :
    30;
unsigned from :
    2;
};

struct Score_t
{
    int pairID;  // identifier of match pair
    int x, y;  // x,y coordinates
    float score;
    string gene1;
    string gene2;
    bool operator< (const Score_t & node) const
    {
        return  (x < node.x || (x == node.x && y < node.y));
    }
};

struct Path_t
{
    float score;
    int rc;  // sum of row and column of last entry
    int sub;
};
//////for MCScanX_orthomcl////////////////
struct ortho_stat
{
    int all_num;
    int syn_num;
};
extern map<string, ortho_stat> cmp_sp;
/////////////////////////////////////////

extern map<string, Gene_feat> gene_map;
//extern map<string, double>blast_map;
extern vector<Blast_record> match_list;
extern vector<Seg_feat> seg_list;
extern map<string, int> mol_pairs;
extern geneSet allg;
//extern vector<more_feat> gene_more;
//extern map<string, geneSet > chr_map;

/***** CONSTANTS *****/
// match bonus
extern int MATCH_SCORE;
extern int MATCH_SIZE;
// gap extension penalty
extern int GAP_PENALTY;
// length for a single gap in basepairs
extern int GAP_SIZE;
// The filter window for linking locally repetitive hits
extern int OVERLAP_WINDOW;
// significance cut-off
extern double E_VALUE;
// maximum gaps allowed
extern int MAX_GAPS;
// reference genome
//extern string PIVOT;
// intergenic distance
//extern int UNIT_DIST;
// segment extension limit
//extern int EXTENSION_DIST;
// alignment significance score
extern int CUTOFF_SCORE;
extern int IN_SYNTENY;
extern int N_PROXIMAL;
extern int e_mode;
// direction in the 2d dynamic matrix
enum { DIAG, UP, LEFT, DEL };

/***** Helper functions (Some from James Kent library) *****/
void progress(const char *format, ...)
/* Print progress message */
__attribute__((format(printf, 1, 2)));

void err(const char *format, ...)
/* Print error message but do not exit */
__attribute__((format(printf, 1, 2)));

void warn(const char *format, ...)
/* Print error message but do not exit */
__attribute__((format(printf, 1, 2)));

void errAbort(const char *format, ...)
/* Print error message to stderr and exit */
__attribute__((noreturn, format(printf, 1, 2)));

long clock1000();
/* A millisecond clock. */

void uglyTime(const char *label, ...)
/* Print label and how long it's been since last call.  Call with
 * a NULL label to initialize. */
__attribute__((format(printf, 1, 2)));

FILE *mustOpen(const char *fileName, const char *mode);
/* Open a file or die */

#endif
