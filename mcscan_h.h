#ifndef __MCSCAN_H
#define __MCSCAN_H

#include "struct.h"
extern void read_orthomcl(const char *prefix_fn, bool gff_flag=true);
extern void read_gff(const char *prefix_fn);
extern void feed_dag (const string &mol_pair);

extern void msa_main(const char *prefix_fn);

extern void print_align(FILE* fw);

/***** Instantiate all data *****/
map<string, Gene_feat> gene_map;
vector<Blast_record> match_list;
vector<Seg_feat> seg_list;
map<string, int> mol_pairs;
geneSet allg;
map<string, ortho_stat> cmp_sp;

/***** CONSTANTS *****/
int MATCH_SCORE;
int MATCH_SIZE;
int GAP_PENALTY;
int GAP_SIZE;
int OVERLAP_WINDOW;
int UNIT_DIST;
double E_VALUE;
int MAX_GAPS;
string PIVOT;
int EXTENSION_DIST;
int CUTOFF_SCORE;
int IN_SYNTENY;
int e_mode;

#endif
