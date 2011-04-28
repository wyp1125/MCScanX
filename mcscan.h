#ifndef __MCSCAN_H
#define __MCSCAN_H

#include "struct.h"
#include <sys/stat.h>
// read_data
extern void read_blast(const char *prefix_fn, bool gff_flag=true);
//extern void read_mcl(const char *prefix_fn);
extern void read_gff(const char *prefix_fn);
//extern void feed_pog();
extern void feed_dag (const string &mol_pair);
//extern void read_cfg();

// pog
//extern void POG_main(FILE *fw);
extern void msa_main();

// out_utils
extern void print_align(FILE* fw);
//extern void print_align_mcl(FILE *fw);

/***** Instantiate all data *****/
map<string, Gene_feat> gene_map;
//map<string, master_feat> master_map;
//map<string, double>blast_map;
vector<Blast_record> match_list;
vector<Seg_feat> seg_list;
map<string, int> mol_pairs;
//map<string, geneSet > chr_map;

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

#endif
