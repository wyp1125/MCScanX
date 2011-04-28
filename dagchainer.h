#ifndef __DAGCHAINER_H
#define __DAGCHAINER_H

#include  "struct.h"

void dag_main(vector<Score_t> &score, const string &mol_pair);

// segment id, plus one when found a pairwise alignment
int ali_ct = 0;
int Best_g = -1, Best_i, Best_j;
int Max_Y;

/* calculation procedure in permutation.cc */
extern double ln_perm(int n, int r);
extern double ln_comb(int n, int k);

#endif
