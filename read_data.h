#ifndef __READ_DATA_H
#define __READ_DATA_H

#include "struct.h"
#include <iostream>
#include <fstream>
#include <sstream>


void read_blast(const string &prefix_fn);
void read_mcl(const string &prefix_fn);
void read_gff(const string &prefix_fn, int gff_flag=1);

void feed_dag(const string &mol_pair);

// dagchainer
extern void dag_main(vector<Score_t>& score, const string &mol_pair);

// use gene name to search its node, mol, mid
vector<Score_t> score;

#endif
