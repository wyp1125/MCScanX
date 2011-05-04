/*
 * Author: Haibao Tang <bao@uga.edu> May 10, 2007
 *
 * Data input module for mcscan, contains several procedures
 * Read blast output file, formatted by -m8, is a bunch of hits
 * Read MCL cluster, which is retrieved by clustering the above blast file
 * Read GFF file, which includes the chromosome, position information

 *Modified by Yupeng Wang, Mar 31, 2011
 *The part of reading MCL was removed
 */


#include "read_data.h"

// incremental sorting y coord
static bool cmp_y (const Score_t& t1, const Score_t& t2)
{
    return t1.y < t2.y ||
           (t1.y == t2.y && t1.x < t2.x);
}

// incremental sorting e-value
static bool cmp_ev (const Score_t& t1, const Score_t& t2)
{
    return t1.score < t2.score;
}

// filter the blast -m8 output by the following threshold:
// lexically sorted, gene #1 < gene #2
// non-self blast match
// both be present in the mcl output file and in the same group
void read_blast(const char *prefix_fn, bool gff_flag=true)
{

    char fn[LABEL_LEN], g1[LABEL_LEN], g2[LABEL_LEN];
    sprintf(fn,"%s.blast",prefix_fn);
    ifstream in(fn);
    int i;
    int total_num=0;
    string line,word,geneids,gene1,gene2;
    double evalue;
    map<string, double>blast_map;
    map<string, double>::iterator it;
    cout<<"Reading BLAST file and pre-processing"<<endl;
    while (!in.eof())
    {
        getline(in,line);
        if (line=="")
            break;
        istringstream test(line);
        getline(test,gene1,'\t');
//    gene1=word;
        getline(test,gene2,'\t');
//    gene2=word;
        getline(test,word,'\t');
        getline(test,word,'\t');
        getline(test,word,'\t');
        getline(test,word,'\t');
        getline(test,word,'\t');
        getline(test,word,'\t');
        getline(test,word,'\t');
        getline(test,word,'\t');
        getline(test,word,'\t');
        istringstream double_iss(word);
        double_iss>>evalue;
        i=gene1.compare(gene2);
        if (i==0)
        {
            continue;
        }
        else if (i<0)
        {
            geneids=gene1+"&"+gene2;
        }
        else
        {
            geneids=gene2+"&"+gene1;
        }
        it = blast_map.find(geneids);
        if (it==blast_map.end())
        {
            blast_map[geneids]=evalue;
        }
        else
        {
            if (evalue<it->second)
            {
                it->second=evalue;
            }
        }
        total_num++;
    }
    in.close();

    double score;
    Blast_record br;
    int pair_id = 0;
    //int total_num = 0;
    map<string, Gene_feat>::iterator it1, it2;
    Gene_feat *gf1, *gf2;
    cout<<"Generating BLAST list"<<endl;
    for (it=blast_map.begin();it!=blast_map.end();it++)
    {
        istringstream test(it->first);
        getline(test,gene1,'&');
        getline(test,gene2,'&');
        //total_num++;
        it1 = gene_map.find(gene1);
        it2 = gene_map.find(gene2);
        if (it1==gene_map.end() || it2==gene_map.end()) continue;
        gf1 = &(it1->second), gf2 = &(it2->second);
        if (gf1->mol.empty() || gf2->mol.empty()) continue;
        if (IN_SYNTENY==1 && gf1->mol.substr(0,2)!=gf2->mol.substr(0,2)) continue;
        if (IN_SYNTENY==2 && gf1->mol.substr(0,2)==gf2->mol.substr(0,2)) continue;
/////////////bug here/////////////////////////////////////////////////////////////
        i=gf1->mol.compare(gf2->mol);
//////////////////////////////////////////////////////////////////////////////////
        if (i<0)
        {
            br.gene1=gene1;
            br.gene2=gene2;
            br.mol_pair = gf1->mol+"&"+gf2->mol;
        }
        else if (i==0)
        {
            if (gf1->mid<=gf2->mid)
            {
                br.gene1=gene1;
                br.gene2=gene2;
            }
            else
            {
                br.gene1=gene2;
                br.gene2=gene1;
            }
            br.mol_pair = gf1->mol+"&"+gf2->mol;
        }
        else
        {
            br.gene1=gene2;
            br.gene2=gene1;
            br.mol_pair = gf2->mol+"&"+gf1->mol;
        }
//////////////////////////////////////////////////////////////////////////////////
        mol_pairs[br.mol_pair]++;

        br.pair_id = pair_id++;
        br.score = it->second;
        match_list.push_back(br);

    }

    int selected_num = match_list.size();
    progress("%d matches imported (%d discarded)",
             selected_num, total_num - selected_num);
}

void read_gff(const char *prefix_fn)
{
    char fn[LABEL_LEN], gn[LABEL_LEN], mol[LABEL_LEN];
    int end5, end3;
    Gene_feat gf;
    //master_feat mf;
    sprintf(fn, "%s.gff", prefix_fn);
    FILE *fp = mustOpen(fn, "r");

    while (fscanf(fp, "%s%s%d%d",
                  &mol[0], &gn[0], &end5, &end3) == 4)
    {
        gf.mol = string(mol);
        gf.name = string(gn);
        gf.mid = end5;
        //gf.depth = 0;
        gene_map[gf.name] = gf;
        //mf.n=gf;
        //mf.depth=0;
        //master_map[gf.name]=mf;
    }

    fclose(fp);
}

static void filter_matches_x ()
{
    // match_bin is a list of records that are potentially repetitive
    vector<Score_t> match_bin, score_cpy;
    vector<Score_t>::const_iterator it, prev_rec;

    sort(score.begin(), score.end());
    prev_rec = it = score.begin();
    it++;
    match_bin.push_back(*(prev_rec));
    for (; it != score.end(); it++)
    {
        // scan whether it has a linking window with previous one
        if ((prev_rec->x != it->x) ||
                (it->y - prev_rec->y) > OVERLAP_WINDOW)
        {
            // record last match_bin, take only least e-value
            score_cpy.push_back(*min_element(match_bin.begin(),
                                             match_bin.end(), cmp_ev));
            // start a new match_bin
            match_bin.clear();
        }
        match_bin.push_back(*it);
        prev_rec = it;
    }
    // don't forget the last match_bin
    score_cpy.push_back(*min_element(match_bin.begin(),
                                     match_bin.end(), cmp_ev));
    match_bin.clear();

    // copy into score
    score.clear();
    score = score_cpy;
    score_cpy.clear();
}

static void filter_matches_y ()
{
    // match_bin is a list of records that are potentially repetitive
    vector<Score_t> match_bin, score_cpy;
    vector<Score_t>::const_iterator it, prev_rec;

    sort(score.begin(), score.end(), cmp_y);
    prev_rec = it = score.begin();
    it++;
    match_bin.push_back(*(prev_rec));
    for (; it != score.end(); it++)
    {
        // scan whether it has a linking window with previous one
        if ((prev_rec->y != it->y) ||
                (it->x - prev_rec->x) > OVERLAP_WINDOW)
        {
            // record last match_bin, take only least e-value
            score_cpy.push_back(*min_element(match_bin.begin(),
                                             match_bin.end(), cmp_ev));
            // start a new match_bin
            match_bin.clear();
        }
        match_bin.push_back(*it);
        prev_rec = it;
    }
    // don't forget the last match_bin
    score_cpy.push_back(*min_element(match_bin.begin(),
                                     match_bin.end(), cmp_ev));
    match_bin.clear();

    // copy into score
    score.clear();
    score = score_cpy;
    score_cpy.clear();
}

// feed into dagchainer
void feed_dag(const string &mol_pair)
{
    // two additional filters will be applied here
    // best hsp (least e-value)
    // non-repetitive in a window of 50kb region
    vector<Blast_record>::const_iterator it;
    Score_t cur_score;

    for (it = match_list.begin(); it < match_list.end(); it++)
    {
        if (it->mol_pair != mol_pair) continue;

        cur_score.pairID = it->pair_id;
        cur_score.x = gene_map[it->gene1].mid;
        cur_score.y = gene_map[it->gene2].mid;
        cur_score.gene1=it->gene1;
        cur_score.gene2=it->gene2;
        cur_score.score = MATCH_SCORE;

        score.push_back(cur_score);
    }

    // sort by both axis and remove redundant matches within
    // a given window length (default 50kb)
    filter_matches_x();
    filter_matches_y();

    dag_main(score, mol_pair);
}

