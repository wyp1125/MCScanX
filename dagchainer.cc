/*
 * Author: A. L. Delcher and modified by bhaas.
 *
 * File: dagchainer.cpp
 * Last Modified: 7 November 2003
 * Second modification: Haibao Tang <bao@uga.edu> May 10, 2007
 *
 * Do DP on dag of matches to get chains of matches
 * Changes made in the I/O section, the core chaining algorithm remains intact
 * 
 * Modified by Yupeng Wang <wyp1125@uga.edu> March 31, 2011
 * GAP_PENATY can be set by users
 */

#include "dagchainer.h"

// check whether an alignment overlap (tandem alignment)
static bool check_overlap(vector<int>& xx, vector<int>& yy)
{
    int xmin = *min_element(xx.begin(), xx.end());
    int xmax = *max_element(xx.begin(), xx.end());
    int ymin = *min_element(yy.begin(), yy.end());
    int ymax = *max_element(yy.begin(), yy.end());

    return xmin <= ymax && ymin <= xmax;
}

static void retrieve_pos(int pid, int *pos1, int *pos2)
/* returns pos1, pos2 for the blast pair */
{
    Blast_record *match_rec = &match_list[pid];

    *pos1 = gene_map[match_rec->gene1].mid;
    *pos2 = gene_map[match_rec->gene2].mid;
}

static bool is_significant(Seg_feat *sf, vector<Score_t>& score)
/* test if a syntenic block is significant, see description in permutation.cc */
{
    /* see formula in permutation.cc, unknowns are m, N, L1, L2, l_1i l_2i*/
    int s1_a, s1_b, s2_a, s2_b, m, N=0, L1, L2, i;
    double l1, l2, summation=0;

    /* get the start and stop coordinates on each syntenic segment */
    s1_a = sf->s1->mid, s1_b = sf->t1->mid;
    s2_a = sf->s2->mid, s2_b = sf->t2->mid;

    /* calculate m, number of anchor points */
    m = sf->pids.size();

    /* calculate N, number of matches in the defined region*/
    vector<Score_t>::const_iterator it;
    for (it=score.begin(); it!=score.end(); it++)
    {
        if (it->x >=s1_a && it->x <=s1_b && it->y >=s2_a && it->y <=s2_b)
            N++;
    }

    /* calculate l1, l2, distance between successive anchor points */
    int l1_pos1, l1_pos2, l2_pos1, l2_pos2;
    retrieve_pos(sf->pids[0], &l1_pos1, &l2_pos1);
    for (i=1; i<m; i++)
    {
        retrieve_pos(sf->pids[i], &l1_pos2, &l2_pos2);
        l1 = fabs(l1_pos2 - l1_pos1);
        l2 = fabs(l2_pos2 - l2_pos1);
        l1_pos1 = l1_pos2;
        l2_pos1 = l2_pos2;

        summation += log(l1)+log(l2);
    }

    /* calculate L1, L2, respective length of the matching region */
    L1 = s1_b - s1_a, L2 = s2_b - s2_a;

    /* this is the formula */
    sf->e_value = exp(M_LN2 + ln_perm(N, m) + \
                      summation - (m-1)*(log(L1)+log(L2)));

    return sf->e_value < E_VALUE;
}

static bool Descending_Score(const Path_t &a, const Path_t &b)
{
    return a.score > b.score ||
           (a.score == b.score && a.rc > b.rc);
}

// whether a mol_pair is self comparison, e.g. "at1&at1"
static bool check_self (const string &s)
{
    int pos = s.find('&');
    return s.substr(0, pos) == s.substr(pos+1);
}

static void print_chains(vector<Score_t>& score, const string &mol_pair)
/* Find and output highest scoring chains in score treating it as a DAG*/
{
    vector<float> path_score;
    vector<int> from, ans;
    vector<Path_t> high;
    vector<int> xx, yy;
    Path_t  p;
    bool done;
    int i, j, m, n, s, pid, num_gaps;
    int del_x, del_y;
    double x;
    bool is_self = check_self(mol_pair);

    sort(score.begin(), score.end());

    do
    {
        done = true;
        n = score.size();
        path_score.resize(n);
        from.resize(n);
        for (i=0; i<n; i++)
        {
            path_score[i] = score[i].score;
            from[i] = -1;
        }
        for (j=1; j<n; j++)
        {
            for (i=j-1; i>=0; i--)
            {
                del_x = score[j].x - score[i].x - 1;
                del_y = score[j].y - score[i].y - 1;

                if  (del_x >= 0 && del_y >= 0)
                {
                    // if (del_x > EXTENSION_DIST && del_y > EXTENSION_DIST)
                    //     break;
                    // if (del_x > EXTENSION_DIST || del_y > EXTENSION_DIST)
                    //     continue;
                    if (del_x > MAX_GAPS)
                        break;
                    if (del_y > MAX_GAPS)
                        continue;
                    //num_gaps = MAX(del_x, del_y)/UNIT_DIST;
                    num_gaps = MAX(del_x, del_y);
                    x = path_score[i] + score[j].score;
                    /* gap penalty */
                    if (num_gaps > 0)
                        x +=num_gaps*GAP_PENALTY;
                    if  (x > path_score[j])
                    {
                        path_score[j] = x;
                        from[j] = i;
                    }
                }
            }
        }

        high.clear();
        for (i=0; i<n; i++)
        {
            if (path_score[i] >= CUTOFF_SCORE)
            {
                p.score = path_score[i];
                p.sub = i;
                p.rc = score[i].x + score[i].y;
                high.push_back(p);
            }
        }

        sort (high.begin(), high.end(), Descending_Score);

        m = high.size();
        for  (i=0; i<m; i++)
        {
            if  (from[high[i].sub] != -2)
            {
                ans.clear();
                for  (j=high[i].sub; from[j]>=0; j=from[j])
                {
                    ans.push_back(j);
                }
                ans.push_back(j);
                if (from[j] == -2)
                {
                    done = false;
                    break;
                }
                else
                {
                    reverse(ans.begin(), ans.end());
                    s = ans.size();
                    if (is_self)
                    {
                        for  (j=0; j<s; j++)
                        {
                            from[ans[j]] = -2;
                            xx.push_back(score[ans[j]].x);
                            yy.push_back(score[ans[j]].y);
                        }
                    }

                    Seg_feat sf;
                    Blast_record *br;
                    if (!(is_self && check_overlap(xx, yy)))
                    {
                        sf.score = path_score[high[i].sub];
                        for (j=0; j<s; j++)
                        {
                            from[ans[j]] = -2;

                            pid = score[ans[j]].pairID;
                            br = &match_list[pid];
                            sf.pids.push_back(pid);
                        }
                        /* start and stop positions for two sub-segments */
                        br = &match_list[sf.pids.front()];
                        sf.s1 = &gene_map[br->gene1];
                        sf.s2 = &gene_map[br->gene2];
                        br = &match_list[sf.pids.back()];
                        sf.t1 = &gene_map[br->gene1];
                        sf.t2 = &gene_map[br->gene2];

                        /* determine the orientation of the alignment */
                        sf.sameStrand = *(sf.s2) < *(sf.t2);
                        if (!sf.sameStrand) swap(sf.s2, sf.t2);

                        sf.mol_pair = mol_pair;

                        /* significance testing */
                        if (is_significant(&sf, score))
                            seg_list.push_back(sf);
                    }
                    xx.clear(), yy.clear();
                }
            }
        }

        if (!done)
        {
            for (i=j=0; i<n; i++)
            {
                if (from[i] != -2)
                {
                    if (i!=j) score[j] = score[i];
                    j++;
                }
            }
            score.resize(j);
        }
    }
    while (!done);
}

void dag_main(vector<Score_t> &score, const string &mol_pair)
{
    int i, n=score.size();

    // should be sorted by y incremental
    Max_Y = score[n-1].y;
    // forward direction
    print_chains(score, mol_pair);
    // reverse complement the second coordinate set.
    for (i=0; i<n; i++)
        score[i].y = Max_Y - score[i].y + 1;
    // reverse direction
    print_chains(score, mol_pair);

    score.clear();
}

