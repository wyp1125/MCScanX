/*
 * Author: Haibao Tang <bao@uga.edu> May 17, 2007
 *
 * Simple output subroutines, controls what to print in the final .blocks file.
 * Also contains some helper functions for other data structures.
 *
 * Modified by Yupeng Wang, Mar 31, 2011
 * Some functions were dropped in MCScanX
 */

#include "out_utils.h"

void print_align(FILE* fw)
/* print pair-wise alignment */
{
    int i, j, pid;
    int nseg = seg_list.size(), nanchor;
    Seg_feat *s;
    string sp1,sp2,spc;

    print_params(fw);
    //////////////////////////////////////////////////
    set<string> colgenes;
    for (i=0; i<nseg; i++)
    {
        s = &seg_list[i];
        nanchor = s->pids.size();
        for (j=0; j<nanchor; j++)
        {
            pid = s->pids[j];
            colgenes.insert(match_list[pid].gene1);
            colgenes.insert(match_list[pid].gene2);
        }
    }
    fprintf( fw, "############### Statistics ###############\n");
    double temp=100*(double)colgenes.size()/(double)gene_map.size();
    fprintf(fw,"# Number of collinear genes: %d, Percentage: %.2f\n",(int)colgenes.size(),temp);
    fprintf(fw,"# Number of all genes: %d\n", (int)gene_map.size());
    fprintf( fw, "##########################################\n");
//////////////////////////////////////////////////


    for (i=0; i<nseg; i++)
    {
        s = &seg_list[i];
        nanchor = s->pids.size();
        fprintf(fw, "## Alignment %d: score=%.1f e_value=%.2g N=%d %s %s\n",
                i, s->score, s->e_value, nanchor, s->mol_pair.c_str(),
                s->sameStrand?"plus":"minus");
        sp1=(s->s1)->mol.substr(0,2);
        sp2=(s->s2)->mol.substr(0,2);
        spc=sp1+"&"+sp2; 
        cmp_sp[spc].syn_num+=nanchor;
        
        for (j=0; j<nanchor; j++)
        {
            pid = s->pids[j];
            fprintf(fw, "%3d-%3d:\t%s\t%s\t%7.1g\n",
                    i, j, match_list[pid].gene1.c_str(),
                    match_list[pid].gene2.c_str(), match_list[pid].score);
        }
    }
}

void print_params(FILE *fw)
/* print parameters */
{
    fprintf( fw, "############### Parameters ###############\n");
    fprintf( fw, "# MATCH_SCORE: %d\n", MATCH_SCORE );
    fprintf( fw, "# MATCH_SIZE: %d\n", MATCH_SIZE );
    //fprintf( fw, "# UNIT_DIST: %d\n", UNIT_DIST );
    fprintf( fw, "# GAP_PENALTY: %d\n", GAP_PENALTY );
    fprintf( fw, "# OVERLAP_WINDOW: %d\n", OVERLAP_WINDOW );
    //fprintf( fw, "# EXTENSION_DIST: %d\n", EXTENSION_DIST );
    fprintf( fw, "# E_VALUE: %lg\n", E_VALUE );
    fprintf( fw, "# MAX GAPS: %d\n", MAX_GAPS );
    //fprintf( fw, "##########################################\n\n");
}

