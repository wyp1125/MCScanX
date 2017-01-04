/*
 * Author: Haibao Tang <bao@uga.edu> May 10, 2007
 * Main entry point for the executable mcscan
 *
 * Modified by Yupeng Wang, Mar 31, 2011
 * The parameter MAX_GAPS was added. IN_SYNTENY has 3 choices and is functional.
*/

#include "mcscan_h.h"
#include <iostream>
using namespace std;

static bool IS_PAIRWISE;
static bool BUILD_MCL;
static char prefix_fn[LABEL_LEN];

static void print_banner()
{
    ;
}

static void init_opt()
{
    // match bonus, final score=MATCH_SCORE-GAPS^GAP_PENALTY
    MATCH_SCORE = 50;
    // the number of genes required to call synteny, sometimes more
    MATCH_SIZE = 5;
    // gap extension penalty
    GAP_PENALTY = -1;
    // alignment significance
    E_VALUE = 1e-5;
    // maximum gaps allowed
    MAX_GAPS =25;
    OVERLAP_WINDOW=5;
    e_mode=0;
    IS_PAIRWISE = false;
    IN_SYNTENY = 2;
}

static void print_help(const char *prg)
{
     progress("[Usage] %s prefix_fn [options]\n"
             " -k  MATCH_SCORE, final score=MATCH_SCORE+NUM_GAPS*GAP_PENALTY\n"
             "     (default: %d)\n"
             " -g  GAP_PENALTY, gap penalty (default: %d)\n"
             " -s  MATCH_SIZE, number of genes required to call collinear blocks\n"
             "     (default: %d)\n"
             " -e  E_VALUE, alignment significance (default: %lg)\n"
             " -m  MAX_GAPS, maximum gaps allowed (default: %d)\n"
             " -w  OVERLAP_WINDOW, maximum distance (# of genes) to collapse BLAST matches (default: %d)\n"
             " -a  only builds the pairwise blocks (.collinearity file)\n"
             " -b  patterns of collinear blocks. 0:intra- and inter-species (default); 1:intra-species; 2:inter-species\n"
             " -c  whether to consider homology scores. 0:not consider (default); 1: lower preferred; 2: higher preferred\n"
             " -h  print this help page\n",
             prg, MATCH_SCORE, GAP_PENALTY, MATCH_SIZE, E_VALUE, MAX_GAPS, OVERLAP_WINDOW);
    exit(1);

}

static void read_opt(int argc, char *argv[])
{
    int c;
    opterr = 0;

    if (argc < 2) print_help(argv[0]);

    while ((c = getopt(argc, argv, "k:g:s:e:b:m:w:c:ah")) != -1)
        switch (c)
        {
        case 'k':
            MATCH_SCORE = atoi(optarg);
            break;
        case 'g':
            GAP_PENALTY = atoi(optarg);
            break;
        case 's':
            MATCH_SIZE = atoi(optarg);
            break;
        case 'e':
            E_VALUE = atof(optarg);
            break;
        case 'b':
            IN_SYNTENY = atoi(optarg);
            break;
        case 'm':
            MAX_GAPS = atoi(optarg);
            break;
        case 'w':
            OVERLAP_WINDOW = atoi(optarg);
            break;
        case 'c':
            e_mode = atoi(optarg);
            break;
        case 'a':
            IS_PAIRWISE = true;
            break;
        case '?':
            if (optopt=='k' || optopt=='s' || optopt=='g' || optopt=='e' || optopt=='b' || optopt=='m' || optopt=='w' || optopt=='c')
                errAbort("Option -%c requires an argument.", optopt);
            else if (isprint (optopt))
                errAbort("Unknown option `-%c'.", optopt);
            else
                errAbort("Unknown option character `\\x%x'.", optopt);
        default:
            print_help(argv[0]);
            break;
        }

    if (optind==argc) errAbort("Please enter your input file");
    else strcpy(prefix_fn, argv[optind]);
    CUTOFF_SCORE = MATCH_SCORE*MATCH_SIZE;
}
void fill_allg()
{
    Gene_feat *gf1;
    map<string, Gene_feat>::iterator it;
    for (it=gene_map.begin(); it!=gene_map.end(); it++)
    {
        gf1 = &(it->second);
        allg.insert(gf1);
    }
    int i=0;
    geneSet::const_iterator it77=allg.begin();   
    for (; it77!=allg.end(); it77++)
    {
    (*it77)->gene_id=i;
    i++;
    }
}
void show_stat()
{
    map<string, ortho_stat>::iterator it19;
    cout<<"Print statistics:"<<endl;
    cout<<"Species\t# of collinear homolog pairs\t# of homolog pairs\tPercentage"<<endl;
    for(it19=cmp_sp.begin();it19!=cmp_sp.end();it19++)
    {
    double temp=100.0*(double)it19->second.syn_num/(double)it19->second.all_num;
    cout<<it19->first<<"\t"<<it19->second.syn_num<<"\t"<<it19->second.all_num<<"\t"<<temp<<endl;
    }
}
int main(int argc, char *argv[])
{
    /* Start the timer */
    uglyTime(NULL);

    print_banner();
    char align_fn[LABEL_LEN], block_fn[LABEL_LEN];
    FILE *fw;

    init_opt();
    read_opt(argc, argv);

    read_gff(prefix_fn);
    read_orthomcl(prefix_fn);

    sprintf(align_fn, "%s.collinearity", prefix_fn);
    fw = mustOpen(align_fn, "w");

    progress("%d pairwise comparisons", (int) mol_pairs.size());
    fill_allg();
    map<string, int>::const_iterator ip;
    for (ip=mol_pairs.begin(); ip!=mol_pairs.end(); ip++)
    {
        if (ip->second >= MATCH_SIZE) feed_dag(string(ip->first));
    }

    progress("%d alignments generated", (int) seg_list.size());
    print_align(fw);
    fclose(fw);
    uglyTime("Pairwise collinear blocks written to %s", align_fn);

    if(IS_PAIRWISE)
    {
    show_stat();
    return 0;
    }
    msa_main(prefix_fn);
    show_stat();
    uglyTime("Done!");

    return 0;
}

