/*
 * Author: Haibao Tang (Mar. 7, 2008)
 *
 * Statistical evaluation for syntenic block, see Wang et al. 2006
 *
 * Statistical inference of chromosomal homology based on gene colinearity
 * and applications to Arabidopsis and rice, BMC Bioinformatics (7) pg. 447
 *
 * The expectation of a pairwise syntenic block is
 *
 * E = 2*P(N,m)*Multiply((l_1i*l_2i)/(L1*L2))
 *
 * where N is the number of where N is the number of matching gene pairs
 * (by BLASTP or BLAT, etc.) between two chromosomal regions; m is the
 * number of matching gene pairs in the identified colinearity pattern;
 * L1 and L2 are respective lengths of the two chromosomal regions;
 * and l1i and l2i are distances between two adjacent collinear gene pairs
 * in the syntenic block. The expectation multiplies by two since there are
 * two possible orientation configurations between two collinear segments.
 *
 * Permutations code reference:
 * http://www.ciphersbyritter.com/JAVASCRP/PERMCOMB.HTM
 *
 */

#include "permutation.h"

static double fact(int x)
/* returns factorial of x */
{
    double ans = 1;
    while (x > 1) ans *= x--;
    return ans;
}

static double ln_fact(int x)
/* ln(x!) using Stirling's formula, see Knuth I: 111 */
{
    double dx = x, invx, invx2, invx3, invx5, invx7, sum;

    if (x < 12) return log(fact(x));
    else
    {
        invx = 1 / dx;
        invx2 = invx * invx;
        invx3 = invx2 * invx;
        invx5 = invx3 * invx2;
        invx7 = invx5 * invx2;

        sum = ((dx + 0.5) * log(dx)) - dx;
        sum += log(2 * M_PI) / 2;
        sum += invx / 12 - invx3 / 360;
        sum += invx5/ 1260 - invx7 / 1680;

        return sum;
    }
}

double ln_perm(int n, int r)
/* natural log of permutation number */
{
    if (r > n || r <= 0) return 0;
    return ln_fact(n) - ln_fact(n-r);
}

double ln_comb(int n, int k)
/* natural log of combination number */
{
    if (k <= 0 || k >= n) return 0;
    return ln_fact(n) - ln_fact(k) - ln_fact(n-k);
}

