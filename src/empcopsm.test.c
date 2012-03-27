/*
  Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 3 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, see <http://www.gnu.org/licenses/>.
*/


/*****************************************************************************

  Multivariate serial independence test based on the empirical
  copula process

  Ivan Kojadinovic, December 2007

*****************************************************************************/

#include <R.h>
#include <Rmath.h>
#include "set.utils.h"
#include "empcop.stat.h"

/*****************************************************************************

  Array J

******************************************************************************/

void J_sm(int n, int p, int q, const double U[], const int B[], double *J)
{
  int i, j, k, l, m, np = n + p - 1;

  m=0;
  for (j=0;j<p;j++)
    for (l=0;l<n;l++)
      for (i=0;i<n;i++)
	{
	  J[m] = 1.0;
	  for (k=0;k<q;k++)
	    J[m] *= 1.0 - fmax2(U[np * k + B[i + j]], U[np * k + B[l + j]]);
	  m++;
	}
}

/*****************************************************************************

  Bootstrap/permutation of the MAn, up to subsets of cardinality p containing 1
  and of In
  n + p - 1: sample size
  N: number of repetitions
  p: number of lags + 1
  m: max. card. of A
  MA0: bootstrap values of MAn under serial independence (N repetitions)
  I0: bootstrap values of In under serial independence (N repetitions)
  subset: subsets of {1,...,p} in binary notation (int) whose card. is
  between 2 and m in "natural" order and that contain 1
  subset_char: similar, for printing

******************************************************************************/

void bootstrap_serial(int *n, int *N, int *p, int *q, double *U, int *m,
		      double *MA0, double *I0, int *subset, char **subset_char,
		      int *pe)
{
  int i, k, np = *n + *p - 1, p1[1], m1[1], sb[1];
  int *B = Calloc(np, int);
  double *J = Calloc((*n) * (*n) * (*p), double);
  double *K = Calloc((*n) * (*p), double);
  double *L = Calloc(*p, double);

  /* for the random permutation */
  int t1, t2;

  /* number of subsets 2^{p-1}-1 */
  *p1 = *p - 1;
  *m1 = *m - 1;
  *sb = (int)sum_binom(*p1,*m1);

  /* partial power set in "natural" order */
  k_power_set(p1, m1, subset);

  /* shift to the left and add the element 1 */
  for (i=0;i<*sb;i++)
    subset[i] = (subset[i] << 1) + 1;

  /* convert partial power set to char for printing */
  k_power_set_char(p, sb, subset, subset_char);

  /* Simulate the distribution of TA under independence */
  GetRNGstate();

  /* N repetitions */
  for (k=0;k<*N;k++) {

    if ((*pe > 0) && ((k+1) % (*pe) == 0))
      Rprintf("Simulation iteration %d\n",k+1);

    /* identity row selection */
    for (i=0;i<np;i++)
      B[i] = i;
    /* random permutation */
    for (i=np-1;i>=0;i--)
      {
	t1 = B[i];
	t2 = (int)( (i + 1) *  unif_rand());
	B[i] = B[t2];
	B[t2] = t1;
      }

    /* compute arrays J, K, L */
    J_sm(*n, *p, *q, U, B, J);
    K_array(*n, *p, J, K);
    L_array(*n, *p, K, L);

    /* for subsets i of cardinality greater than 1,
       containing 1, and lower than m */
    /* subset[i+1]: i+1 because do not need {1}; truncated after */
    for (i=0;i<*sb-1;i++)
      MA0[k + (*N) * i] = M_A_n(*n, *p, J, K, L, subset[i+1]);

    /* global statistic */
    I0[k] = I_n(*n, *p, J, K, L);
  }
  PutRNGstate();

  Free(B);

  Free(J);
  Free(K);
  Free(L);
}

/*****************************************************************************

  Compute MAn for subsets of cardinality 2 to p containing 1
  U: ranks/np (pseudo-obs)
  n + p - 1: sample size
  p: number of lags + 1
  q: dimension of the time series
  m: max. cardinality of subsets of {1,...,p}
  MA0: simulated values of TA under independence (size: N * (sb - 1))
  where sb = sum_binom(p-1,m-1)
  N: number of repetitions (nrows MA0, fisher0, tippett0)
  subset: subsets of {1,...,p} in binary notation (int) whose card. is
  between 2 and m in "natural" and order containing 1
  MA: test statistics (size: sum.bin. - 1)
  pval: corresponding p-values (size = sum.bin. - 1)
  fisher: pvalue à la Fisher
  tippett: pvalue à la Tippett

******************************************************************************/

void empirical_copula_test_rv_serial(double *U, int *n, int *p, int *q, int *m, double *MA0,
				     double *I0, int *N, int *subset, double *MA, double *I,
				     double *pval, double *fisher, double *tippett, double *Ipval)
{
  int i, j, k, count, sb = (int)sum_binom(*p-1,*m-1), np = *n + *p - 1;
  double *fisher0 = Calloc(*N, double);
  double *tippett0 = Calloc(*N, double);
  double *J = Calloc((*n) * (*n) * (*p), double);
  double *K = Calloc((*n) * (*p), double);
  double *L = Calloc(*p, double);
  int *B = Calloc(np, int);
  double pvalue;

  /* identity row selection */
  for (i=0;i<np;i++)
      B[i] = i;

  /* compute W à la Fisher and à la Tippett from MA0*/
  for (k=0;k<*N;k++)
    {
      fisher0[k] = 0.0;
      tippett0[k] = 1.0;
      for (i=0;i<sb-1;i++)
	{
	  /* p-value */
	  count = 0;
	  for (j=0;j<*N;j++)
	    if (MA0[j + (*N) * i] >= MA0[k + (*N) * i])
	      count ++;
	  pvalue = (double)(count + 0.5)/(*N + 1.0);
	  fisher0[k] -= 2*log(pvalue);
	  tippett0[k] = fmin2(tippett0[k],pvalue);
	}
    }

  /* compute W from the current data */
  *fisher = 0.0;
  *tippett = 1.0;

  /* compute arrays J, K, L */
  J_sm(*n, *p, *q, U, B, J);
  K_array(*n, *p, J, K);
  L_array(*n, *p, K, L);

  /* for subsets i of cardinality greater than 2 containing 1 */
  for (i=0;i<sb-1;i++)
    {
      MA[i] = M_A_n(*n, *p, J, K, L, subset[i]);

      /* p-value */
      count = 0;
      for (k=0;k<*N;k++)
	if (MA0[k + (*N) * i] >= MA[i])
	  count ++;
      pval[i] = (double)(count + 0.5)/(*N + 1.0);

      *fisher -= 2*log(pval[i]);
      *tippett = fmin2(*tippett,pval[i]);
    }

  /* p-values of the Fisher and Tippett statistics */
  count = 0;
  for (k=0;k<*N;k++)
    if (fisher0[k] >= *fisher)
      count ++;
  *fisher = (double)(count + 0.5)/(*N + 1.0);

  count = 0;
  for (k=0;k<*N;k++)
    if (tippett0[k] <= *tippett)
      count ++;
  *tippett = (double)(count + 0.5)/(*N + 1.0);

  /* compute In from the current data and the corresponding pvalue*/
  *I = I_n(*n, *p, J, K, L);
  count = 0;
  for (k=0;k<*N;k++)
    if (I0[k] >= *I)
      count ++;
  *Ipval = (double)(count + 0.5)/(*N + 1.0);

  Free(fisher0);
  Free(tippett0);
  Free(J);
  Free(K);
  Free(L);
  Free(B);
}


