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

  Independence test among random vectors based on the empirical
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

void J_m(int n, int p, const int b[], const double U[], const int R[],
	 double *J)
{
  int m=0;
  for (int k=0; k < p; k++)
    for (int l=0; l < n; l++)
      for (int i=0; i < n; i++)
	{
	  J[m] = 1.0;
	  for (int j= b[k]; j < b[k+1]; j++)
	    J[m] *= 1.0 - fmax2(U[n * j + R[n * k + i]],
				U[n * j + R[n * k + l]]);
	  m++;
	}
}

/*****************************************************************************

  Bootstrap of the MAn, up to subsets of cardinality p
  and of In
  n: sample size
  N: number of repetitions
  p: dimension of data
  m: max. card. of A
  MA0: bootstrap values of MAn under independence (N repetitions)
  I0: bootstrap values of In under independence (N repetitions)
  subset: subsets of {1,...,p} in binary notation (int) whose card. is
  between 2 and m in "natural" order
  subset_char: similar, for printing

******************************************************************************/

void bootstrap(int *n, int *N, int *p, int *b, double *U, int *m,
	       double *MA0, double *I0, int *subset, char **subset_char,
	       int *pe)
{
  int i, j, k, sb[1];
  int *R = Calloc((*n) * (*p), int);
  double *J = Calloc((*n) * (*n) * (*p), double);
  double *K = Calloc((*n) * (*p), double);
  double *L = Calloc(*p, double);

  /* for the random permutation */
  int t1, t2;

  /* number of subsets */
  *sb = (int)sum_binom(*p,*m);

  /* partial power set in "natural" order */
  k_power_set(p, m, subset);

  /* convert partial power set to char for printing */
  k_power_set_char(p, sb, subset, subset_char);

  /* Simulate the distribution of TA under independence */
  GetRNGstate();

  /* N repetitions */
  for (k=0;k<*N;k++) {

    if ((*pe > 0) && ((k+1) % (*pe) == 0))
      Rprintf("Simulation iteration %d\n",k+1);

    /* generate row selection within the blocks */
    /* for (j=0;j<*p;j++)
      for (i=0;i<*n;i++)
      R[(*n) * j + i] = (int)( unif_rand() * (*n) ); */

    /* generate row selection */
    for (j=0;j<*p;j++)
      {
	for (i=0;i<*n;i++)
	  R[i + (*n) * j] = i;

	/* random permutation */
	for (i=*n-1;i>=0;i--)
	  {
	    t1 = R[j * (*n) + i];
	    t2 = (int)((i + 1) *  unif_rand());
	    R[j * (*n) + i] = R[j * (*n) + t2];
	    R[j * (*n) + t2] = t1;
	  }
      }

    /* compute arrays J, K, L */
    J_m(*n, *p, b, U, R, J);
    K_array(*n, *p, J, K);
    L_array(*n, *p, K, L);

    /* for subsets i of cardinality greater than 1 and lower than m */
    for (i=*p+1;i<*sb;i++)
      MA0[k + (*N) * (i - *p - 1)] =  M_A_n(*n, *p, J, K, L, subset[i]);
    /* global statistic */
    I0[k] = I_n(*n, *p, J, K, L);
  }
  PutRNGstate();

  Free(R);

  Free(J);
  Free(K);
  Free(L);
}

/*****************************************************************************

  Compute MAn for subsets of cardinality 2 to p
  U: pseudo-observations
  n: sample size
  p: number of random vectors
  b: vector of vector dimensions (p blocks)
  m: max. cardinality of subsets of {1,...,p}
  MA0: simulated values of TA under independence (size: N * (sb - p - 1))
  N: number of repetitions (nrows MA0, fisher0, tippett0)
  subset: subsets of {1,...,p} in binary notation (int) whose card. is
  between 2 and m in "natural" order
  MA: test statistics (size: sum.bin. - p - 1)
  pval: corresponding p-values (size = sum.bin. - p - 1)
  fisher: pvalue à la Fisher
  tippett: pvalue à la Tippett

******************************************************************************/

void empirical_copula_test_rv(double *U, int *n, int *p, int *b, int *m, double *MA0,
			      double *I0, int *N, int *subset, double *MA, double *I,
			      double *pval, double *fisher, double *tippett, double *Ipval)
{
  int i, j, k, count, sb = (int)sum_binom(*p,*m);
  double *fisher0 = Calloc(*N, double);
  double *tippett0 = Calloc(*N, double);
  double *J = Calloc((*n) * (*n) * (*p), double);
  double *K = Calloc((*n) * (*p), double);
  double *L = Calloc(*p, double);
  int *R = Calloc((*n) * (*p), int);
  double pvalue;

  /* generate identity selection within the blocks */
  /* used for the computation of the test statistics */
  for (j=0;j<*p;j++)
    for (i=0;i<*n;i++)
      R[(*n) * j + i] = i;

  /* compute W à la Fisher and à la Tippett from MA0*/
  for (k=0;k<*N;k++) {
    fisher0[k] = 0.0;
    tippett0[k] = 1.0;
    for (i=0;i<sb-*p-1;i++) {
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
  J_m(*n, *p, b, U, R, J);
  K_array(*n, *p, J, K);
  L_array(*n, *p, K, L);

  /* for subsets i of cardinality greater than 1 */
  for (i=0;i<sb-*p-1;i++)
    {
      MA[i] =  M_A_n(*n, *p, J, K, L, subset[i]);

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
  Free(R);
}


