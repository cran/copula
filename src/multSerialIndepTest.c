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

/**
 * @file   multSerialIndepTest.c
 * @author Ivan Kojadinovic
 * @date   December 2007
 *
 * @brief Multivariate serial independence test based on the empirical
 *        copula process -- see AISM paper
 *
 */


#include <R.h>
#include <Rmath.h>
#include "set_utils.h"
#include "indepTests.h"

/// Temporary array J
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

/**
 * Bootstrap/permutation of the MAn, up to subsets of cardinality p
 * containing 1 and of In
 *
 * @param n n+p-1 is the available sample size
 * @param N number of repetitions
 * @param p number of lags + 1
 * @param q number of columns of U
 * @param U pseudo-observations
 * @param m max. card. of A
 * @param MA0 bootstrap values of MAn under serial independence (N repetitions)
 * @param I0 bootstrap values of In under serial independence (N repetitions)
 * @param subset subsets of {1,...,p} in binary notation (int) whose card. is
 *               between 2 and m in "natural" order and that contain 1
 * @param subset_char similar, for printing
 * @param verbose display progress bar if > 0
 * @author Ivan Kojadinovic
 */
void bootstrap_serial(int *n, int *N, int *p, int *q, double *U, int *m,
		      double *MA0, double *I0, int *subset, char **subset_char,
		      int *verbose)
{
  int i, k, np = *n + *p - 1, p1[1], m1[1], sb[1];
  size_t max_size = (size_t)-1,// C99 has SIZE_MAX
      n_ = (size_t)(*n);
  double J_size = ((double)n_) * ((double)n_) * (*p);
  if(J_size > max_size)
      error(_("** bootstrap_serial(): n or p too large: n^2*p = %12.0g > %12.0g = max(size_t)\n"),
	    J_size, (double)max_size);

  int *B = Calloc(np, int);
  double *J = Calloc((size_t) J_size, double);
  double *K = Calloc(n_ * (*p), double);
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

    if (*verbose)
      progressBar(k, *N, 70);
  }
  PutRNGstate();

  Free(B);

  Free(J);
  Free(K);
  Free(L);
}

/**
 * Compute the MAn for subsets of cardinality 2 up to p containing 1, as
 * well as the global statistic I
 *
 * @param U pseudo-observations
 * @param n n+p-1 is the available sample size
 * @param p number of lags + 1
 * @param q dimension of the time series
 * @param m max. cardinality of subsets of {1,...,p}
 * @param MA0 simulated values of MA under independence (size: N * (sb - 1))
 * @param I0 simulated values of I under independence
 * @param N number of repetitions (nrows MA0, fisher0, tippett0)
 * @param subset subsets of {1,...,p} in binary notation (int) whose card. is
 *               between 2 and m in "natural" and order containing 1
 * @param MA test statistics (size: sum.bin. - 1)
 * @param I global test statistic
 * @param pval p-values for the MA (size = sum.bin. - 1)
 * @param fisher pvalue à la Fisher
 * @param tippett pvalue à la Tippett
 * @param Ipval pvalue for I
 * @author Ivan Kojadinovic
 */
void empirical_copula_test_rv_serial(double *U, int *n, int *p, int *q, int *m, double *MA0,
				     double *I0, int *N, int *subset, double *MA, double *I,
				     double *pval, double *fisher, double *tippett, double *Ipval)
{
  int i, j, k, count, sb = (int)sum_binom(*p-1,*m-1), np = *n + *p - 1;
  size_t max_size = (size_t)-1,// C99 has SIZE_MAX
      n_ = (size_t)(*n);
  double J_size = ((double)n_) * ((double)n_) * (*p);
  if(J_size > max_size)
      error(_("** empirical_copula_t.r.s(): n or p too large: n^2*p = %12.0g > %12.0g = max(size_t)\n"),
	    J_size, (double)max_size);

  double *fisher0 = Calloc(*N, double);
  double *tippett0 = Calloc(*N, double);
  double *J = Calloc((size_t) J_size, double);
  double *K = Calloc(n_ * (*p), double);
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


