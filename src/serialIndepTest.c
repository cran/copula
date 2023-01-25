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
 * @file   serialIndepTest.c
 * @author Ivan Kojadinovic
 * @date   December 2007
 *
 * @brief Serial independence test  based on the empirical
 *        copula process as proposed by Christian Genest and Bruno
 *        Rémillard (2004), Test 13:2, pages 335-369.
 *
 */


#include <R.h>
#include <Rmath.h>
#include "set_utils.h"
#include "indepTests.h"


/// Temporary array J
void J_s(int n, int p, const double U[], double *J)
{
  int m=0;
  for (int j=0;j<p;j++)
    for (int l=0;l<n;l++)
      for (int i=0;i<n;i++)
	J[m++] = 1.0 - fmax2(U[i + j], U[l + j]);
}

/**
 * Simulate the distribution of the TAs, up to subsets of cardinality p
 * containing 1, and of the global statistic
 *
 * @param n n+p-1 is the sample size
 * @param N number of repetitions
 * @param p number of lags + 1
 * @param m max. card. of A
 * @param TA0 values of TAs under serial independence (N repetitions)
 * @param G0 values of the global stat. under serial independence (N repetitions)
 * @param subset subsets of {1,...,p} in binary notation (int) whose card. is
 *               between 2 and m in "natural" order and that contain 1
 * @param subset_char similar, for printing
 * @param fisher0 p-values à la Fisher
 * @param tippett0 p-values à la Tippett
 * @param verbose display progress bar if > 0
 * @author Ivan Kojadinovic
 */
void simulate_empirical_copula_serial(int *n, int *N, int *p, int *m,
				      double *TA0, double *G0, int *subset,
				      char **subset_char, double *fisher0,
				      double *tippett0, int *verbose)
{
  int i, j, k, np = *n + *p - 1, p1[1], m1[1], sb[1], count, index;
  size_t max_size = (size_t)-1,// C99 has SIZE_MAX
      n_ = (size_t)(*n);
  double J_size = ((double)n_) * ((double)n_) * (*p);
  if(J_size > max_size)
      error(_("** simulate_emp.cop.serial(): n or p too large: n^2*p = %12.0g > %12.0g = max(size_t)\n"),
	    J_size, (double)max_size);

  double *U = Calloc(np, double);
  double *J = Calloc((size_t) J_size, double);
  double *K = Calloc(n_ * (*p), double);
  double *L = Calloc(*p, double);
  double pvalue, u;

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

    /* generate data */
    for (i=0;i<np;i++)
      U[i] = (i+1)/(double)np;

    /* permutation = random ranks */
    for (i=np-1;i>=0;i--)
      {
	u = U[i];
	index = (int)((i + 1) *  unif_rand());
	U[i] = U[index];
	U[index] = u;
      }

    /* compute arrays J, K, L */
    J_s(*n, *p, U, J);
    K_array(*n, *p, J, K);
    L_array(*n, *p, K, L);

    /* for subsets i of cardinality greater than 1 and lower than m */
    /* subset[i+1]: i+1 because do not need {1}; truncated after */
    for (i=0;i<*sb-1;i++)
      TA0[k + (*N) * i] =  M_A_n(*n, *p, J, K, L, subset[i+1]);

    /* global stat under independence*/
    G0[k] = I_n(*n, *p, J, K, L);

    if (*verbose)
      progressBar(k, *N, 70);
  }

  PutRNGstate();

  /* compute W à la Fisher and à la Tippett from TA0*/
  for (k=0;k<*N;k++) {
    fisher0[k] = 0.0;
    tippett0[k] = 1.0;
    for (i=0;i<*sb-1;i++) {
      /* p-value */
      count = 0;
      for (j=0;j<*N;j++)
	if (TA0[j + (*N) * i] >= TA0[k + (*N) * i])
	  count ++;
      pvalue = (double)(count + 0.5)/(*N + 1.0);
      fisher0[k] -= 2*log(pvalue);
      tippett0[k] = fmin2(tippett0[k],pvalue);
    }
  }

  Free(U);

  Free(J);
  Free(K);
  Free(L);
}

/**
 * Compute the TAs for subsets of cardinality 2 to p containing 1,
 * as well as In
 *
 * @param U univariate pseudo-observations
 * @param n n+p-1 is the sample size
 * @param p number of lags + 1
 * @param m max. cardinality of subsets of {1,...,p}
 * @param TA0 simulated values of TAs under serial independence (size: N * (sb - 1))
 * @param G0 simulated values of G under serial independence
 * @param N number of repetitions (nrows TA0, fisher0, tippett0)
 * @param subset subsets of {1,...,p} in binary notation (int) whose card. is
 *               between 2 and m in "natural" order and that contain 1
 * @param TA test statistics (size: sum.bin. - 1)
 * @param G global statistic
 * @param pval p-values corresponding to the TAs (size = sum.bin. - 1)
 * @param fisher pvalue à la Fisher
 * @param tippett pvalue à la Tippett
 * @param globpval pvalue of G
 * @param fisher0 pvalues à la Fisher under the null
 * @param tippett0 pvalues à la Tippett under the null
 * @author Ivan Kojadinovic
 */
void empirical_copula_test_serial(double *U, int *n, int *p, int *m, double *TA0, double *G0,
				  int *N, int *subset, double *TA, double *G, double *pval,
				  double *fisher, double *tippett, double *globpval,
				  double *fisher0, double *tippett0)
{
  int i, k, count, sb = (int)sum_binom(*p-1,*m-1);
  size_t max_size = (size_t)-1,// C99 has SIZE_MAX
      n_ = (size_t)(*n);
  double J_size = ((double)n_) * ((double)n_) * (*p);
  if(J_size > max_size)
      error(_("** emp.cop.test_serial(): n or p too large: n^2*p = %12.0g > %12.0g = max(size_t)\n"),
	    J_size, (double)max_size);

  double *J = Calloc((size_t) J_size, double);
  double *K = Calloc(n_ * (*p), double);
  double *L = Calloc(*p, double);

  /* compute W and T from the current data */
  *fisher = 0.0;
  *tippett = 1.0;

  /* compute arrays J, K, L */
  J_s(*n, *p, U, J);
  K_array(*n, *p, J, K);
  L_array(*n, *p, K, L);

  /* for subsets i of cardinality greater than 1 */
  for (i=0;i<sb-1;i++)
    {
      TA[i] =  M_A_n(*n, *p, J, K, L, subset[i]);

      /* p-value */
      count = 0;
      for (k=0;k<*N;k++)
	if (TA0[k + (*N) * i] >= TA[i])
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

  /* compute global stat from the current data and the corresponding p-value*/
  *G = I_n(*n, *p, J, K, L);
  count = 0;
  for (k=0;k<*N;k++)
    if (G0[k] >= *G)
      count ++;
  *globpval = (double)(count + 0.5)/(*N + 1.0);

  Free(J);
  Free(K);
  Free(L);
}
