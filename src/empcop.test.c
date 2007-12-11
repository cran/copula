/*#################################################################
##   Copula R package by Jun Yan and Ivan Kojadinovic Copyright (C) 2007
##
##   Copyright (C) 2007 Ivan Kojadinovic <ivan@stat.auckland.ac.nz>
##
##   This program is free software; you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation; either version 2 of the License, or
##   (at your option) any later version.
##
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License along
##   with this program; if not, write to the Free Software Foundation, Inc.,
##   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
##
#################################################################*/

/*****************************************************************************

  Multivariate independence test  based on the empirical 
  copula process as proposed by Christian Genest and Bruno 
  Rémillard (2004), Test 13:2, pages 335-369.

  Ivan Kojadinovic, May 2007, modified Dec 2007

*****************************************************************************/

#include <R.h>
#include <Rmath.h>
#include "set.utils.h"

/*****************************************************************************

  Computation of the Cramer-von Mises statistics used in the independence
  test of Genest and Remillard based on the empirical copula for a
  subset A of {1,...,p}
  
******************************************************************************/

double global_stat(int n, int p, int *R)
{
  int i,j,l;
  double stat, sum, prod, part1, part2;
  
  /* first term */
  sum = 0.0;
  for (i=0;i<n;i++)
    for (l=0;l<n;l++) 
      {
	prod = 1.0;
	for (j=0;j<p;j++)
	  prod *= 1.0 - fmax2(R[i + n*j], R[l + n*j])/n;
	sum += prod;
      }
  part1 = sum/n;
  
  /* second term */
  sum = 0.0;
  for (i=0;i<n;i++)
    {
      prod = 1.0;
      for (j=0;j<p;j++)
	prod *= (n * (n - 1.0) - R[i + n*j] * (R[i + n*j] - 1.0) )
	  /(2.0 * n * n);
      sum += prod;
    }
  part2 = sum;

  stat = part1 - 2.0 * part2 
    + n * R_pow_di((n - 1.0) * (2.0 * n - 1.0) / (6.0 * n * n) ,p);

  return stat;
}

double Dn(int n, int s, int t)
{
  return (2.0 * n + 1.0) * (n + 1.0) / (6.0 * n * n) 
    + (s * (s - 1.0) + t * (t - 1.0)) / (2.0 * n * n) 
    - fmax2(s,t) / n;
}

double empirical_copula_subset(int n, int p, int *R, int A)
{
  int i,j,k;
  double TA = 0.0, prod;
  
  for (i=0;i<n;i++)
    for (k=0;k<n;k++) {
      prod = 1.0;
      for (j=0;j<p;j++)
	if (1<<j & A)
	  prod *= Dn(n, R[i + n*j], R[k + n*j]);
      TA += prod;
    }
  return TA/n;
}

/*****************************************************************************

  Simulate the distribution of TA, up to subsets of cardinality p
  and of the global statistic
  n: sample size
  N: number of repetitions
  p: dimension of data
  m: max. card. of A
  TA0: values of TA under independence (N repetitions) 
  G0: values of the global stat. under independence (N repetitions)
  subset: subsets of {1,...,p} in binary notation (int) whose card. is 
  between 2 and m in "natural" order
  subset_char: similar, for printing
  
******************************************************************************/

void simulate_empirical_copula(int *n, int *N, int *p, int *m, double *TA0, 
			       double *G0, int *subset, char **subset_char)
{
  int i, j, k, t1, t2, sb = (int)sum_binom(*p,*m);
  int *R = (int *)R_alloc((*n) * (*p), sizeof(int));
  
  /* partial power set in "natural" order */ 
  k_power_set(p, m, subset);

  /* convert partial power set to char for printing */
  k_power_set_char(p,m, subset, subset_char);
  
  /* Simulate the distribution of TA under independence */

  GetRNGstate();
  
  /* N repetitions */
  for (k=0;k<*N;k++) { 

    Rprintf("Simulation iteration %d\n",k+1);
    
    /* generate data */
    for (j=0;j<*p;j++) 
      {	
	for (i=0;i<*n;i++)
	  R[i + (*n) * j] = i+1;
	
	/* permutation = random ranks in column j */
	for (i=*n-1;i>=0;i--)
	  {
	    t1 = R[j * (*n) + i]; 
	    t2 = (int)(i *  unif_rand());
	    R[j * (*n) + i] = R[j * (*n) + t2];
	    R[j * (*n) + t2] = t1;
	  }
      }
    
    /* for subsets i of cardinality greater than 1 and lower than m */
    for (i=*p+1;i<sb;i++) 
      TA0[k + (*N) * (i - *p - 1)] = empirical_copula_subset(*n, *p, R, subset[i]);

    /* global stat under independence*/
    G0[k] = global_stat(*n, *p, R);
  }
  PutRNGstate();
}

/*****************************************************************************

  Compute TA for subsets of cardinality 2 to p
  R: ranks
  n: sample size
  p: number of variables
  m: max. cardinality of subsets of {1,...,p}
  TA0: simulated values of TA under independence (size: N * (sb - p - 1))
  N: number of repetitions (nrows TA0, fisher0, tippett0)
  subset: subsets of {1,...,p} in binary notation (int) whose card. is 
  between 2 and m in "natural" order
  TA: test statistics (size: sum.bin. - p - 1) 
  pval: corresponding p-values (size = sum.bin. - p - 1)
  fisher: pvalue à la Fisher
  tippett: pvalue à la Tippett
  
******************************************************************************/

void empirical_copula_test(int *R, int *n, int *p, int *m, double *TA0, double *G0,
			   int *N, int *subset, double *TA, double *G, double *pval,
			   double *fisher, double *tippett, double *globpval)			  
{
  int i, j, k, count, sb = (int)sum_binom(*p,*m);
  double *fisher0 = (double *)R_alloc(*N, sizeof(double));
  double *tippett0 = (double *)R_alloc(*N, sizeof(double));
  double pvalue;

  /* compute W à la Fisher and à la Tippett from TA0*/ 
  for (k=0;k<*N;k++) {
    fisher0[k] = 0.0;
    tippett0[k] = 1.0;
    for (i=*p+1;i<sb;i++) {
      /* p-value */
      count = 1;
      for (j=0;j<*N;j++) 
	if (TA0[j + (*N) * (i - *p - 1)] > TA0[k + (*N) * (i - *p - 1)])
	  count ++;
      pvalue = (double)(count + 0.5)/(*N + 1.0);
      fisher0[k] -= 2*log(pvalue);
      tippett0[k] = fmin2(tippett0[k],pvalue);
    }
  }

  /* compute W from the current data */
  *fisher = 0.0;
  *tippett = 1.0;
  
  /* for subsets i of cardinality greater than 1 */
  for (i=*p+1;i<sb;i++) 
    {
      TA[i - *p - 1] = empirical_copula_subset(*n, *p, R, subset[i - *p - 1]);
      
      /* p-value */
      count = 1;
      for (k=0;k<*N;k++) 
	if (TA0[k + (*N) * (i - *p - 1)] > TA[i - *p - 1])
	  count ++;
      pval[i - *p - 1] = (double)(count + 0.5)/(*N + 1.0);
      
      *fisher -= 2*log(pval[i - *p - 1]);
      *tippett = fmin2(*tippett,pval[i - *p - 1]);
    }

  /* p-values of the Fisher and Tippett statistics */
  count = 1;
  for (k=0;k<*N;k++) 
    if (fisher0[k] > *fisher)
      count ++;
  *fisher = (double)(count)/(*N + 1.0);

  count = 1;
  for (k=0;k<*N;k++) 
    if (tippett0[k] < *tippett)
      count ++;
  *tippett = (double)(count)/(*N + 1.0);


  /* compute global stat from the current data and the corresponding p-value*/ 
  *G = global_stat(*n, *p, R); 
  count = 1;
  for (k=0;k<*N;k++) 
    if (G0[k] > *G)
      count ++;
  *globpval = (double)(count)/(*N + 1.0);
}


