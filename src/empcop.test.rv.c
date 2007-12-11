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

  Independence test among random vectors based on the empirical 
  copula process 

  Ivan Kojadinovic, December 2007

*****************************************************************************/

#include <R.h>
#include <Rmath.h>
#include "set.utils.h"

/*****************************************************************************

  Computation of the Cramer-von Mises statistic derived from the independence 
  copula process.
  
******************************************************************************/

double I_n(int n, int p, int *b, double *U)
{
  int i,j,k,l;
  double In, sum, sum2, prod, prod2, part1, part2, part3;

  /* first term */
  sum = 0.0;
  for (i=0;i<n;i++)
    for (l=0;l<n;l++)
      {
	prod = 1.0;
	for (k=0;k<p;k++)
	  for (j=b[k];j<b[k+1];j++)
	    prod *= 1.0 - fmax2(U[n * j + i], U[n * j + l]);
	sum += prod;
      }
  part1 = sum/n;
    
  /* second term */
  sum = 0.0;
  for (i=0;i<n;i++)
    {
      prod = 1.0;
      for (k=0;k<p;k++)
	{
	  sum2 = 0.0;
	  for (l=0;l<n;l++)
	    {
	      prod2 = 1.0;
	      for (j=b[k];j<b[k+1];j++)
		prod2 *= 1.0 - fmax2(U[n * j + i], U[n * j + l]);
	      sum2 += prod2;	      
	    }
	  prod *= sum2;
	}
      sum += prod;
    }
  part2 = 2.0 * sum / R_pow_di(n,p);
  

  /* third term */
  prod = 1.0;
  for (k=0;k<p;k++)
    {
      sum = 0.0;
      for (i=0;i<n;i++)
	for (l=0;l<n;l++)
	  { 
	    prod2 = 1.0;
	    for (j=b[k];j<b[k+1];j++)
	      prod2 *= 1.0 - fmax2(U[n * j + i], U[n * j + l]);
	    sum += prod2;
	  }
      prod *= sum;
    }
  part3 = prod / R_pow_di(n,2*p-1);

  In = part1 - part2 + part3;
 
  return In;
}


/*****************************************************************************

  Computation of the Cramer-von Mises statistics used in the independence
  test. One statistic per subset A of {1,...,p}, |A|>1.
  
******************************************************************************/

void compute_part23(int n, int p, int *b, double *U, double *part23)
{
  int i,j,k,m;
  double prod;

  for (i=0;i<n;i++)
    for (k=0;k<p;k++)
      {
	part23[n * k + i] = 0.0;
	for (m=0;m<n;m++)
	  {
	    prod = 1.0;
	    for (j=b[k];j<b[k+1];j++)
	      prod *= 1.0 - fmax2(U[n * j + i], U[n * j + m]);
	    part23[n * k + i] += prod;
	  }
	part23[n * k + i] /= n;
      }
}

/************************************************************************/

void compute_part4(int n, int p, int *b, double *U, double *part4)
{
  int j,k,m,q;
  double prod;

  for (k=0;k<p;k++)
    {
      part4[k] = 0.0;
      for (m=0;m<n;m++)
	for (q=0;q<n;q++)
	  {
	    prod = 1.0;
	    for (j=b[k];j<b[k+1];j++)
	      prod *= 1.0 - fmax2(U[n * j + q], U[n * j + m]);
	    
	    part4[k] += prod;
	  }
      part4[k] /= n * n;
    }
}

/************************************************************************/

double M_A_n(int n, int p, int *b, double *U, double *part23, double *part4, 
	     int A)
{
  int i,j,k,l;
  double MAn, prod, part1;
  
  MAn = 0.0;
  for (i=0;i<n;i++)
    for (l=0;l<n;l++)
      {
	prod = 1.0;
	for (k=0;k<p;k++)
	  if (1<<k & A)
	    {
	      part1 = 1.0;
	      for (j=b[k];j<b[k+1];j++)
		part1 *= 1.0 - fmax2(U[n * j + i], U[n * j + l]);
	      
	      prod *= part1 - part23[n * k + i] - part23[n * k + l] 
		+ part4[k];
	    }
	MAn += prod;
      }

  return MAn/n;
}

/*****************************************************************************

  Bootstrap version of the Cramer-von Mises statistic used in the independence
  test 
  
******************************************************************************/

double boostrap_I_n(int n, int p, int *b, double *U, int *R)
{
  int i,j,k,l;
  double In, sum, sum2, prod, prod2, part1, part2, part3;

  /* first term */
  sum = 0.0;
  for (i=0;i<n;i++)
    for (l=0;l<n;l++)
      {
	prod = 1.0;
	for (k=0;k<p;k++)
	  for (j=b[k];j<b[k+1];j++)
	    prod *= 1.0 - fmax2(U[n * j + R[n * k + i]], U[n * j + R[n * k + l]]);
	sum += prod;
      }
  part1 = sum/n;
    
  /* second term */
  sum = 0.0;
  for (i=0;i<n;i++)
    {
      prod = 1.0;
      for (k=0;k<p;k++)
	{
	  sum2 = 0.0;
	  for (l=0;l<n;l++)
	    {
	      prod2 = 1.0;
	      for (j=b[k];j<b[k+1];j++)
		prod2 *= 1.0 - fmax2(U[n * j + R[n * k + i]], U[n * j + R[n * k + l]]);
	      sum2 += prod2;	      
	    }
	  prod *= sum2;
	}
      sum += prod;
    }
  part2 = 2.0 * sum / R_pow_di(n,p);
  

  /* third term */
  prod = 1.0;
  for (k=0;k<p;k++)
    {
      sum = 0.0;
      for (i=0;i<n;i++)
	for (l=0;l<n;l++)
	  { 
	    prod2 = 1.0;
	    for (j=b[k];j<b[k+1];j++)
	      prod2 *= 1.0 - fmax2(U[n * j + R[n * k + i]], U[n * j + R[n * k + l]]);
	    sum += prod2;
	  }
      prod *= sum;
    }
  part3 = prod / R_pow_di(n,2*p-1);

  In = part1 - part2 + part3;
 
  return In;
}

/*****************************************************************************/

double bootstrap_M_A_n(int n, int p, int *b, double *U, int *R, double *part23,
		       double *part4, int A)
{
  int i,j,k,l;
  double MAn, prod, part1;
  
  MAn = 0.0;
  for (i=0;i<n;i++)
    for (l=0;l<n;l++)
      {
	prod = 1.0;
	for (k=0;k<p;k++)
	  if (1<<k & A)
	    {
	      part1 = 1.0;
	      for (j=b[k];j<b[k+1];j++)
		part1 *= 1.0 - fmax2(U[n * j + R[n * k + i]], 
				     U[n * j + R[n * k + l]]);


	      prod *= part1 - part23[n * k + R[n * k + i]] 
		- part23[n * k + R[n * k + l]] + part4[k];
	    }
	MAn += prod;
      }

  return MAn/n;
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
	       double *MA0, double *I0, int *subset, char **subset_char)
{
  int i, j, k, sb = (int)sum_binom(*p,*m);
  int *R = (int *)R_alloc((*n) * (*p), sizeof(int));
  double *part23 = (double *)R_alloc((*n) * (*p), sizeof(double));
  double *part4 = (double *)R_alloc(*p, sizeof(double));
  
  /* partial power set in "natural" order */ 
  k_power_set(p, m, subset);

  /* convert partial power set to char for printing */
  k_power_set_char(p,m, subset, subset_char);
  
  /* Simulate the distribution of TA under independence */
  GetRNGstate();
  
  /* N repetitions */
  for (k=0;k<*N;k++) { 

    Rprintf("Bootstrap iteration %d\n",k+1);
    
    /* generate row selection within the blocks */
    for (j=0;j<*p;j++) 
      for (i=0;i<*n;i++)  
	R[(*n) * j + i] = (int)( unif_rand() * (*n) );
    
    /* compute part23 and part4 of the test statistics */
    compute_part23(*n, *p, b, U, part23);
    compute_part4(*n, *p, b, U, part4);

    /* for subsets i of cardinality greater than 1 and lower than m */
    for (i=*p+1;i<sb;i++) 
      MA0[k + (*N) * (i - *p - 1)] = bootstrap_M_A_n(*n, *p, b, U, R, part23,
						     part4, subset[i]);
    /* global statistic */ 
    I0[k] = boostrap_I_n(*n, *p, b, U, R); 
  }
  PutRNGstate();
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
  double *fisher0 = (double *)R_alloc(*N, sizeof(double));
  double *tippett0 = (double *)R_alloc(*N, sizeof(double));
  double *part23 = (double *)R_alloc((*n) * (*p), sizeof(double));
  double *part4 = (double *)R_alloc(*p, sizeof(double));
  double pvalue;

  /* compute W à la Fisher and à la Tippett from MA0*/ 
  for (k=0;k<*N;k++) {
    fisher0[k] = 0.0;
    tippett0[k] = 1.0;
    for (i=*p+1;i<sb;i++) {
      /* p-value */
      count = 1;
      for (j=0;j<*N;j++) 
	if (MA0[j + (*N) * (i - *p - 1)] > MA0[k + (*N) * (i - *p - 1)])
	  count ++;
      pvalue = (double)(count + 0.5)/(*N + 1.0);
      fisher0[k] -= 2*log(pvalue);
      tippett0[k] = fmin2(tippett0[k],pvalue);
    }
  }

  /* compute W from the current data */
  *fisher = 0.0;
  *tippett = 1.0;
  
  /* compute part4 of the test statistics */
  compute_part23(*n, *p, b, U, part23);
  compute_part4(*n, *p, b, U, part4);

  /* for subsets i of cardinality greater than 1 */
  for (i=*p+1;i<sb;i++) 
    {
      MA[i - *p - 1] = M_A_n(*n, *p, b, U, part23, part4, 
			     subset[i - *p - 1]);
      
      /* p-value */
      count = 1;
      for (k=0;k<*N;k++) 
	if (MA0[k + (*N) * (i - *p - 1)] > MA[i - *p - 1])
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

  /* compute In from the current data and the corresponding pvalue*/ 
  *I = I_n(*n, *p, b, U); 
  count = 1;
  for (k=0;k<*N;k++) 
    if (I0[k] > *I)
      count ++;
  *Ipval = (double)(count)/(*N + 1.0);
}


