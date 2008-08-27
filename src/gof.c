/*#################################################################################
##
##   R package Copula by Jun Yan and Ivan Kojadinovic Copyright (C) 2008
##
##   This file is part of the R package copula.
##
##   The R package copula is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package copula is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with the R package copula. If not, see <http://www.gnu.org/licenses/>.
##
#################################################################################*/


/*****************************************************************************

  Goodness-of-fit for copulas 

  Ivan Kojadinovic, April 2008

*****************************************************************************/

#include <R.h>
#include <Rmath.h>
#include "set.utils.h"

/***********************************************************************
  
  Computes the empirical copula
  U contains the pseudo-obs that define the empirical copula
  V[k + m * j], j=1...p, is the value at which the empirical copula 
  is to be computed
  m is the number of lines of V

***********************************************************************/

double empcop(int n, int p, double *U, double *V, int m, int k)
{
  int i,j, ind;
  double ec = 0.0;
  
  for (i=0;i<n;i++)
    {
      ind = 1;
      for (j=0;j<p;j++)
	ind *= (U[i + n * j] <= V[k + m * j]);
      ec += (double)ind;
    }
  return ec/(double)n;
}

/***********************************************************************
  
  Computes the Cramer-von Mises test statistic
  U contains the pseudo-obs that define the empirical copula
  FU are observations simulated from a copula fitted from U
  m is the size of the sample used to approximate the test statistic

  Two versions: with empcop approximation, without
 
***********************************************************************/

void cramer_vonMises_approx(int *n, int *p, double *U, double *FU, int *m, 
			    double *stat)
{
  int i;
  double s = 0.0, diff;
  
  for (i=0;i<*n;i++)
    {
      diff = empcop(*n,*p,U,U,*n,i) - empcop(*m,*p,FU,U,*n,i);
      s += diff * diff;
    }
  
  *stat = s;
}

void cramer_vonMises(int *n, int *p, double *U, double *Ctheta, 
		     double *stat)
{
  int i;
  double s = 0.0, diff;
  
  for (i=0;i<*n;i++)
    {
      diff = empcop(*n,*p,U,U,*n,i) - Ctheta[i];
      s += diff * diff;
    }
  
  *stat = s;
}

void cramer_vonMises_2(int *p, double *U, int *n, double *V, int *m, 
		       double *Ctheta, double *stat)
{
  int i;
  double s = 0.0, diff;
  
  for (i=0;i<*m;i++)
    {
      diff = empcop(*n,*p,U,V,*m,i) - Ctheta[i];
      s +=  diff * diff;
    }
  
  *stat = s * (*n) / (*m);
}

/***********************************************************************
  
  Computes the Cramer-von Mises test statistics from the 
  Mobius decomposition of goodness-of-fit process
  U contains the pseudo-obs that define the empirical copula
  FU are observations simulated from a copula fitted from U
  m is the size of the sample used to approximate the test statistic
  nsubsets is the number of subset statistic that will be computed
  subset are the subsets in "natural order" (pairs, triplets, etc)

  Two versions: with empcop approximation, without
 
***********************************************************************/

double Mobius_empcop(int n, int p, int A, double *U, double *V, 
		     int m, int k)
{
  int i,j, ind;
  double ec = 0.0;
  
  for (i=0;i<n;i++)
    {
      ind = 1;
      for (j=0;j<p;j++)
	if (1<<j & A)
	  ind *= (U[i + n * j] <= V[k + m * j]);
      ec += (double)ind;
    }
  return ec/(double)n;
}

void Mobius_cramer_vonMises_approx(int *n, int *p, int *nsubsets, 
				   int *subset, double *U, double *FU, 
				   int *m, double *stat)
{
  int i,j,k;
  double *diff = Calloc((*n) * (*nsubsets), double);
  double sum;
  
  /* compute the differences */
  for (i=0;i<*n;i++)
      for (j=0;j<*nsubsets;j++) 
	diff[i + j * (*n)] = Mobius_empcop(*n,*p,subset[j],U,U,*n,i) 
	  - Mobius_empcop(*m,*p,subset[j],FU,U,*n,i);


  /* compute the Mobius statistics */
  for (k=0;k<*nsubsets;k++)
    {
      stat[k] = 0.0;
      for (i=0;i<*n;i++)
	{
	  sum = 0.0;
	  for (j=0;j<*nsubsets;j++)
	    if ((subset[j] & subset[k]) == subset[j])
	      /* can be improved */
	      sum += R_pow_di(-1,card(subset[j])) * diff[i + j * (*n)];
	  stat[k] += sum * sum; 
	}
    }

  /* compute diff for the global statistic if necessary */
  if (*nsubsets != ((int)R_pow_di(2.0,*p) - *p - 1)) 
    for (i=0;i<*n;i++) 
      diff[i + (*nsubsets - 1) * (*n)] = empcop(*n,*p,U,U,*n,i) - empcop(*m,*p,FU,U,*n,i);

  /* compute the global statistic */
  stat[*nsubsets] = 0.0;
  for (i=0;i<*n;i++)
    stat[*nsubsets] += diff[i + (*nsubsets - 1) * (*n)] * diff[i + (*nsubsets - 1) * (*n)]; 
  
  Free(diff); 
}

void Mobius_cramer_vonMises(int *n, int *p, int *nsubsets, int *subset, 
			    double *U, double *Ctheta, double *stat)
{
  int i,j,k;
  double *diff = Calloc((*n) * (*nsubsets), double);
  double sum;
  
  /* compute the differences */
  for (i=0;i<*n;i++)
    for (j=0;j<*nsubsets;j++) 
      diff[i + j * (*n)] = Mobius_empcop(*n,*p,subset[j],U,U,*n,i) 
	- Ctheta[i + j * (*n)];

 
  /* compute the Mobius statistics */
  for (k=0;k<*nsubsets;k++)
    {
      stat[k] = 0.0;
      for (i=0;i<*n;i++)
	{
	  sum = 0.0;
	  for (j=0;j<*nsubsets;j++)
	    if ((subset[j] & subset[k]) == subset[j])
	      /* can be improved */
	      sum += R_pow_di(-1,card(subset[j])) * diff[i + j * (*n)];
	  stat[k] += sum * sum; 
	}
    }
  
  /* compute diff for the global statistic if necessary */
  if (*nsubsets != ((int)R_pow_di(2.0,*p) - *p - 1)) 
    for (i=0;i<*n;i++) 
      diff[i + (*nsubsets - 1) * (*n)] = empcop(*n,*p,U,U,*n,i) - Ctheta[i + (*nsubsets - 1) * (*n)];
  
  /* compute the global statistic */
  stat[*nsubsets] = 0.0;
  for (i=0;i<*n;i++)
    stat[*nsubsets] += diff[i + (*nsubsets - 1) * (*n)] * diff[i + (*nsubsets - 1) * (*n)]; 
    
  Free(diff); 
}

/***********************************************************************
  
  Computes pvalues from Mobius statistics 
  stat is N + 1 by nsubsets + 1 array
  the last line contains the observed statistics
  the last column the global statistics (Genest et al.)
  the first N lines the statistics under H0
  pval array of nsubsets + 3; last three values 
  are Fisher, Tippett and global pvalues 
 
***********************************************************************/

void Mobius_pvalues(int *nsubsets, int *N, double *stat, double *pval)
{
  int i,j,k,count;
  double *fisher = Calloc(*N+1, double);
  double *tippett = Calloc(*N+1, double);
  double pvalue;

  /* compute W à la Fisher and à la Tippett from stat */ 
  for (k=0;k<=*N;k++) 
    {
      fisher[k] = 0.0;
      tippett[k] = 1.0;
      for (i=0;i<=*nsubsets;i++) 
	{
	  /* p-value */
	  count = 0;
	  for (j=0;j<*N;j++) 
	    if (stat[j + (*N+1) * i] >= stat[k + (*N+1) * i])
	      count ++;
	  pvalue = (double)(count + 0.5)/(*N + 1.0);
	  /* save the pvalues of the observed statistics */
	  if (k == *N)
	    pval[i] = pvalue;
	  if (i < *nsubsets)
	    {
	      fisher[k] -= 2*log(pvalue);
	      tippett[k] = fmin2(tippett[k],pvalue);
	    }
	}
    }
  
  /* p-values of the Fisher and Tippett statistics */
  count = 0;
  for (k=0;k<*N;k++) 
    if (fisher[k] >= fisher[*N])
      count ++;
  pval[*nsubsets+1] = (double)(count + 0.5)/(*N + 1.0);

  count = 0;
  for (k=0;k<*N;k++) 
    if (tippett[k] <= tippett[*N])
      count ++;
  pval[*nsubsets+2] = (double)(count + 0.5)/(*N + 1.0);

  Free(fisher);
  Free(tippett);

}

/***********************************************************************
  
  Influence matrices for multiplier CLT / OLD VERSION
  u0 is m x p
  u is n x p
  Ctheta is n 
  der is n x (p+1)
  s0 is N
 
***********************************************************************/

void simulate(int *p, double *u0, int *m, double *u, int *n, double *Ctheta, 
	      double *der, double *influ, int *N, double *s0) 
{
  int i, j, k, l, ind;
  double *influ_mat = Calloc((*m) * (*n), double);
  double *random = Calloc(*m, double);
  double process;
  /*double mean;*/

  /* centered influence matrix for the nonparametric part */
  /* centering commented out */
  for (j = 0; j < *n; j++) { /* the jth column */
    /*mean = 0.0;*/
    for (i = 0; i < *m; i++) {
	influ_mat[i + j * (*m)] = 0.0;
	ind = 1;
	for (k = 0; k < *p; k++) {
	  ind *= (u0[i + k * (*m)] <= u[j + k * (*n)]);
	  influ_mat[i + j * (*m)] -= der[j + k * (*n)] 
	    * ((u0[i + k * (*m)] <= u[j + k * (*n)]) - u[j + k * (*n)]);
	}
	influ_mat[i + j * (*m)] += ind - Ctheta[j];
	/*mean += influ_mat[i + j * (*m)];*/
    }
    /*mean /= (*m); 
    for (i = 0; i < *m; i++)
    influ_mat[i + j * (*m)] -= mean;*/
  }

  /* combine with the parametric influence which may already be centered in R */
  for (i = 0; i < *m; i++) {
    for (j = 0; j < *n; j++) {
      influ_mat[i + j * (*m)] -= influ[i] * der[j + (*p) * (*n)];
      influ_mat[i + j * (*m)] /= sqrt(*m);
    }
  }

  GetRNGstate();

  /* generate N approximate realizations */
  for (l=0;l<*N;l++)
    {
      /*if ((l+1) % 100 == 0)
	Rprintf("Simulation iteration %d\n",l+1);*/

      /* generate m variates */
      for (i=0;i<*m;i++)
	random[i] = norm_rand(); /*(unif_rand() < 0.5) ? -1.0 : 1.0 ;*/
      
      /* realization number l */
      s0[l] = 0.0;
      for (j=0;j<*n;j++)
	{ 
	  process = 0.0;
	  for (i=0;i<*m;i++)
	    process += random[i] * influ_mat[i + j * (*m)];

	  s0[l] += process * process;
	}
      s0[l] /= *n; 
    }

  PutRNGstate();

  Free(influ_mat);
  Free(random);
}


/***********************************************************************

  Influence matrice and simulation (multiplier CLT)
  u0 is m x p
  u is n x p
  Ctheta is n 
  der is n x p
  influ is n x m
  s0 is N
 
***********************************************************************/


void multiplier(int *p, double *u0, int *m, double *u, int *n, double *Ctheta, 
		double *der, double *influ, int *N, double *s0) 
{
  int i, j, k, l, ind;
  double *influ_mat = Calloc((*m) * (*n), double);
  double *random = Calloc(*m, double);
  double process;

  /* influence matrix */
  for (j = 0; j < *n; j++) { /* the jth column */
    for (i = 0; i < *m; i++) {
	influ_mat[i + j * (*m)] = 0.0;
	ind = 1;
	for (k = 0; k < *p; k++) {
	  ind *= (u0[i + k * (*m)] <= u[j + k * (*n)]);
	  influ_mat[i + j * (*m)] -= der[j + k * (*n)] 
	    * ((u0[i + k * (*m)] <= u[j + k * (*n)]) - u[j + k * (*n)]);
	}
	influ_mat[i + j * (*m)] += ind - Ctheta[j] - influ[j + i * (*n)]; /* influ transposed ! */
	influ_mat[i + j * (*m)] /= sqrt(*m);
    }
  }

  GetRNGstate();

  /* generate N approximate realizations */
  for (l=0;l<*N;l++)
    {
      /*if ((l+1) % 100 == 0)
	Rprintf("Simulation iteration %d\n",l+1);*/

      /* generate m variates */
      for (i=0;i<*m;i++)
	random[i] = norm_rand(); /*(unif_rand() < 0.5) ? -1.0 : 1.0 ;*/
      
      /* realization number l */
      s0[l] = 0.0;
      for (j=0;j<*n;j++)
	{ 
	  process = 0.0;
	  for (i=0;i<*m;i++)
	    process += random[i] * influ_mat[i + j * (*m)];

	  s0[l] += process * process;
	}
      s0[l] /= *n; 
    }

  PutRNGstate();

  Free(influ_mat);
  Free(random);
}

/***********************************************************************/
