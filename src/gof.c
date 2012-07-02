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
 * @file   gof.c
 * @author Ivan Kojadinovic
 * @date   2008
 *
 * @brief  Goodness-of-fit tests for copulas and EV copulas
 *         based on the parametric bootstrap and on a multiplier
 *         bootstrap
 *
 */


#include <R.h>
#include <Rmath.h>

#include "An.h"
#include "gof.h"
#include "empcop.h"

/**
 * Cramer-von Mises test statistic
 *
 * @param n sample size
 * @param p dimension
 * @param U pseudo-observations
 * @param Ctheta values of the fitted copula at U
 * @param stat value of the test statistic
 * @author Ivan Kojadinovic
 */
void cramer_vonMises(int *n, int *p, double *U, double *Ctheta,
		     double *stat) {
  double s = 0.;
  for(int k = 0; k < *n; k++) {
    double diff = multCn(U, *n, *p, U, *n, k, 0.0) - Ctheta[k];
    s += diff * diff;
  }
  *stat = s;
}

/**
 * Cramer-von Mises test statistic (grid version)
 *
 * @param p dimension
 * @param U pseudo-observations
 * @param n sample size
 * @param V grid
 * @param m grid size
 * @param Ctheta values of the fitted copula at V
 * @param stat value of the test statistic
 * @author Ivan Kojadinovic
 */
void cramer_vonMises_grid(int *p, double *U, int *n, double *V, int *m,
			  double *Ctheta, double *stat) {
  double s = 0.;
  for (int k = 0; k < *m; k++) {
    double diff = multCn(U, *n, *p, V, *m, k, 0.0) - Ctheta[k];
    s +=  diff * diff;
  }
  *stat = s * (*n) / (*m);
}

/**
 * Multiplier bootstrap for the GoF testing
 *
 * @param p dimension
 * @param U pseudo-observations (n x p)
 * @param n sample size
 * @param G grid (g x p)
 * @param g grid size
 * @param influ influence matrix (g x n)
 * @param N number of multiplie replications
 * @param s0 N replications of the test statistic
 * @author Ivan Kojadinovic
 */
void multiplier(int *p, double *U, int *n, double *G, int *g,
		double *influ, int *N, double *s0)
{
  int i, j, k, l, ind;
  double *influ_mat = Calloc((*n) * (*g), double);
  double *random = Calloc(*n, double);
  double *v1 = Calloc(*p, double);
  double *v2 = Calloc(*p, double);
  double *der = Calloc(*p, double);
  double mean, process, invsqrtn = 1.0/sqrt(*n);

  /* influence matrix */
  for (j = 0; j < *g; j++) /* loop over the grid points */
    {
      /* derivatives wrt args */
      for (k = 0; k < *p; k++)
	{
	  v1[k] = G[j + k * (*g)];
	  v2[k] = v1[k];
	}
      for (k = 0; k < *p; k++)
	{
	  v1[k] += invsqrtn;
	  v2[k] -= invsqrtn;
	  der[k] = der_multCn(U, *n, *p, v1, v2, 2 * invsqrtn);
	  v1[k] -= invsqrtn;
	  v2[k] += invsqrtn;
	}

      for (i = 0; i < *n; i++) /* loop over the data */
	{
	  influ_mat[i + j * (*n)] = 0.0;
	  ind = 1;
	  for (k = 0; k < *p; k++)
	    {
	      ind *= (U[i + k * (*n)] <= G[j + k * (*g)]);
	      influ_mat[i + j * (*n)] -= der[k] * (U[i + k * (*n)] <= G[j + k * (*g)]);
	    }
	  influ_mat[i + j * (*n)] += ind; /* - influ[j + i * (*g)];*/
	  influ[j + i * (*g)] *= invsqrtn;
	  influ_mat[i + j * (*n)] *= invsqrtn;
	}
    }

  GetRNGstate();

  /* generate N approximate realizations */
  for (l=0;l<*N;l++)
    {
      /* generate n variates */
      mean = 0.0;
      for (i=0;i<*n;i++)
	{
	  random[i] = norm_rand(); /*(unif_rand() < 0.5) ? -1.0 : 1.0 ;*/
	  mean += random[i];
	}
      mean /= *n;

      /* realization number l */
      s0[l] = 0.0;
      for (j=0;j<*g;j++)
	{
	  process = 0.0;
	  for (i=0;i<*n;i++)
	    process += (random[i] - mean) * influ_mat[i + j * (*n)]
	      - random[i] * influ[j + i * (*g)];
	  s0[l] += process * process;
	}
      s0[l] /= *g;
    }

  PutRNGstate();

  Free(influ_mat);
  Free(random);
  Free(v1);
  Free(v2);
  Free(der);
}

/// Goodness-of-fit tests for extreme-value copulas

/**
 * Cramer-von Mises test statistic based on the Pickands estimator
 * used in the GOF for EV copulas
 * stat = n/m \sum_{i=0}^{m-1} [An(i/m) - A_{theta_n}(i/m)]^2
 *
 * @param n sample size
 * @param m grid size
 * @param S unit Fréchet pseudo-observations
 * @param T unit Fréchet pseudo-observations
 * @param Atheta values of the fitted A at the grid points
 * @param stat value of the test statistic
 * @author Ivan Kojadinovic
 */
void cramer_vonMises_Pickands(int *n, int *m, double *S,
			      double *T, double *Atheta,
			      double *stat) {
  double t, Ac, Au, dc, du,
    invA0 = biv_invAP(*n, S, T, 0.0);

  stat[0] = 0.0; stat[1] = 0.0;
  for (int i = 0; i < *m; i++) {
      t = (double)i/(double)(*m);
      Au = biv_invAP(*n, S, T, t);
      Ac = Au - invA0 + 1.0; // correction
      du = 1 / Au - Atheta[i];
      dc = 1 / Ac - Atheta[i];
      stat[0] += dc * dc;
      stat[1] += du * du;
    }
  stat[0] = stat[0] * (double)(*n) / (double)(*m);
  stat[1] = stat[1] * (double)(*n) / (double)(*m);
}

/**
 * Cramer-von Mises test statistic based on the CFG estimator
 * used in the GOF for EV copulas
 * stat = n/m \sum_{i=0}^{m-1} [An(i/m) - A_{theta_n}(i/m)]^2
 *
 * @param n sample size
 * @param m grid size
 * @param S unit Fréchet pseudo-observations
 * @param T unit Fréchet pseudo-observations
 * @param Atheta values of the fitted A at the grid points
 * @param stat value of the test statistic
 * @author Ivan Kojadinovic
 */
void cramer_vonMises_CFG(int *n, int *m, double *S,
			 double *T, double *Atheta,
			 double *stat) {
  double t, Au, Ac, dc, du,
    logA0 = biv_logACFG(*n, S, T, 0.0);

  stat[0] = 0.0; stat[1] = 0.0;
  for (int i = 0; i < *m; i++) {
      t = (double) i / (double) (*m);
      Au = biv_logACFG(*n, S, T, t);
      Ac = Au - logA0; // correction
      dc = exp(Ac) - Atheta[i];
      du = exp(Au) - Atheta[i];
      stat[0] += dc * dc;
      stat[1] += du * du;
    }
  stat[0] = stat[0] * (double)(*n)/(double)(*m);
  stat[1] = stat[1] * (double)(*n)/(double)(*m);
}

/**
 * Wrapper for the Cramer-von Mises test statistics
 * used in the GOF for EV copulas
 * stat = n/m \sum_{i=0}^{m-1} [An(i/m) - A_{theta_n}(i/m)]^2
 *
 * @param n sample size
 * @param m grid size
 * @param S unit Fréchet pseudo-observations
 * @param T unit Fréchet pseudo-observations
 * @param Atheta values of the fitted A at the grid points
 * @param stat value of the test statistic
 * @param CFG if > 0 then CFG, else Pickands
 * @author Ivan Kojadinovic
 */
void cramer_vonMises_Afun(int *n, int *m, double *S,
			  double *T, double *Atheta,
			  double *stat, int *CFG)
{
  if (*CFG)
    cramer_vonMises_CFG(n, m, S, T, Atheta, stat);
  else
    cramer_vonMises_Pickands(n, m, S, T, Atheta, stat);
}
