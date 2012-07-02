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
 * @file   An.c
 * @author Ivan Kojadinovic
 * @date   2009-2012
 *
 * @brief Rank-based versions of the Pickands and CFG estimators
 *        of the Pickands dependence function. Bivariate and 
 *        multivariate versions.
 *
 */

#include <R.h>
#include <Rmath.h>

#include "An.h"


/**
 * Inverse of the rank-based version of the Pickands estimator
 * Bivariate version
 *
 * @param n sample size
 * @param S unit Fréchet pseudo-obs
 * @param T unit Fréchet pseudo-obs
 * @param t argument
 * @return value at t
 * @author Ivan Kojadinovic
 */
double biv_invAP(int n, double *S, double *T, double t) {
  double At = 0.0;
  if (t > 0.0 && t < 1.0) {
    for (int i = 0; i < n; i++) {
      double Sit = S[i] / (1.0 - t), 
	Tit = T[i] / t, 
	Eit = MIN(Sit, Tit);
      At += Eit;
    }
  }
  else if (t <= 0.0)  /* t set to 0 */
    for (int i = 0; i < n; i ++)
      At += S[i];
  else  /* t >= 1.0 */ /* t set to 1 */
    for (int i = 0; i < n; i++)
      At += T[i];
  return At / n;
}

/**
 * Rank-based version of the Pickands estimator interfaced in R
 * Bivariate version
 *
 * @param n sample size
 * @param S unit Fréchet pseudo-obs
 * @param T unit Fréchet pseudo-obs
 * @param t vector argument
 * @param m length of t and A
 * @param corrected if non zero, return the corrected version
 * @param A array of values at t
 * @author Ivan Kojadinovic
 */
void biv_AP(int *n, double *S, double *T, double *t, int *m,
		int *corrected, double *A) {
  if (*corrected) {
    double invA0 = biv_invAP(*n, S, T, 0.0);
    for (int i = 0; i < *m; i++)
      A[i] = 1.0 / (biv_invAP(*n, S, T, t[i])
		    - invA0 + 1.0);
  }
  else
    for (int i = 0; i < *m; i++)
      A[i] = 1.0 / biv_invAP(*n, S, T, t[i]);
}

/* Euler's constant */
#define EULER 0.5772156649015328606065120900824

/**
 * Log of the rank-based version of the CFG estimator
 * Bivariate version
 *
 * @param n sample size
 * @param S unit Fréchet pseudo-obs
 * @param T unit Fréchet pseudo-obs
 * @param t argument
 * @return value at t
 * @author Ivan Kojadinovic
 */
double biv_logACFG(int n, double *S, double *T, double t) {
  double At = 0.0;
  if (0. < t && t < 1.) {
    for (int i = 0; i < n ; i++) {
      double
	Sit = S[i] / (1.0 - t),
	Tit = T[i] / t,
	Eit = MIN(Sit, Tit);
      At += log(Eit);
    }
  }
  else if (t <= 0.0) /* t set to 0 */
    for (int i = 0; i < n; i++)
      At += log(S[i]);
  else /* t == 1.0 */ /* t set to 1 */
    for (int i = 0; i < n; i++)
      At += log(T[i]);

  return -EULER - At / n;
}

/**
 * Rank-based version of the CFG estimator interfaced in R
 * Bivariate version
 *
 * @param n sample size
 * @param S unit Fréchet pseudo-obs
 * @param T unit Fréchet pseudo-obs
 * @param t vector argument
 * @param m length of t and A
 * @param corrected if non zero, return the corrected version
 * @param A array of values at t
 * @author Ivan Kojadinovic
 */
void biv_ACFG(int *n, double *S, double *T, double *t, int *m,
	   int *corrected, double *A) {
  if (*corrected) {
    double
      logA0 = biv_logACFG(*n, S, T, 0.0);
    for (int i = 0; i < *m; i++)
      A[i] = exp(biv_logACFG(*n, S, T, t[i]) - logA0);
  }
  else
    for (int i = 0; i <* m; i++)
      A[i] = exp(biv_logACFG(*n, S, T, t[i]));
}



/**
 * Utility function: computes the xi at line l of w
 *
 * @param U pseudo-observations
 * @param n sample size
 * @param d dimension
 * @param w matrix with m lines and d columns
 * @param m number of lines of m
 * @param l line at which to compute the xi
 * @param xw the xi at line l of w
 * @author Ivan Kojadinovic
 */
void x_w(const double U[], int n, int d, const double w[],
	 int m, int l, double *xw) {
  int j;
  for (int i = 0; i < n; i++) {
    /// find first non null w[l + m * j]
    for (j = 0; j < d; j++)
      if (w[l + m * j] > 0.0) {
	xw[i] = -log(U[i + n * j]) / w[l + m * j];
	break;
      }
    /// find smallest term
    for (int k = j + 1; k < d; k++)
      if (w[l + m * k] > 0.0) {
	double y = -log(U[i + n * k]) / w[l + m * k];
	if (y < xw[i])
	  xw[i] = y;
      }
  }
}


/// Utility function: computes the inverse of An^P
double invAP(const double xw[], int n) {
  double invAw = 0.0;
  for (int i = 0; i < n; i++)
    invAw += xw[i];
  return invAw / n;
}

/// Utility function: computes the log of An^CFG
double logACFG(const double xw[], int n) {
  double logAw = 0.0;
  for (int i = 0; i < n; i++)
    logAw += log(xw[i]);
  return - logAw / n - EULER;
}

/**
 * Nonparametric estimators of the Pickands dependence function
 * interfaced in R; multivariate version
 * See the article 'Nonparametric estimation of multivariate
 * extreme-value copulas' by Gudendorf and Segers on the arXiv
 *
 * @param U pseudo-observations
 * @param n sample size
 * @param d dimension
 * @param w matrix with m lines and d columns
 * @param m number of lines of m
 * @param AP the m values of AP
 * @param ACFG the m values of ACFG
 * @param AHT the m values of AHT
 * @author Ivan Kojadinovic
 */
void mult_A(double *U, int *n, int *d, double *w, int *m,
	    double *AP, double *ACFG, double *AHT) {

  double *xw = Calloc(*n, double);
  double *xw0 = Calloc(*n, double);

  /// for corrections
  for (int i = 0; i < *n; i++)
    xw0[i] = log( (*n + 1.0) / (i + 1.0) );

  /// for every line of w
  for (int i = 0; i < *m; i++) {
    /// compute the xi
    x_w(U, *n, *d, w, *m, i, xw);

    /// for correction
    double invAP0 = invAP(xw0, *n); // for corrections
    double invAPi = invAP(xw, *n);

    /// the rank-based estimates at line i of w
    AP[i] = 1.0 / ( invAPi - invAP0 + 1.0);
    AHT[i] = invAP0 / invAPi;
    ACFG[i] = exp( logACFG(xw, *n) - logACFG(xw0, *n));
  }

  Free(xw);
  Free(xw0);
}
