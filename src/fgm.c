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
 * @file   fgm.c
 * @author Ivan Kojadinovic
 * @date   2007
 *
 * @brief  Functions for the FGM copula family
 *         Need to be properly tested
 *
 */

#include <R.h>
#include <Rmath.h>

#include "copula.h"

/**
 * Verify that the pdf of the FGM is positive at each vertex of [0,1]^p
 * Based on the polynomial used in the FGM copula
 *
 * @param p dimension
 * @param alpha parameters in natural order (d+1 zeros added in front)
 * @param valid is 1 if pdf positive at each vertex, 0 otherwise
 * @author Ivan Kojadinovic
 */
void validity_fgm(int *p, double *alpha, int *valid) {
  int i, j, k;
  double prod, val;
  double *alpha_bin = Calloc(1<<*p, double);
  int *subset = Calloc(1<<*p, int);

  k_power_set(p, p, subset);
  natural2binary(p, alpha, subset, alpha_bin);

  /* for each vertex of [0,1]^p */
  for (k=0;k<(1<<*p);k++) {
    val = 1.0;
    for (i=1;i<(1<<*p);i++)
      if (card(i) > 1) {
	prod = alpha_bin[i];
	for (j=0;j<*p;j++)
	  if ((1<<j) & i)
	    prod *= 1.0 - 2.0 * (((1<<j) & k) == (1<<j));
	val += prod;
      }
    if (val < 0) {
      *valid = 0;
      return;
    }
  }
  *valid = 1;

  Free(alpha_bin);
  Free(subset);
}

/**
 * Function to compute the conditional FGM copula
 * Based on the polynomial used in the FGM copula
 *
 * @param u array of standard uniform variates (length p)
 * @param p dimension
 * @param S a subset
 * @param alpha parameter values in binary order
 * @return a double (see function rfgm below)
 * @author Ivan Kojadinovic
 */
static
double B(double *u, int p, int S, double *alpha) {
  int i, j;
  double prod, b = 1.0;

  for (i=1;i<(1<<p);i++)
    if (card(i) > 1 && ((S & i) == i)) {
      prod = alpha[i];
      for (j=0;j<p;j++)
	if ((1<<j) & S)
	  prod *= 1.0 - 2.0 * u[j];
      b += prod;
    }
  return b;
}

/**
 * Function to compute the conditional FGM copula
 * Based on the polynomial used in the FGM copula
 *
 * @param u array of standard uniform variates (length p)
 * @param p dimension
 * @param S a subset
 * @param m (see function rfgm below)
 * @param alpha parameter values in binary order
 * @return a double (see function rfgm below)
 * @author Ivan Kojadinovic
 */
static
double A(double *u, int p, int S, int m, double *alpha) {
  int i, j;
  double prod, a = 0.0;

  for (i=1;i<(1<<p);i++)
    if ((S & i) == i) {
      prod = alpha[i + (1<<m)];
      for (j=0;j<p;j++)
	if ((1<<j) & S)
	  prod *= 1.0 - 2.0 * u[j];
      a += prod;
    }
  return a;
}

/**
 * Random generation for the FGM copula; see Nelsen
 * Reference: Armstong and Galli (2002), working paper
 * http://www.cerna.ensmp.fr/Documents/MA-AG-WPCopula.pdf
 *
 * @param p dimension
 * @param alpha parameters in natural order (d+1 zeros added in front)
 * @param n sample size
 * @param x the random variates
 * @author Ivan Kojadinovic
 */
void rfgm(int *p, double *alpha, int *n, double *x) {
  int i,j,S;
  double a,b;
  double *alpha_bin = Calloc(1<<*p, double);
  int *subset = Calloc(1<<*p, int);

  k_power_set(p, p, subset);
  natural2binary(p, alpha, subset, alpha_bin);

  GetRNGstate();

  for (i=0;i<*n;i++) {
    x[i * (*p)] = unif_rand();
    S = (1<<0);
    for (j=1;j<*p;j++) {
      b = B(x + i * (*p), *p, S, alpha_bin);
      a = A(x + i * (*p), *p, S, j, alpha_bin);
      if (fabs(a) < 1e-16)
	x[i * (*p) + j] = unif_rand();
      else
	x[i * (*p) + j] =
	  (a+b - sqrt((a+b)*(a+b) - 4.*a*b * unif_rand())) / (2.*a) ;
      S += (1<<j);
    }
  }

  PutRNGstate();

  Free(alpha_bin);
  Free(subset);
}


