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
 * @file   empcop.c
 * @author Ivan Kojadinovic
 * @date   Sat Jul 21 17:25:20 2012
 *
 * @brief  Bivariate and multivariate versions of the empirical copula
 *         and related estimators of the partial derivatives
 *
 */

#include <Rmath.h>

/**
 * Bivariate empirical copula; used by exchTest and evTestA
 *
 * @param U pseudo-obs
 * @param V pseudo-obs
 * @param n sample size
 * @param u value at which to evalute the empirical copula
 * @param v value at which to evalute the empirical copula
 * @return value of the empirical copula at (u,v)
 * @author Ivan Kojadinovic
 */
double bivCn(const double U[], const double V[], int n, double u, double v) {
  double sumind = 0.0;
  for (int i = 0; i < n; i++)
      sumind += (U[i] <= u) * (V[i] <= v);
  return sumind / (double)n;
}

/**
 * Estimator of the first partial derivative of the unknown copula
 * See CJS paper 2011 -- used by exchTest and evTestA
 *
 * @param U pseudo-obs
 * @param V pseudo-obs
 * @param n sample size
 * @param u value at which to evalute the empirical derivative
 * @param v value at which to evalute the empirical derivative
 * @return value of the empirical derivative at (u,v)
 * @author Ivan Kojadinovic
 */
double der1bivCn(const double U[], const double V[], int n, double u, double v) {
  double invsqrtn = 1.0 / sqrt(n);
  if (u < invsqrtn)
    u = invsqrtn;
  else if (u > 1.0 - invsqrtn)
    u = 1.0 - invsqrtn;
  return (bivCn(U, V, n, u + invsqrtn, v) - bivCn(U, V, n, u - invsqrtn, v))
    / (2.0 * invsqrtn);
}

/**
 * Estimator of the second partial derivative of the unknown copula
 * See CJS paper 2011 -- used by exchTest and evTestA
 *
 * @param U pseudo-obs
 * @param V pseudo-obs
 * @param n sample size
 * @param u value at which to evalute the empirical derivative
 * @param v value at which to evalute the empirical derivative
 * @return value of the empirical derivative at (u,v)
 * @author Ivan Kojadinovic
 */
double der2bivCn(const double U[], const double V[], int n, double u, double v) {
  double invsqrtn = 1.0 / sqrt(n);
  if (v < invsqrtn)
    v = invsqrtn;
  else if (v > 1.0 - invsqrtn)
    v = 1.0 - invsqrtn;
  return (bivCn(U, V, n, u, v + invsqrtn) - bivCn(U, V, n, u, v - invsqrtn))
    / (2.0 * invsqrtn);
}

/**
 * Multivariate empirical copula
 *
 * @param U pseudo-observations
 * @param n sample size
 * @param p dimension of the pseudo-observations
 * @param V is vector representing a matrix of dimension m x p
 * @param m "number of lines" of V
 * @param k "line" of V at which to compute the empirical copula
 * @param o offset (usually 0.0)
 * @return the value of the empirical copula at V[k + m * j], j=1...p
 * @author Ivan Kojadinovic
 */
double multCn(const double U[], int n, int p, const double V[], int m, int k, double o) {
    double sumind = 0.0;
    for (int i = 0; i < n; i++) {
	int ind = 1;
	for (int j = 0; j < p; j++)
	    ind *= (U[i + n * j] <= V[k + m * j]);
	sumind += (double)ind;
    }
    return sumind / (n + o);
}

/**
 * Estimator of the derivative of the unknown copula by finite differences
 * See the CJS 2011 paper -- used by evTestC and by the multiplier gof tests
 *
 * @param U pseudo-observations
 * @param n sample size
 * @param p dimension of the pseudo-observations
 * @param u vector of length p
 * @param v vector of length p
 * @param denom denominator (of the order 2/sqrt(n))
 * @return the value of the empirical derivative
 * @author Ivan Kojadinovic
 */
double der_multCn(const double U[], int n, int p,
		  const double u[], const double v[], double denom) {
    return (multCn(U, n, p, u, 1, 0, 0.0) - multCn(U, n, p, v, 1, 0, 0.0)) / denom;
}

/** 
 * Computes the empirical copula at all the lines of V 
 * Called in R by .C
 * 
 * @param U pseudo-observations 
 * @param n sample size
 * @param p dimension of the pseudo-observations 
 * @param V is vector representing a matrix of dimension m x p 
 * @param m "number of lines" of V
 * @param ec values of the empirical copula at the lines of V 
 * @author Ivan Kojadinovic
 */
void RmultCn(double *U, int *n, int *p, double *V, int *m, double *ec) {
   for (int i = 0; i < *m; i++)
     ec[i] = multCn(U, *n, *p, V, *m, i, 0.0);
}
