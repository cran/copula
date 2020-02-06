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
 * Wrapper of the beta df for the empirical beta copula
 * See Segers, Sibuya, Tsukahara, JMVA 2016+
 * u assumed to be rank/(n+1)
 **/
double emp_beta_cop(double u, double v, int n) {
    return pbeta(v, u * (n + 1), (n + 1) * (1.0 - u), TRUE, FALSE);
}

/**
 * For the empirical checkerboard (multilinear extension) copula
 * See Segers, Sibuya, Tsukahara, JMVA 2016+
 * u assumed to be rank/(n+1)
 **/
double emp_mult_lin_cop(double u, double v, int n) {
    return fmin2(fmax2(n * v - (n + 1) * u + 1.0 , 0.0), 1.0);
}

/* For the classical empirical copula */
double emp_cop(double u, double v, int n) {
    if (u <= v)
	return 1.0;
    else
	return 0.0;
}

/**
 * Multivariate empirical copula, empirical beta copula or
 * empirical checkerboard copula
 * See Segers, Sibuya, Tsukahara, JMVA 2016+
 *
 * @param U pseudo-observations (rank/(n+1))
 * @param n sample size
 * @param p dimension of the pseudo-observations
 * @param V is vector representing a matrix of dimension m x p
 * @param m "number of lines" of V
 * @param k "line" of V at which to compute the empirical ... copula
 * @param o offset (usually 0.0)
 * @param f function determing which version is computed
 * @return the value of the empirical copula at V[k + m * j], j=1...p
 * @author Ivan Kojadinovic
 */
double Cn_f(const double U[], int n, int p, const double V[], int m, int k,
	      double offset, double (*f)(double, double, int)) {
    double sumprod = 0.0;
    for (int i = 0; i < n; i++) {
	double prod = 1.0;
	for (int j = 0; j < p; j++)
	    prod *= f(U[i + n * j], V[k + m * j], n);
	sumprod += prod;
    }
    return sumprod / (n + offset);
}

/* The classical multivariate empirical copula */
double multCn(const double U[], int n, int p, const double V[], int m, int k,
	      double offset) {
    return Cn_f(U, n, p, V, m, k, offset, emp_cop);
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
 * @param offset offset to be used for scaling
 * @param type: 1 = empirical, 2 = beta, 3 = checkerboard
 * @author Ivan Kojadinovic
 */
void Cn_C(double *U, int *n, int *p, double *V, int *m, double *ec, double *offset,
	  int *type) { // called via .C() from ../R/empCopula.R
    switch(*type) {

    case 2 : /* empirical beta copula */
	for (int i = 0; i < *m; i++)
	    ec[i] = Cn_f(U, *n, *p, V, *m, i, *offset, emp_beta_cop);
	return;

    case 3 : /* empirical checkerboard copula */
	for (int i = 0; i < *m; i++)
	    ec[i] = Cn_f(U, *n, *p, V, *m, i, *offset, emp_mult_lin_cop);
	return;

    default : /* classical empirical copula */
	for (int i = 0; i < *m; i++)
	    ec[i] = multCn(U, *n, *p, V, *m, i, *offset);
    }
}
