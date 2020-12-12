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

#include "nacopula.h"


/**
 * Sample V from a Sibuya(alpha) distribution with cdf F(n) = 1-1/(n*B(n,1-alpha)),
 * n in IN, with Laplace-Stieltjes transform 1-(1-exp(-t))^alpha via the
 * algorithm of Hofert (2011).
 * Note: The caller of this function must use GetRNGstate() and PutRNGstate().
 *
 * @param alpha parameter theta0/theta1 in (0,1]
 * @param gamma_1_a Gamma(1-alpha)
 * @return a random variate from F
 * @author Marius Hofert, Martin Maechler
 */
double rSibuya(double alpha, double gamma_1_a /**< == Gamma(1 - alpha) */){
    /**< FIXME(MM): (for alpha not too close to 1): re-express using 1-U */
    double U = unif_rand();
    if(U <= alpha)
	return 1.;
    else { /**< alpha < U < 1 */
	const double xMax = 1./DBL_EPSILON; // ==> floor(x) == ceil(x)  for x >= xMax
	double Ginv = pow((1-U)*gamma_1_a, -1./alpha), fGinv = floor(Ginv);
	if(Ginv > xMax) return fGinv; /* else */
	if(1-U < 1./(fGinv*beta(fGinv, 1.-alpha))) return ceil(Ginv); /* else */
	return fGinv;
    }
}


/**
 * Sample an n-fold sum of i.i.d. V ~ Sibuya(alpha). Sum-version of rSibuya.
 * Note: The caller of this function must use GetRNGstate() and PutRNGstate().
 *
 * @param alpha parameter theta0/theta1 in (0,1]
 * @param gamma_1_a Gamma(1-alpha)
 * @return n-fold sum of i.i.d. V ~ F
 * @author Marius Hofert
 */
double rSibuya_sum(R_xlen_t n, double alpha, double gamma_1_a /**< == Gamma(1 - alpha) */){
    double n_sum_V = 0.;
    for(R_xlen_t i = 0; i < n; i++)
	n_sum_V += rSibuya(alpha, gamma_1_a);
    return n_sum_V;
}


/**
 * Generate a vector of variates from a Sibuya(alpha) distribution.
 *
 * @param V vector of random variates from F (result)
 * @param n length of the vector V
 * @param alpha parameter theta0/theta1 in (0,1]
 * @return none
 * @author Marius Hofert, Martin Maechler
 */
void rSibuya_vec(double* V, R_xlen_t n, double alpha){
    if(n >= 1) {
	double gamma_1_a = gammafn(1.-alpha);
	GetRNGstate();

	for(R_xlen_t i=0; i < n; i++)
	    V[i] = rSibuya(alpha, gamma_1_a);

	PutRNGstate();
    }
    return;
}


/**
 * Generate a vector of variates from a Sibuya(alpha) distribution. Bridge to R.
 *
 * @param n sample size
 * @param alpha parameter theta0/theta1 in (0,1]
 * @return vector of random variates
 * @author Martin Maechler
 */
SEXP rSibuya_vec_c(SEXP n_, SEXP alpha_){
    R_xlen_t n;
#ifdef LONG_VECTOR_SUPPORT
    double dn = asReal(n_);
    if (ISNAN(dn) || dn < 0 || dn > R_XLEN_T_MAX)
	error(_("invalid 'n'"));
    n = (R_xlen_t) dn;
#else
    n = asInteger(n_);
    if (n == NA_INTEGER || n < 0)
	error(_("invalid 'n'"));
#endif
    double alpha = asReal(alpha_);
    SEXP res = PROTECT(allocVector(REALSXP, n));
    if(n >= 1) rSibuya_vec(REAL(res), n, alpha);
    UNPROTECT(1);
    return(res);
}
