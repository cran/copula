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
 * Sample a Log(p) distribution with the algorithm "LK" of Kemp (1981).
 * Note: The caller of this function must use GetRNGstate() and PutRNGstate().
 *
 * @param p in (0,1)
 * @param Ip = 1 - p_ (possibly more accurate)
 * @return a random variate from Log(p)
 * @author Marius Hofert, Martin Maechler
 */
double rLog(double p, double Ip) {
    if(p <= 0. ||  p > 1.) {
	error("rLog(): p must be inside (0,1)");
	return -1.; /**< -Wall */
    }
    else if(Ip <= 0. || Ip >= 1.) {
	error("rLog(): Ip must be inside (0,1)");
	return -1.; /**< -Wall */
    }
    else {
	double U=unif_rand();
	if(U > p) {
	    return 1.;
	}
	else {
	    double Q, logQ;
	    if(p <= 0.5) {
		Q = - expm1(log1p(- p) * unif_rand()); /* = 1-(1-p)^unif */
		/**
		 * == 1. - exp(log1p(- p) * unif_rand())
		 * == 1. - pow(1. - p, unif_rand())
		 */
		logQ = log(Q);
	    } else { // p > 0.5  <==> Ip < 0.5
		double iQ = pow(Ip, unif_rand()); /* = (1-p)^unif */
		Q = 1. - iQ;
		logQ = log1p(-iQ);
	    }
	    return(U < Q*Q
		   ? floor(1. + log(U)/logQ)
		   : ((U > Q) ? 1. : 2.));
	}
    }
}


/**
 * Generate a vector of variates from a Log(p) distribution with the algorithm
 * "LK" of Kemp (1981).
 *
 * @param n_ sample size
 * @param p_ parameter p in (0,1)
 * @param Ip_ = 1 - p_ (possibly more accurate)
 * @return vector of random variates from Log(p)
 * @author Martin Maechler
 */
SEXP rLog_vec_c(SEXP n_, SEXP p_, SEXP Ip_) {
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
    double p = asReal(p_), Ip = asReal(Ip_);
    SEXP res = PROTECT(allocVector(REALSXP, n));
    double* X = REAL(res);

    GetRNGstate();

    for(R_xlen_t i=0; i < n; i++)
	X[i] = rLog(p, Ip);

    PutRNGstate();
    UNPROTECT(1);
    return res;
}
