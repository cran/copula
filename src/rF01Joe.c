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
 * Sample V01 ~ F01 with Laplace-Stieltjes transform ((1-(1-exp(-t))^alpha))^V0
 * Used, for example, for sampling F01 for Joe and for sampling F01 for Frank.
 * Note: The caller of this function must use GetRNGstate() and PutRNGstate().
 *
 * @param V0 parameter V0
 * @param alpha parameter theta0/theta1 in (0,1]
 * @param gamma_1_a Gamma(1-alpha)
 * @param approx largest number of summands before asymptotics is used
 * @return a random variate from F01
 * @author Marius Hofert, Martin Maechler
 */
double rF01Joe(double V0, double alpha, double gamma_1_a /**< == Gamma(1 - alpha) */,
	       int approx){
    if(V0 > approx) /**< approximation */
	return pow(V0,1./alpha)*rstable0(alpha); /**< rstable0() in retstable.c; */
    /* generates S(alpha, 1, (cos(alpha*pi/2))^{1/alpha}, I_{alpha == 1}; 1) */
    else /**< sample sum */
	return rSibuya_sum((int) V0, alpha, gamma_1_a);
}


/**
 * Generate a vector of variates V01 ~ F01 with Laplace-Stieltjes transform
 * ((1-(1-exp(-t))^alpha))^V0. Vectorized version of rF01Joe. Used, for example, to draw
 * several variates from rF01Joe.
 *
 * @param V01 vector of random variates from F01 (result)
 * @param V0 vector of random variates from F0
 * @param n length of the vector V0
 * @param alpha parameter theta0 in (0,1]
 * @param approx largest number of summands before asymptotics is used
 * @return none
 * @author Marius Hofert
 */
void rF01Joe_vec(double* V01, const double *V0, R_xlen_t n, double alpha, double approx) {
    double gamma_1_a = gammafn(1. - alpha);
    GetRNGstate();

    for(R_xlen_t i=0; i < n; i++)
	V01[i] = rF01Joe(V0[i], alpha, gamma_1_a, (int)approx);

    PutRNGstate();
    return;
}


/**
 * Generate a vector of variates V01 ~ F01 with Laplace-Stieltjes transform
 * ((1-(1-exp(-t))^alpha))^V0. Bridge to R. Used, for example, to draw several variates
 * from rF01Joe.
 *
 * @param V0_ vector of random variates from F0
 * @param alpha_ parameter alpha = theta0/theta1 in (0,1]
 * @param approx_ largest number of summands before asymptotics is used
 * @return vector of random variates V01
 * @author Marius Hofert
 */
SEXP rF01Joe_vec_c(SEXP V0_, SEXP alpha_, SEXP approx_){
    double *V0 = REAL(V0_);
    R_xlen_t n = xlength(V0_);
    double alpha = asReal(alpha_);
    int approx = asInteger(approx_);
    SEXP res = PROTECT(allocVector(REALSXP, n));
    if(n >= 1) rF01Joe_vec(REAL(res), V0, n, alpha, approx);
    UNPROTECT(1);
    return res;
}

