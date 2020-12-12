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
 * Sample V01 ~ F01 with Laplace-Stieltjes transform
 * Used, for example, for sampling F01 for Joe and for sampling F01 for Frank.
 * Note: The caller of this function must use GetRNGstate() and PutRNGstate().
 *       The commented part also includes the algorithm of Hofert (2011)
 *       "Efficiently sampling nested Archimedean copulas", but does not improve
 *       run time.
 *
 * @param V0 parameter V0
 * @param theta_0 parameter theta0 in (0,infinity)
 * @param theta_1 parameter theta1 in [theta0, infinity)
 * @param p0 parameter 1-exp(-theta0)
 * @param p1 parameter 1-exp(-theta1)
 * @param gamma_1_a Gamma(1-alpha)
 * @param rej method switch: if V0*theta_0*p0^(V0-1) > rej a rejection
 *        from F01 of Joe is applied (otherwise, the sum is
 *        sampled via a logarithmic envelope for the summands)
 * @param approx largest number of summands before asymptotics is used
 * @return a random variate from F01
 * @author Marius Hofert, Martin Maechler
 */
double rF01Frank(double V0, double theta0, double theta1, double p0, double p1,
		 double gamma_1_a, double rej, int approx)
{
    double alpha = theta0 / theta1, iAlpha = (theta1-theta0)/theta1, U, V;
    if(V0*theta0*pow(p0,V0-1.) > rej) { /**< sample V01 via standard rejection from F01 for Joe */
	do {
	    U = unif_rand();
	    V = rF01Joe(V0, alpha,gamma_1_a, approx);
	} while(U > pow(p1, V));
    } else { /**< sample V01 as the V0-fold sum where the summands are sampled
		via rejection with a logarithmic envelope */
	double Ip = exp(-theta1);
	V = 0.;
	double X;
	/* if(theta0 <= theta1){ */
	for(int j=0; j < (int) V0; j++){ /**< sample V01 as a sum */
	    do {
		U = unif_rand();
		X = rLog(p1,Ip);
	    } while (U*(X-alpha) > 1./beta(X, iAlpha)); /**< X is now one summand of the sum of V01 */
	    V += X;
	}
	/* }else{ */
	/*     for(int j=0; j < (int) V0; j++){ /\**< sample V01 as a sum *\/ */
	/* 	do { */
	/* 	    U = unif_rand(); */
	/* 	    X = rSibuya(alpha,gamma_1_a); */
	/* 	} while (U > pow(p1,X-1.)); /\**< X is now one summand of the sum of V01*\/ */
	/* 	V += X; */
	/*     } */
	/* } */
    }
    return V;
}


/**
 * Generate a vector of variates V01 ~ F01 with Laplace-Stieltjes transform
 * ((1-(1-exp(-t)*(1-e^(-theta1)))^alpha)/(1-e^(-theta0)))^V0. Vectorized version
 * of rF01Frank.
 *
 * @param V01 vector of random variates from F01 (result)
 * @param V0 vector of random variates from F0
 * @param n length of the vectors V0 and V01
 * @param theta0 parameter theta0 in (0,infinity)
 * @param theta1 parameter theta1 in [theta0, infinity)
 * @param rej method switch, see rF01Frank
 * @param approx largest number of summands before asymptotics is used
 * @return none
 * @author Marius Hofert
 */
void rF01Frank_vec(double *V01, const double *V0, R_xlen_t n, double theta0,
		   double theta1, double rej, int approx){
    double p0 = -expm1(-theta0), p1 = - expm1(-theta1),
	iAlpha = (theta1 - theta0) / theta1, gamma_1_a = gammafn(iAlpha);
    GetRNGstate();

    for(R_xlen_t i=0; i < n; i++)
	V01[i] = rF01Frank(V0[i], theta0, theta1, p0, p1, gamma_1_a, rej, approx);

    PutRNGstate();
    return;
}


/**
 * Generate a vector of variates V01 ~ F01 with Laplace-Stieltjes transform
 * ((1-(1-exp(-t)*(1-e^(-theta1)))^alpha)/(1-e^(-theta0)))^V0. Bridge to R. Used,
 * for example, to draw several variates from rF01Frank.
 *
 * @param V0_ vector of random variates from F0
 * @param theta_0_ parameter theta0 in (0,infinity)
 * @param theta_1_ parameter theta1 in [theta0, infinity)
 * @param rej_ method switch, see rF01Frank
 * @param approx_ largest number of summands before asymptotics is used
 * @return vector of random variates V01
 * @author Martin Maechler, Marius Hofert
 */
SEXP rF01Frank_vec_c(SEXP V0_, SEXP theta_0_, SEXP theta_1_, SEXP rej_, SEXP approx_){
    double *V0 = REAL(V0_);
    R_xlen_t n = xlength(V0_);
    double theta_0 = asReal(theta_0_), theta_1 = asReal(theta_1_), rej = asReal(rej_);
    int approx = asInteger(approx_);
    SEXP res = PROTECT(allocVector(REALSXP, n));
    if(n >= 1) rF01Frank_vec(REAL(res), V0, n, theta_0, theta_1, rej, approx);
    UNPROTECT(1);
    return res;
}
