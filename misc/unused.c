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

/**< lang5 will probably be defined in future versions of R */
#ifndef lang5
static SEXP lang5(SEXP s, SEXP t, SEXP u, SEXP v, SEXP w)
{
    PROTECT(s);
    s = LCONS(s, list4(t, u, v, w));
    UNPROTECT(1);
    return s;
}
#endif


/**
 * This is really not useful, as the time is spent in the call to R's rstable1()
 * Keep it as a didactical example of how to call R from C
 * @param V vector of variables V01 (result)
 * @param V0 vector of variates V0
 * @param alpha parameter theta0/theta1 in (0,1]
 * @param n number of variates V0, thus V01
 * @return none
 * @author Martin Maechler
*/
void C_retstable_R(double *V, const double V0[], double alpha, int n)
{
    double c0 = cos(M_PI_2 * alpha);
    SEXP one, alpha_R;
    one     = PROTECT(ScalarReal(1.));
    alpha_R = PROTECT(ScalarReal(alpha));

    GetRNGstate();
    for(int i = 0; i < n; i++) {
	/**< find m := optimal number of summands, using asymptotic formula */
	double m = fmax(1., floor(0.5 + V0[i]));
	/**
	 * Previously:
	 * double m;
	 * if(V0[i] <= 1)
	 *     m= 1.;
	 * else {
	 *     double fV= floor(V0[i]), cV= ceil(V0[i]);
	 *     if(pow(exp(-V0[i]), 1./cV - 1./fV) <= cV/fV)
	 * 	m= fV;
	 *     else
	 * 	m= cV;
	 * }
	*/
	/**< apply standard rejection m times */
	SEXP gamma_, rstableCall;
	PROTECT(gamma_ = ScalarReal(pow(c0*V0[i]/m, 1./alpha)));
	PROTECT(rstableCall = lang5(install("rstable1"), one, /**< alpha= */
			alpha_R, /**<beta= */ one, /**< gamma= */ gamma_));
	for(int k= 0; k < m; k++) {
	    double V01_k, U;
	    /**< apply standard rejection */
	    do {
		/**
		* sample from the distribution corresponding to the
		* Laplace-Stieltjes transform exp(-(V0/m)*t^alpha) via R
		*/
		V01_k = asReal(eval(rstableCall, R_GlobalEnv));
		U = unif_rand();

	    } while (U > exp(- V01_k));

	    V[i] += V01_k;
	}
	UNPROTECT(2);
    } /**< end for */
    PutRNGstate();
    UNPROTECT(1);
    return;
}

