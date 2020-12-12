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
 * Sample S_0 ~ S(alpha, 1, (cos(alpha*pi/2))^{1/alpha}, I_{alpha == 1}; 1)
 * with Laplace-Stieltjes transform exp(-t^alpha) via the algorithm as described
 * in Chambers, Mallows, Stuck (1976) (S_0 = S(alpha,1) in this reference), see
 * Nolan's book for the parametrization.
 *
 * @param alpha parameter in (0,1]
 * @return S_0
 * @author Marius Hofert, Martin Maechler
 */
double rstable0(double alpha){
	/**< S(1, 1, 0, 1; 1) with Laplace-Stieltjes transform exp(-t) */
	if(alpha == 1.) return 1.;

	/**< alpha in (0,1) */
	double U = unif_rand();
	double W;
	do { W = exp_rand(); } while (W == 0.);
	return pow(A_(M_PI*U,alpha)/pow(W,1.-alpha),1./alpha);
}

/**
 * Sample S ~ S(alpha, 1, 1, (1/cos(alpha*pi/2))*I_{alpha == 1}; 1)
 * with Laplace-Stieltjes transform exp(-(1/cos(alpha*pi/2))*t^alpha),
 * see Nolan's book for the parametrization.
 *
 * @param alpha parameter in (0,1]
 * @return S
 * @author Marius Hofert, Martin Maechler
 */
double rstable(double alpha){
    if(alpha == 1.)
	return 1./0.; /**< Inf or NaN */
    else
	return pow(cos(M_PI_2*alpha), -1./alpha) * rstable0(alpha);
}

/**
 * Vectorize rstable. Generate a vector of random variates
 * S ~ S(alpha, 1, 1, (1/cos(alpha*pi/2))*I_{alpha == 1}; 1)
 * with Laplace-Stieltjes transform exp(-(1/cos(alpha*pi/2))*t^alpha),
 * see Nolan's book for the parametrization.
 *
 * @param S vector of stable random variates (result)
 * @param n length of the vector S
 * @param alpha parameter in (0,1]
 * @return none
 * @author Marius Hofert, Martin Maechler
 */
void rstable_vec(double S[], const R_xlen_t n, const double alpha){
    if(n >= 1) {
	double cf = pow(cos(M_PI_2*alpha), -1./alpha);
	GetRNGstate();
	for(R_xlen_t i=0; i < n; i++)
	    S[i] = cf * rstable0(alpha);
	PutRNGstate();
    }
    return;
}

/**
 * Generate a vector of random variates
 * S ~ S(alpha, 1, 1, (1/cos(alpha*pi/2))*I_{alpha == 1}; 1)
 * with Laplace-Stieltjes transform exp(-(1/cos(alpha*pi/2))*t^alpha),
 * see Nolan's book for the parametrization.
 *
 * @param n sample size
 * @param alpha parameter in (0,1]
 * @return vector of random variates S
 * @author Martin Maechler
 */
SEXP rstable_c(SEXP n, SEXP alpha)
{
    int nn = asInteger(n);
    SEXP res = PROTECT(allocVector(REALSXP, nn));

    rstable_vec(REAL(res), nn, asReal(alpha));
    UNPROTECT(1);
    return(res);
}

/**
 * Sample St ~ \tilde{S}(alpha, 1, (cos(alpha*pi/2)*V_0)^{1/alpha},
 *			 V_0*I_{alpha = 1}, h*I_{alpha != 1}; 1)
 * with Laplace-Stieltjes transform exp(-V_0((h+t)^alpha-h^alpha)),
 * see Nolan's book for the parametrization, via the "fast rejection algorithm",
 * see Hofert (2012).
 *
 * @param St vector of random variates (result)
 * @param V0 vector of random variates V0
 * @param h parameter in [0,infinity)
 * @param alpha parameter in (0,1]
 * @param n length of St
 * @return none
 * @author Marius Hofert, Martin Maechler
 */
void retstable_MH(double *St, const double V0[], double h, double alpha, R_xlen_t n){
    /**
     * alpha == 1 => St corresponds to a point mass at V0 with Laplace-Stieltjes
     * transform exp(-V0*t)
    */
    if(alpha == 1.) {
	for(R_xlen_t i = 0; i < n; i++)
	    St[i] = V0[i];
	return;
    }

    GetRNGstate();

    for(R_xlen_t i = 0; i < n; i++) { /**< for each of the n required variates */

	/**< find m := optimal number of summands, using asymptotic formula */
	int m;
	m = imax2(1, (int)round(V0[i] * pow(h,alpha)));

	/**< apply standard rejection m times and sum the variates St_k */
	double c = pow(V0[i]/m, 1./alpha);
	St[i] = 0.; /**< will be the result after the summation */
	for(int k = 0; k < m; k++){
	    double St_k, U;
	    /**
	     * apply standard rejection for sampling
	     * \tilde{S}(alpha, 1, (cos(alpha*pi/2)*V_0/m)^{1/alpha},
	     *           (V_0/m)*I_{alpha = 1}, h*I_{alpha != 1}; 1) with
	     * Laplace-Stieltjes transform exp(-(V_0/m)*((h+t)^alpha-h^alpha))
	    */
	    do {
		/**
		 * sample St_k ~ S(alpha, 1, (cos(alpha*pi/2)*V_0/m)^{1/alpha},
		 *		   (V_0/m)*I_{\alpha = 1}; 1) with
		 * Laplace-Stieltjes transform exp(-(V_0/m)*t^alpha)
		*/
		St_k = c * rstable0(alpha);
		U = unif_rand();
	    } while (U > exp(- h * St_k));
	    /**
	     * on acceptance, St_k ~ \tilde{S}(alpha, 1, (cos(alpha*pi/2)
	     *				       *V_0/m)^{1/alpha},
	     *				       (V_0/m)*I_{alpha = 1},
	     *				       h*I_{alpha != 1}; 1) with
	     * Laplace-Stieltjes transform exp(-(V_0/m)*((h+t)^alpha-h^alpha))
	    */
	    St[i] += St_k;
	}

    } /**< end for */

    PutRNGstate();
    return;
}

/**
 * Fast and accurate evaluation of sinc(x) := sin(x)/x, including the limit x=0
 *
 * @param x any (double precision) number
 * @return sinc(x)
 * @author Martin Maechler (2010-04-28)
*/
double sinc_MM(double x) {
    double ax = fabs(x);
    if(ax < 0.006) {
	if(x == 0.) return 1;
	double x2 = x*x;
	if(ax < 2e-4)
	     return 1. - x2/6.;
	else return 1. - x2/6.*(1 - x2/20.);
    }
    /**< else */
    return sin(x)/x;
}

/** to be called from R */
SEXP sinc_c(SEXP x_) {
    R_xlen_t n = XLENGTH(PROTECT(x_ = isReal(x_) ?
				 Rf_duplicate(x_) : coerceVector(x_, REALSXP)));
    SEXP r_ = allocVector(REALSXP, n); // the result
    double *x = REAL(x_), *r = REAL(r_);

    for(R_xlen_t i=0; i < n; i++)
	r[i] = sinc_MM(x[i]);

    UNPROTECT(1);
    return r_;
}

/**
 * Evaluation of Zolotarev's function, see Devroye (2009), to the power 1-alpha.
 * The 3-arg. version allows more precision for  alpha ~=~ 1
 *
 * @param x argument
 * @param alpha parameter in (0,1]
 * @return sin(alpha*x)^alpha * sin((1-alpha)*x)^(1-alpha) / sin(x)
 * @author Martin Maechler (2010-04-28)
 */
#define _A_3(_x, _alpha_, _I_alpha)				\
    pow(_I_alpha* sinc_MM(_I_alpha*_x), _I_alpha) *		\
    pow(_alpha_ * sinc_MM(_alpha_ *_x), _alpha_ ) / sinc_MM(_x)

// called in retstable_LD()
double A_(double x, double alpha) {
  double Ialpha = 1.-alpha;
  return _A_3(x, alpha, Ialpha);
}

/** to be called from R---see experiments in ../tests/retstable-ex.R */
SEXP A__c(SEXP x_, SEXP alpha, SEXP I_alpha) {
    R_xlen_t n = XLENGTH(PROTECT(x_ = isReal(x_) ?
				 Rf_duplicate(x_) : coerceVector(x_, REALSXP)));
    double alp = asReal(alpha), I_alp = asReal(I_alpha);
    if(fabs(alp + I_alp - 1.) > 1e-12)
	error("'I_alpha' must be == 1 - alpha more accurately");
    SEXP r_ = allocVector(REALSXP, n); /**< the result */
    double *x = REAL(x_), *r = REAL(r_);

    for(R_xlen_t i=0; i < n; i++)
	r[i] = _A_3(x[i], alp, I_alp);

    UNPROTECT(1);
    return r_;
}

/**
 * Evaluation of B(x)/B(0), see Devroye (2009).
 *
 * @param x argument
 * @param alpha parameter in (0,1]
 * @return sinc(x) / (sinc(alpha*x)^alpha * sinc((1-alpha)*x)^(1-alpha))
 * @author Martin Maechler (2010-04-28)
 */
double BdB0(double x,double alpha) {
    double Ialpha = 1.-alpha;
    double den = pow(sinc_MM(alpha*x),alpha) * pow(sinc_MM(Ialpha*x),Ialpha);
    return sinc_MM(x) / den;
}

/**
 * Sample St ~ \tilde{S}(alpha, 1, (cos(alpha*pi/2)*V_0)^{1/alpha},
 *			 V_0*I_{alpha = 1}, h*I_{alpha != 1}; 1)
 * with Laplace-Stieltjes transform exp(-V_0((h+t)^alpha-h^alpha)),
 * see Nolan's book for the parametrization, via double rejection,
 * see Devroye (2009).
 *
 * @param St vector of random variates (result)
 * @param V0 vector of random variates V0
 * @param h parameter in [0,infinity)
 * @param alpha parameter in (0,1]
 * @param n length of St
 * @return none
 * @author Marius Hofert, Martin Maechler
 */
void retstable_LD(double *St, const double V0[], double h, double alpha, R_xlen_t n)
{
   /**
     * alpha == 1 => St corresponds to a point mass at V0 with Laplace-Stieltjes
     * transform exp(-V0*t)
   */
   if(alpha == 1.) {
       for(R_xlen_t i = 0; i < n; i++)
	   St[i] = V0[i];
       return;
   }

   // compute variables not depending on V0
   const double c1 = sqrt(M_PI_2);
   const double c2 = 2.+c1;
   double Ialpha = 1. - alpha, // "FIXME" : allow Ialpha to be argument
       b = Ialpha/alpha,
       h_a = pow(h, alpha);

   for(R_xlen_t i = 0; i < n; i++) { // for each of the n required variates

	/**< set lambda for our parameterization */
       double lambda_alpha = h_a*V0[i]; // < Marius Hofert: work directly with lambda^alpha (numerically more stable for small alpha)

	/**
	 * Apply the algorithm of Devroye (2009) to draw from
	 * \tilde{S}(alpha, 1, (cos(alpha*pi/2))^{1/alpha}, I_{alpha = 1},
	 * 	     lambda*I_{alpha != 1};1) with Laplace-Stieltjes transform
	 * exp(-((lambda+t)^alpha-lambda^alpha))
	*/
	double gamma = lambda_alpha*alpha*Ialpha;
	double sgamma = sqrt(gamma);
	double c3 = c2* sgamma;
	double xi = (1. + M_SQRT2 * c3)/M_PI; /**< according to John Lau */
	double psi = c3*exp(-gamma*M_PI*M_PI/8.)/M_SQRT_PI;
	double w1 = c1*xi/sgamma;
	double w2 = 2.*M_SQRT_PI * psi;
	double w3 = xi*M_PI;
	double X, c, E;
	do {
	    double U, z, Z;
	    do {
		double V = unif_rand();
		if(gamma >= 1) {
		    if(V < w1/(w1+w2)) U = fabs(norm_rand())/sgamma;
		    else{
			double W_ = unif_rand();
			U = M_PI*(1.-W_*W_);
		    }
		}
		else{
		    double W_ = unif_rand();
		    if(V < w3/(w2+w3)) U = M_PI*W_;
		    else U = M_PI*(1.-W_*W_);
		}
		double W = unif_rand();
		double zeta = sqrt(BdB0(U,alpha));
		z = 1/(1-pow(1+alpha*zeta/sgamma,-1/alpha)); /**< Marius Hofert: numerically more stable for small alpha */
		/**< compute rho */
		double rho = M_PI*exp(-lambda_alpha*(1. - 1./(zeta*zeta))) /
		    ((1.+c1)*sgamma/zeta + z);
		double d = 0.;
		if(U >= 0 && gamma >= 1) d += xi*exp(-gamma*U*U/2.);
		if(U > 0 && U < M_PI) d += psi/sqrt(M_PI-U);
		if(U >= 0 && U <= M_PI && gamma < 1) d += xi;
		rho *= d;
		Z = W*rho;
	    } while( !(U < M_PI && Z <= 1.)); /* check rejection condition */

	    double
		a = pow(_A_3(U,alpha,Ialpha), 1./Ialpha),
		m = pow(b/a,alpha)*lambda_alpha,
		delta = sqrt(m*alpha/a),
		a1 = delta*c1,
		a3 = z/a,
		s = a1+delta+a3;
	    double V_ = unif_rand(), N_ = 0., E_ = 0.; // -Wall
	    if(V_ < a1/s) {
		N_ = norm_rand();
		X = m-delta*fabs(N_);
	    } else {
		if(V_ < (a1+delta)/s) /**< according to John Lau */
		    X = m+delta*unif_rand();
		else {
		    E_ = exp_rand();
		    X = m+delta+E_*a3;
		}
	    }
	    E = -log(Z);
	    /**< check rejection condition */
		c = a*(X-m)+exp((1/alpha)*log(lambda_alpha)-b*log(m))*(pow(m/X,b)-1); /**< Marius Hofert: numerically more stable for small alpha */
	    if(X < m) c -= N_*N_/2.;
	    else if(X > m+delta) c -= E_;

	} while (!(X >= 0 && c <= E));
	/**
	 * Transform variates from the distribution corresponding to the
	 * Laplace-Stieltjes transform exp(-((lambda+t)^alpha-lambda^alpha))
	 * to those of the distribution corresponding to the Laplace-Stieltjes
	 * transform exp(-V_0((h+t)^alpha-h^alpha)).
	*/
	St[i] = exp(1/alpha*log(V0[i])-b*log(X)); /**< Marius Hofert: numerically more stable for small alpha */

    } /**< end for */
    return;
}

/**
 * Sample St ~ \tilde{S}(alpha, 1, (cos(alpha*pi/2)*V_0)^{1/alpha},
 *			 V_0*I_{alpha = 1}, h*I_{alpha != 1}; 1)
 * with Laplace-Stieltjes transform exp(-V_0((h+t)^alpha-h^alpha)),
 * see Nolan's book for the parametrization, via the algorithms of
 * Devroye (2009) and Hofert (2012).
 *
 * @param V0_ R vector of type numeric(n) containing the random variates V0
 * @param h R parameter of type numeric(1)
 * @param alpha R parameter in (0,1] of type numeric(1)
 * @param method R "string" of type character(1)
 * @return R vector of type numeric(n) containing the random variates
 * @author Martin Maechler
 */
SEXP retstable_c(SEXP V0_, SEXP h, SEXP alpha, SEXP method)
{
    R_xlen_t n = LENGTH(PROTECT(V0_ = isReal(V0_) ?
				Rf_duplicate(V0_) : coerceVector(V0_, REALSXP)));
    const char* meth_ch = CHAR(STRING_ELT(method, 0));
    enum method_kind { MH, LD
    } meth_typ = ((!strcmp(meth_ch, "MH")) ? MH :
		  ((!strcmp(meth_ch, "LD")) ? LD :
		   -1));
    SEXP St; // the result
    PROTECT(St = allocVector(REALSXP, n));
    switch(meth_typ) {
    case MH:
	retstable_MH(REAL(St), REAL(V0_), asReal(h), asReal(alpha), n);
	break;
    case LD:
	retstable_LD(REAL(St), REAL(V0_), asReal(h), asReal(alpha), n);
	break;
    default:
	error(_("retstable_c(): invalid '%s'"), "method");
    }

    UNPROTECT(2);
    return(St);
}
