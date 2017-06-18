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


/* Kemp, A. W., 1981, Applied Statistics 30(3) p. 249--253 */

#include <math.h>
#include <R.h>
#include <Rmath.h>
#include "copula.h"

/**
 * Algorithm LS of Kemp (1981)
 *
 * @param alpha in (0,1)
 * @return random variate following the logarithmic (series) distribution Log(alpha)
 */
static
int rlogseries_LS ( double alpha ) {
  double t = - alpha / log(1 - alpha);
  double u = unif_rand(), p = t;
  int x = 1;
  while (u > p) {
    u = u - p;
    x = x + 1;
    p = p * alpha * (x - 1) / x;
  }
  return (x);
}

/**
 * Algorithm LK of Kemp (1981)
 *
 * @param alpha in (0,1)
 * @return random variate following the logarithmic (series) distribution Log(alpha)
 */
static
int rlogseries_LK ( double alpha ) {
  double h = log(1 - alpha);
  double u2 = unif_rand(), x = 1.0, u1, q;
  if ( u2 > alpha ) return (int) x;
  u1 = unif_rand();
  q = 1 - exp(u1 * h);
  if ( u2 < q * q ) return (int) (1 + log(u2) / log(q));
  if ( u2 > q ) return 1;
  return 2;
}

/**
 * R wrapper
 *
 * @param n sample size
 * @param alpha parameters in (0,1)
 * @return val return value(s)
 */
void rlogseries_R (int *n, double *alpha, int *val) {
  int i;
  double thres = 0.95;// FIXME : thres should be *argument*
  // alpha < 0.95  <==>  h = log(1 - alpha) > log(0.05) = -2.996

  GetRNGstate();
  for (i = 0; i < *n; i++)
      val[i] = ( alpha[i] < thres ) ?
	  rlogseries_LS ( alpha[i] ) :
	  rlogseries_LK ( alpha[i] );
  PutRNGstate();
}

//_____________ FIXME ____________ rather add these to ./rLog.c ______________
//                                                       ======

// From R's sources in  src/nmath/dpq.h :
// log(1 - exp(x))  in more stable form than log1p(- R_D_qIv(x)) :
#define R_Log1_Exp(x)   ((x) > -M_LN2 ? log(-expm1(x)) : log1p(-exp(x)))


/**
 * Algorithm LK of Kemp (1981) -- using h = log(1 - alpha) as parameter
 *
 * @param h in (-oo, 0)
 * @return random variate following the logarithmic (series) distribution Log(1 - e^h)
 * @author Martin Maechler
 */
static
double rlogseries_LK_ln1p (double h) { // h == log(1 - alpha) < 0 <==> alpha = 1 - e^h = -expm(h)
  double alpha = -expm1(h);
  double u2 = unif_rand();
  if (u2 > alpha)
      return 1.;
  double u1 = unif_rand(),
      h_ = u1 * h,
      q = - expm1(h_);
  if ( u2 < q * q ) {
      double log_q = R_Log1_Exp(h_);
      return (log_q == 0.) ? R_PosInf : 1. + floor(log(u2) / log_q);
  } else // u2 >= q*q
      return ( u2 > q ) ? 1. : 2.;
}


/**
 * R wrapper
 *
 * @param n sample size
 * @param h in (-oo, 0),  h = log(1 - alpha) -- a scalar --
 * @return val return value(s) -- double because 'int' may overflow
 * @author Martin Maechler
 */
void rlogseries_R_ln1p (int *n, double *h, double *val) {
    // h == log(1 - alpha) < 0 <==> alpha = 1 - e^h = -expm1(h)
    // alpha < 0.95 <==>  1 - e^h  < 0.95  <==>  e^h > 0.05  <==> h > log(0.05) = -2.995732
  int i;

  GetRNGstate();
  if(*h > -3.) { // <==> alpha < 0.9502
      double alpha = -expm1(*h);
      for (i = 0; i < *n; i++)
	  val[i] = rlogseries_LS(alpha);
  }
  else // h <= -3
      for (i = 0; i < *n; i++)
	  val[i] = rlogseries_LK_ln1p(*h);
  PutRNGstate();
}

