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
 * @param alpha in (0,1)x
 * @return random variate following the logarithmic (series) distribution Log(alpha)
 * @author Marius Hofert
 */
static
int rlogseries_LS ( double alpha ) {
  double t = - alpha / log(1 - alpha);
  double u = runif(0.0, 1.0), p = t;
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
 * @author Marius Hofert
 */
static
int rlogseries_LK ( double alpha ) {
  double h = log(1 - alpha);
  double u2 = runif(0.0, 1.0), x = 1.0, u1, q;
  if ( u2 > alpha ) return (int) x;
  u1 = runif(0.0, 1.0);
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
 * @author Marius Hofert
 */
void rlogseries_R (int *n, double *alpha, int *val) {
  int i;
  double thres = 0.95;// FIXME : thres should be *argument*

  GetRNGstate();
  for (i = 0; i < *n; i++)
      val[i] = ( alpha[i] < thres ) ?
	  rlogseries_LS ( alpha[i] ) :
	  rlogseries_LK ( alpha[i] );
  PutRNGstate();
}
