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
 * @file   indepTest_utils.c
 * @author Ivan Kojadinovic
 * @date   December 2007
 *
 * @brief  Utility functions for multivariate and "vectorial"
 *         tests of independence and serial independence based
 *         on the empirical copula process
 *
 */

#include <R.h>
#include "indepTests.h"


/// Temporary array K
void K_array(int n, int p, const double J[], double *K)
{
  int m=0, n2 = n * n;
  for (int j=0; j < p; j++)
    for (int i=0; i < n; i++)
      {
	K[m] = 0.0;
	for (int l=0; l < n; l++)
	  K[m] +=  J[n2 * j + n * l + i];
	K[m++] /= (double) n;
      }
}

/// Temporary array L
void L_array(int n, int p, const double K[], double *L)
{
  for (int j=0; j < p; j++)
    {
      L[j] = 0.0;
      for (int i=0; i < n; i++)
	L[j] += K[n * j + i];
      L[j] /= (double)n;
    }
}

/**
 * Computation of the global Cramer-von Mises statistic.
 *
 * @param n sample size
 * @param p dimension
 * @param J temporary array
 * @param K temporary array
 * @param L temporary array
 * @return the global Cramer-von Mises statistic
 * @author Ivan Kojadinovic
 */
double I_n(int n, int p, double *J, double *K, double *L)
{
  int i, j, l, n2 = n * n;
  double In, sum, prod, part1, part2, part3;

  /* first term */
  sum = 0.0;
  for (i = 0; i < n; i++)
    for (l = 0; l < n; l++)
      {
	prod = 1.0;
	for (j = 0; j < p; j++)
	  prod *= J[n2 * j + n * l + i];
	sum += prod;
      }
  part1 = sum / (double) (n);

  /* second term */
  sum = 0.0;
  for (i = 0; i < n;i++) {
    prod = 1.0;
    for (j = 0; j < p; j++)
      prod *= K[j * n + i]; /* K(i, j) */
    sum += prod;
  }
  part2 = 2.0 * sum;


  /* third term */
  prod = 1.0;
  for (j = 0; j < p; j++) prod *= L[j];
  part3 = prod * (double) n;

  In = part1 - part2 + part3;

  return In;
}

/**
 * Computes the Cramer-von Mises statistics for
 * subset A of {1,...,p}, |A|>1 (A \ni 1?)
 *
 * @param n sample size
 * @param p dimension
 * @param J temporary array
 * @param K temporary array
 * @param L temporary array
 * @param A the (integer coded) subset A
 * @return the Cramer-von Mises statistic for subset A
 * @author Ivan Kojadinovic
 */
double M_A_n(int n, int p, double *J, double *K, double *L, int A)
{
  size_t n2 = n * n;// size_t should not overflow
  double MAn = 0.0;
  for (int i=0; i < n; i++)
    for (int l=0; l < n; l++)
      {
	double prod = 1.0;
	for (int j=0; j < p; j++)
	  if (1<<j & A)
	    prod *= J[n2 * j + n * l + i] - K[n * j + i] - K[n * j + l] + L[j];
	MAn += prod;
      }

  return MAn/(double)n;
}


/**
 * Text progress bar for k=1 to N loops
 *
 * @param k current iteration
 * @param N maximum number of iteration
 * @param w width of the bar in characters
 * @author Ivan Kojadinovic
 */
inline void progressBar(int k, int N, int w) {
  if (k < N-1) {
    // length of one of "block"
    int b = (N <= w) ? 1 : N/w;

    // update only r times
    if ( k % b != 0 ) return;

    /// percentage done
    double percent = k / (double)N;
    int c = (int)(percent * w);

    /// display the bar
    Rprintf("  |");
    for (int i=0; i<c; i++) Rprintf("=");
    for (int i=c; i<w; i++) Rprintf(" ");
    Rprintf("| %3d%%\r", (int)(percent*100) );
    // \r = carriage return  withOUT newline
  }
  else { /// display the bar for last k (= N-1)
    Rprintf("  |");
    for (int i=0; i<w; i++)
      Rprintf("=");
    Rprintf("| 100%%\n");
  }
}
