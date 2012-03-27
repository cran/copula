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


#include <R.h>
#include <Rmath.h>

#include "Anfun.h"

/*****************************************************************************

  Rank-based version of Pickands' estimator (inverse)

*****************************************************************************/

double inv_A_Pickands(int n, double *S, double *T, double t)
{
  int i;
  double At = 0.0, Eit, Sit, Tit;

  if (t > 0.0 && t < 1.0)
    {
      for (i=0;i<n;i++)
	{
	  Sit = S[i] / (1.0 - t);
	  Tit = T[i] / t;
	  Eit = (Sit < Tit ) ? Sit : Tit ;
	  At += Eit;
	}
    }
  else if (t <= 0.0) /* t set to 0 */
    {
      for (i=0;i<n;i++)
	At += S[i];
    }
  else /* t >= 1.0 */ /* t set to 1 */
    {
      for (i=0;i<n;i++)
	At += T[i];
    }
  return At/(double)n;
}

/*****************************************************************************

  Rank-based version of Pickands' estimator interfaced in R
  n: sample size
  t: vector argument
  A: result
  m: length of t and A

*****************************************************************************/

void A_Pickands(int *n, double *S, double *T, double *t, int *m,
		int *corrected, double *A)
{
  int i;
  double invA0, invA1;

  if (*corrected)
    {
      invA0 = inv_A_Pickands(*n, S, T, 0.0);
      invA1 = inv_A_Pickands(*n, S, T, 1.0);
      for (i=0;i<*m;i++)
	A[i] = 1.0/(inv_A_Pickands(*n, S, T, t[i])
		    - (1.0 - t[i]) * (invA0  - 1.0)
		    - t[i] * (invA1 - 1.0));
    }
  else
    for (i=0;i<*m;i++)
      A[i] = 1.0/inv_A_Pickands(*n, S, T, t[i]);

}

/*****************************************************************************

  Rank-based version of the CFG estimator (log)

*****************************************************************************/

/* Euler's constant */
#define EULER 0.5772156649015328606065120900824

double log_A_CFG(int n, double *S, double *T, double t)
{
  int i;
  double At = 0.0;

  if (0. < t && t < 1.)
    {
      for (i=0;i<n;i++) {
	  double
	      Sit = S[i] / (1.0 - t),
	      Tit = T[i] / t,
	      Eit = (Sit < Tit ) ? Sit : Tit ;
	  At += log(Eit);
      }
    }
  else if (t <= 0.0) /* t set to 0 */
    {
      for (i=0;i<n;i++)
	At += log(S[i]);
    }
  else /* t == 1.0 */ /* t set to 1 */
    {
      for (i=0;i<n;i++)
	At += log(T[i]);
    }

  return -EULER - At/(double)n;
}

/*****************************************************************************

  Rank-based version of the CFG estimator interfaced in R
  n: sample size
  t: vector argument
  A: result
  m: length of t and A

*****************************************************************************/

void A_CFG(int *n, double *S, double *T, double *t, int *m,
	   int *corrected, double *A)
{
  int i;
  if (*corrected) {
      double
	  logA0 = log_A_CFG(*n, S, T, 0.0),
	  logA1 = log_A_CFG(*n, S, T, 1.0);

      for (i=0;i<*m;i++)
	  A[i] = exp(log_A_CFG(*n, S, T, t[i])
		     - (1.0 - t[i]) * logA0 - t[i] * logA1);
  }
  else
    for (i=0;i<*m;i++)
      A[i] = exp(log_A_CFG(*n, S, T, t[i]));
}
