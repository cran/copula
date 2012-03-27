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


/***********************************************************************

 Mult

***********************************************************************/

#include <R.h>
#include <Rmath.h>

#define MIN(x,y) ((x)<(y)?(x):(y))
#define MAX(x,y) ((x)>(y)?(x):(y))

/***********************************************************************

 Empirical copula

***********************************************************************/

// FIXME: provide R interface to this !
double empc(double *U, int n, int p, double *u)
{
  int i,j;
  double ind, res = 0.0;

  for (i = 0; i < n; i++)
    {
      ind = 1.0;
      for (j = 0; j < p; j++)
	ind *= (U[i + n * j] <= u[j]);
      res += ind;
    }
  return res/n;
}

/***********************************************************************

 Derivative of the empirical copula

***********************************************************************/

static
double derempc(double *U, int n, int p, double *u, double *v, double denom)
{
  return (empc(U, n, p, u) - empc(U, n, p, v)) / denom;
}

/***********************************************************************



***********************************************************************/

void mult(double *U, int *n, int *p, double *g, int *m,
	  int *N, double *s0, double *s1, int *deriv, int *multi)
{
  double *influ = Calloc((*n) * (*m), double);
  double *indp = Calloc((*n) * (*m), double);
  double *random = Calloc(*n, double);

  double *u = Calloc(*p, double);
  double *v = Calloc(*p, double);
  double *der = Calloc(*p, double);

  double process, mean, d, denom, invsqrtn = 1.0/sqrt(*n),
    tmpu, tmpv, ecu, ecu2;

  int i, j, k, l;

  /* for each point of the grid */
  for (j = 0; j < *m; j++)
    {

      /* temporary arrays */
      for (k = 0; k < *p; k++)
	{
	  u[k] = g[j + k * (*m)];
	  v[k] = u[k];
	}


      /* derivatives */
      for (k = 0; k < *p; k++)
	{
	  /* bins = 2 invsqrtn */
	  if (*deriv == 1)
	    {
	      denom = 2 * invsqrtn;

	      /* u and v */
	      if (u[k] < invsqrtn)
		{
		  tmpu = u[k]; tmpv = v[k];
		  u[k] = 2 * invsqrtn; v[k] = 0.0;
		  der[k] = derempc(U, *n, *p, u, v, denom);
		  u[k] = tmpu; v[k] = tmpv;
		}
	      else if (u[k] > 1 - invsqrtn)
		{
		  tmpu = u[k]; tmpv = v[k];
		  u[k] = 1; v[k] = 1 - 2 * invsqrtn;
		  der[k] = derempc(U, *n, *p, u, v, denom);
		  u[k] = tmpu; v[k] = tmpv;
		}
	      else
		{
		  u[k] += invsqrtn; v[k] -= invsqrtn;
		  der[k] = derempc(U, *n, *p, u, v, denom);
		  u[k] -= invsqrtn; v[k] += invsqrtn;
		}
	    }
	  else if (*deriv == 2)/* bins decreasing near 0 and 1 */
	    {
	      u[k] += invsqrtn; v[k] -= invsqrtn;
	      denom = MIN(u[k], 1.0) - MAX(v[k], 0.0);
	      der[k] = derempc(U, *n, *p, u, v, denom);
	      u[k] -= invsqrtn; v[k] += invsqrtn;
	    }
	  else /* no correction at 0 and 1 */
	    {
	      u[k] += invsqrtn; v[k] -= invsqrtn;
	      denom = 2 * invsqrtn;
	      der[k] = derempc(U, *n, *p, u, v, denom);
	      u[k] -= invsqrtn; v[k] += invsqrtn;
	    }
	}

      /* for each pseudo-obs */
      for (i = 0; i < *n; i++)
	{
	  indp[i + j * (*n)] = 1.0;
	  d = 0.0;
	  for (k = 0; k < *p; k++)
	    {
	      indp[i + j * (*n)] *= (U[i + k * (*n)] <= u[k]);
	      d += der[k] *  (U[i + k * (*n)] <= u[k]);
	    }
	  influ[i + j * (*n)] = (indp[i + j * (*n)] - d) * invsqrtn;
	}
    }

  GetRNGstate();

  /* generate N approximate realizations */
  for (l = 0; l < *N; l++)
    {
      /* generate n variates */
      mean = 0.0;
      for (i=0;i<*n;i++)
	{
	  if (*multi == 1)
	    random[i] = norm_rand();
	  else
	    random[i] = (unif_rand() < 0.5) ? -1.0 : 1.0 ;
	  mean += random[i];
	}
      mean /= *n;

      /* realization number l */
      for (j = 0; j < *m; j++)
	{
	  /* multiplier */
	  process = 0.0;
	  for (i = 0; i < *n; i++)
	    process += (random[i] - mean) * influ[i + j * (*n)];
	  s0[j + l * (*m)] = process;

	  /* new multiplier */
	  for (k = 0; k < *p; k++)
	      u[k] = g[j + k * (*m)];
	  ecu = empc(U, *n, *p, u);

	  process = 0.0;
	  for (i = 0; i < *n; i++)
	      process += (random[i] - mean) * indp[i + j * (*n)];
	      process *= invsqrtn;

	  for (k = 0; k < *p; k++)
	    {
	      tmpu = 0.0;
	      for (i = 0; i < *n; i++)
		tmpu += (random[i] - mean) * (U[i + k * (*n)] <= u[k]);
	      tmpu /= *n;
	      u[k] += -tmpu + ceil( (*n) * u[k] ) / *n - u[k];
	    }
	  ecu2 = empc(U, *n, *p, u);

	  s1[j + l * (*m)] = process + (ecu2 - ecu) / invsqrtn;

	}
    }

  PutRNGstate();

  Free(influ);
  Free(indp);
  Free(random);
  Free(u);
  Free(v);
  Free(der);
}
