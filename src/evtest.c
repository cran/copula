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

 Testing whether a copula belongs to the Extreme-Value class

***********************************************************************/

#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>

#include "copula.h"

/***********************************************************************

 Empirical copula

***********************************************************************/
// FIXME: we have the empirical copula in about every *.c file ... -- use *one* !
double ec(double *U, int n, int p, double *u, double o)
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
  return res/(n + o);
}

/***********************************************************************

 Derivative of the empirical copula

***********************************************************************/

double derec(double *U, int n, int p, double *u, double *v, double denom)
{
  return (ec(U, n, p, u, 0.0) - ec(U, n, p, v, 0.0)) / denom;
}

/***********************************************************************

 Extreme-Value Test

***********************************************************************/

void evtest(double *U, int *n, int *p, double *g, int *m,
	    int *N, double *tg,  int *nt, double *s0, int *der2n,
	    double *o, double *stat)
{
  double *influ = Calloc((*n) * (*m) * (*nt), double);
  double *random = Calloc(*n, double);

  double *u = Calloc(*p, double);
  double *v = Calloc(*p, double);
  double *ut = Calloc(*p, double);
  double *vt = Calloc(*p, double);
  double *der = Calloc(*p, double);
  double *dert = Calloc(*p, double);

  double t, ecterm, process, mean, ind, indt, d, dt, denom,
    invsqrtn = 1.0/sqrt(*n), tmpu, tmpv, diff, ecut, ecu;

  int i, j, k, l, c;

  for (c = 0; c < *nt; c++)
    {
      t = tg[c];

      /* for each point of the grid */
      for (j = 0; j < *m; j++)
	{

	  /* temporary arrays */
	  for (k = 0; k < *p; k++)
	    {
	      u[k] = g[j + k * (*m)];
	      v[k] = u[k];
	      ut[k] =  R_pow(u[k], t);
	      vt[k] =  ut[k];
	    }

	  ecterm = R_pow(ec(U, *n, *p, ut, 0.0), (1 - t)/t) / t;

	  /* derivatives */
	  for (k = 0; k < *p; k++)
	    {
	      /* bins = 2 invsqrtn */
	      if (*der2n)
		{
		  denom = 2.0 * invsqrtn;

		  /* u and v */
		  if (u[k] < invsqrtn)
		    {
		      tmpu = u[k]; tmpv = v[k];
		      u[k] = 2.0 * invsqrtn; v[k] = 0.0;
		      der[k] = derec(U, *n, *p, u, v, denom);
		      u[k] = tmpu; v[k] = tmpv;
		    }
		  else if (u[k] > 1.0 - invsqrtn)
		    {
		      tmpu = u[k]; tmpv = v[k];
		      u[k] = 1.0; v[k] = 1.0 - 2.0 * invsqrtn;
		      der[k] = derec(U, *n, *p, u, v, denom);
		      u[k] = tmpu; v[k] = tmpv;
		    }
		  else
		    {
		      u[k] += invsqrtn; v[k] -= invsqrtn;
		      der[k] = derec(U, *n, *p, u, v, denom);
		      u[k] -= invsqrtn; v[k] += invsqrtn;
		    }

		  /* ut and vt */
		  if (ut[k] < invsqrtn)
		    {
		      tmpu = ut[k]; tmpv = vt[k];
		      ut[k] = 2.0 * invsqrtn; vt[k] = 0.0;
		      dert[k] = derec(U, *n, *p, ut, vt, denom);
		      ut[k] = tmpu; vt[k] = tmpv;
		    }
		  else if (ut[k] > 1.0 - invsqrtn)
		    {
		      tmpu = ut[k]; tmpv = vt[k];
		      ut[k] = 1.0; vt[k] = 1.0 - 2.0 * invsqrtn;
		      dert[k] = derec(U, *n, *p, ut, vt, denom);
		      ut[k] = tmpu; vt[k] = tmpv;
		    }
		  else
		    {
		      ut[k] += invsqrtn; vt[k] -= invsqrtn;
		      dert[k] = derec(U, *n, *p, ut, vt, denom);
		      ut[k] -= invsqrtn; vt[k] += invsqrtn;
		    }
		}
	      else /* bins decreasing near 0 and 1 */
		{
		  u[k] += invsqrtn; v[k] -= invsqrtn;
		  denom = MIN(u[k], 1.0) - MAX(v[k], 0.0);
		  der[k] = derec(U, *n, *p, u, v, denom);
		  u[k] -= invsqrtn; v[k] += invsqrtn;

		  ut[k] += invsqrtn; vt[k] -= invsqrtn;
		  denom = MIN(ut[k], 1.0) - MAX(vt[k], 0.0);
		  dert[k] = derec(U, *n, *p, ut, vt, denom);
		  ut[k] -= invsqrtn; vt[k] += invsqrtn;
		}
	    }

	  /* for each pseudo-obs */
	  for (i = 0; i < *n; i++)
	    {
	      ind = 1.0;
	      indt = 1.0;
	      d = 0.0;
	      dt = 0.0;
	      for (k = 0; k < *p; k++)
		{
		  ind *= (U[i + k * (*n)] <= u[k]);
		  indt *= (U[i + k * (*n)] <= ut[k]);
		  d += der[k] * (U[i + k * (*n)] <= u[k]);
		  dt += dert[k]  * (U[i + k * (*n)] <= ut[k]);
		}
	      influ[i + j * (*n) + c * (*n) * (*m)]
		= (ecterm * (indt - dt) - (ind - d)) * invsqrtn;
	    }
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
	  random[i] = norm_rand(); /*(unif_rand() < 0.5) ? -1.0 : 1.0 ;*/
	  mean += random[i];
	}
      mean /= *n;

      /* realization number l */
      for (c = 0; c < *nt; c++)
	{
	  k = c + l * (*nt);
	  s0[k] = 0.0;
	  for (j = 0; j < *m; j++)
	    {
	      process = 0.0;
	      for (i = 0; i < *n; i++)
		process += (random[i] - mean) * influ[i + j * (*n) + c * (*n) * (*m)];

	      s0[k] += process * process;
	    }
	  s0[k] /= *m;
	}
    }

  PutRNGstate();

  /* test statistic */
  for (c = 0; c < *nt; c++)
    {
      stat[c] = 0.0;
      t = tg[c];
      for (j = 0; j < *m; j++)
	{
	  for (k = 0; k < *p; k++)
	    {
	      u[k] = g[j + k * (*m)];
	      ut[k] =  R_pow(u[k], t);
	    }
	  ecut = ec(U, *n, *p, ut, *o);
	  ecu = ec(U, *n, *p, u, *o);
	  //diff = (R_pow(ecut, 1/t) * (*n + 1 - 1/t) -  ecu * (*n)) / (*n + 0.3);
	  diff = R_pow(ecut, 1/t) - ecu;
	  stat[c] += diff * diff;
	}
      stat[c] = stat[c] * (*n) / (*m);
    }

  Free(influ);
  Free(random);
  Free(u);
  Free(v);
  Free(ut);
  Free(vt);
  Free(der);
  Free(dert);
}

/***********************************************************************

 Extreme-Value Test based on An
 Derivatives based on Cn

***********************************************************************/

double ecop(double *U, double *V, int n, double u, double v)
{
  int i;
  double res = 0.0;

  for (i = 0; i < n; i++)
      res += (U[i] <= u) * (V[i] <= v);
  return res/n;
}

double der1ec(double *U, double *V, int n, double u, double v)
{
  double invsqrtn = 1.0 / sqrt(n);
  return (ecop(U, V, n, u + invsqrtn, v) - ecop(U, V, n, u - invsqrtn, v))
    / (2.0 * invsqrtn);
}

double der2ec(double *U, double *V, int n, double u, double v)
{
  double invsqrtn = 1.0 / sqrt(n);
  return (ecop(U, V, n, u ,v + invsqrtn) - ecop(U, V, n, u ,v - invsqrtn))
     / (2.0 * invsqrtn);
}

void evtestA(double *U, double *V, int *n, double *u, double *v,
	     int *m, int *CFG, int *N, double *s0)
{
  double *influ = Calloc((*n) * (*m), double);

  double *S = Calloc(*n, double);
  double *Sp = Calloc(*n, double);
  double *Sm = Calloc(*n, double);
  double *T = Calloc(*n, double);
  double *Tp = Calloc(*n, double);
  double *Tm = Calloc(*n, double);

  double *random = Calloc(*n, double);

  double pu, pv, sum, d1, d2,  process,
    mean, invsqrtn = 1.0 / sqrt(*n), minTkSi,  minSkTi,
    lb = 1.0 / (*n + 1.0), ub = *n / (*n + 1.0),
    A, cA0, cA1, t, loguv, loguvA, Aterm;

  int i, j, k, l;

  /* temporary arrays */
  for (i = 0; i < *n; i++)
    {
      S[i] = - log(U[i]);
      T[i] = - log(V[i]);
      Sp[i] = - log(MIN(U[i] + invsqrtn, ub));
      Tp[i] = - log(MIN(V[i] + invsqrtn, ub));
      Sm[i] = - log(MAX(U[i] - invsqrtn, lb));
      Tm[i] = - log(MAX(V[i] - invsqrtn, lb));
    }

  /* correction terms */
  if (*CFG)
    {
      cA0 = log_A_CFG(*n, S, T, 0.0);
      cA1 = log_A_CFG(*n, S, T, 1.0);
    }
  else
    {
      cA0 = inv_A_Pickands(*n, S, T, 0.0);
      cA1 = inv_A_Pickands(*n, S, T, 1.0);
    }

  /* for each point of the grid */
  for (j = 0; j < *m; j++)
    {
      loguv = log(u[j] * v[j]);

      pu = loguv / log(u[j]);
      pv = loguv / log(v[j]);

      d1 = der1ec(U, V, *n, u[j], v[j]);
      d2 = der2ec(U, V, *n, u[j], v[j]);

      t = 1.0 / pv;

      if (*CFG)
	{
	  A = exp(log_A_CFG(*n, S, T, t)
		  - (1.0 - t) * cA0 - t * cA1);
	  loguvA = loguv * A;
	  Aterm = exp(loguvA) * loguvA;
	}
      else
	{
	  A = 1.0 / (inv_A_Pickands(*n, S, T, t)
		     - (1.0 - t) * (cA0  - 1.0)
		     - t * (cA1 - 1.0));
	  loguvA = loguv * A;
	  Aterm = exp(loguvA) * loguvA * A;
	}

      for (i = 0; i < *n; i++) /* for each pseudo-obs */
	{
	  sum = 0.0;
	  for (k = 0; k < *n; k++)
	    {
	      minTkSi = MIN(pv * T[k], pu * S[i]);
	      minSkTi = MIN(pu * S[k], pv * T[i]);
	      if (*CFG)
		sum += -log(MIN(pu * Sm[k], minTkSi))
		  + log(MIN(pu * Sp[k], minTkSi))
		  - log(MIN(pv * Tm[k], minSkTi))
		  + log(MIN(pv * Tp[k], minSkTi));
	      else
		sum += MIN(pu * Sm[k], minTkSi) - MIN(pu * Sp[k], minTkSi)
		  + MIN(pv * Tm[k], minSkTi) - MIN(pv * Tp[k], minSkTi);
	    }
	  sum *= invsqrtn / 2.0;

	  if (*CFG)
	    influ[i + j * (*n)] = (U[i] <= u[j]) * (V[i] <= v[j])
	      - d1 * (U[i] <= u[j]) - d2 * (V[i] <= v[j])
	      - Aterm * (-log(MIN(pu * S[i], pv * T[i])) - sum);
	  else
	    influ[i + j * (*n)] = (U[i] <= u[j]) * (V[i] <= v[j])
	      - d1 * (U[i] <= u[j]) - d2 * (V[i] <= v[j])
	      + Aterm * (MIN(pu * S[i], pv * T[i]) - sum);

	  influ[i + j * (*n)] *= invsqrtn;
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
	  random[i] = norm_rand(); /*(unif_rand() < 0.5) ? -1.0 : 1.0 ;*/
	  mean += random[i];
	}
      mean /= *n;

      /* realization number l */
      s0[l] = 0.0;
      for (j = 0; j < *m; j++)
	{
	  process = 0.0;
	  for (i = 0; i < *n; i++)
	    process += (random[i] - mean) * influ[i + j * (*n)];

	  s0[l] += process * process;
	}
      s0[l] /= *m;
    }

  PutRNGstate();


  Free(influ);
  Free(random);

  Free(S);
  Free(T);
  Free(Sp);
  Free(Tp);
  Free(Sm);
  Free(Tm);
}

/***********************************************************************

 Extreme-Value Test based on An
 Derivatives based on An

***********************************************************************/

double intgr(double x, double termUt, double termVt, double powUt, double powVt,
	     double U, double V, double t, double n)
{
  double x1t = R_pow(x,1-t), xt = R_pow(x,t);
  double indUt = (U <= x1t) - (int)(x1t * (n+1)) / n;
  double indVt = (V <= xt) - (int)(xt * (n+1)) / n;
  double xlogx = x * log(x), res = 0.0;
  if (indUt) res +=  termUt * R_pow(x,powUt) * indUt / xlogx;
  if (indVt) res += termVt * R_pow(x,powVt) * indVt / xlogx;
  return res;
}

void vec_intgr(double *x, int n, void *ex)
{
  int i;
  double *arg = ex;
  for (i = 0; i < n; i++)
    x[i] = intgr(x[i], arg[0], arg[1], arg[2], arg[3],
		 arg[4], arg[5], arg[6], arg[7]);
  return;
}

void evtestA_derA(double *U, double *V, int *n, double *u, double *v,
	     int *m, int *CFG, int *N, double *s0)
{
  double *influ = Calloc((*n) * (*m), double);
  double *influ2 = Calloc((*n) * (*m), double);

  double *S = Calloc(*n, double);
  double *T = Calloc(*n, double);

  double *random = Calloc(*n, double);

  /* for numerical integration: begin */
  double result, abserr;
  int last, neval, ier;
  int limit=100;
  double reltol=0.0001;
  double abstol=0.0001;
  int lenw = 4 * limit;
  int *iwork = Calloc(limit, int);
  double *work = Calloc(lenw, double);
  double ex[8], lower = 0.0, upper = 1.0;
  /* lower = 1.0 / (*n + 1), upper = *n / (*n + 1.0); */
  /* for numerical integration: end */

  double pu, pv, d1, d2,  process,
    mean, invsqrtn = 1.0 / sqrt(*n),
    A, cA0, cA1, t, tp, tm, uv, loguv, loguvA,
    Aterm, dAt, Ap, Am, termUt, termVt, powUt, powVt;

  int i, j, l;

  /* temporary arrays */
  for (i = 0; i < *n; i++)
    {
      S[i] = - log(U[i]);
      T[i] = - log(V[i]);
    }

  /* correction terms */
  if (*CFG)
    {
      cA0 = log_A_CFG(*n, S, T, 0.0);
      cA1 = log_A_CFG(*n, S, T, 1.0);
    }
  else
    {
      cA0 = inv_A_Pickands(*n, S, T, 0.0);
      cA1 = inv_A_Pickands(*n, S, T, 1.0);
    }

  /* for each point of the grid */
  for (j = 0; j < *m; j++)
    {
      uv = u[j] * v[j];
      loguv = log(uv);

      pu = loguv / log(u[j]);
      pv = loguv / log(v[j]);

      t = 1.0 / pv;
      tp = t + invsqrtn;
      tm = t - invsqrtn;

      if (*CFG)
	{
	  A = exp(log_A_CFG(*n, S, T, t)
		  - (1.0 - t) * cA0 - t * cA1);
	  loguvA = loguv * A;
	  Aterm = exp(loguvA) * loguvA;

	  Ap = exp(log_A_CFG(*n, S, T, tp)
		   - (1.0 - tp) * cA0 - tp * cA1);
	  Am = exp(log_A_CFG(*n, S, T, tm)
		   - (1.0 - tm) * cA0 - tm * cA1);
	  dAt = (Ap - Am) / (2.0 * invsqrtn);

	  termUt = A - t * dAt;
	  termVt = A + (1.0 - t) * dAt;
	  powUt = A + t - 1.0;
	  powVt = A - t;
	}
      else
	{
	  A = 1.0 / (inv_A_Pickands(*n, S, T, t)
		     - (1.0 - t) * (cA0  - 1.0)
		     - t * (cA1 - 1.0));
	  loguvA = loguv * A;
	  Aterm = exp(loguvA) * loguvA * A;

	  Ap = 1.0 / (inv_A_Pickands(*n, S, T, tp)
		       - (1.0 - tp) * (cA0  - 1.0)
		       - tp * (cA1 - 1.0));
	  Am = 1.0 / (inv_A_Pickands(*n, S, T, tm)
		       - (1.0 - tm) * (cA0  - 1.0)
		       - tm * (cA1 - 1.0));
	  dAt = (Ap - Am) / (2.0 * invsqrtn);

	  termUt = (A - t * dAt) / (A + t - 1.0);
	  termVt = (A + (1.0 - t) * dAt) / (A - t);
	  powUt = (A + t - 1.0) / (1.0 - t);
	  powVt = (A - t) / t;
	}

      /*d1 = (A - t * dAt) * R_pow(uv, A - 1.0 + t);
	d2 = (A + (1.0 - t) * dAt) * R_pow(uv, A - t);*/

      /* tmp */
      d1 = der1ec(U, V, *n, u[j], v[j]);
      d2 = der2ec(U, V, *n, u[j], v[j]);

      for (i = 0; i < *n; i++) /* for each pseudo-obs */
	{
	  if (*CFG)
	    {
	      ex[0] = termUt; ex[1] = termVt;
	      ex[2] = powUt; ex[3] = powVt;
	      ex[4] = U[i]; ex[5] = V[i];
	      ex[6] = t; ex[7] = *n;

	      Rdqags(vec_intgr, (void *)ex, &lower, &upper,
		     &abstol,  &reltol,
		     &result,  &abserr,  &neval,  &ier,
		     &limit,  &lenw, &last,
		     iwork, work);

	      /* if (ier)
		 Rprintf("%lf %lf %d\n", result, abserr, ier); */

	      influ[i + j * (*n)] = (U[i] <= u[j]) * (V[i] <= v[j])
		- d1 * (U[i] <= u[j]) - d2 * (V[i] <= v[j])
		+ Aterm * log(MIN(pu * S[i], pv * T[i]));
	      influ2[i + j * (*n)] =  Aterm * result;
	    }
	  else
	    influ[i + j * (*n)] = (U[i] <= u[j]) * (V[i] <= v[j])
	      - d1 * (U[i] <= u[j]) - d2 * (V[i] <= v[j])
	      + Aterm * (MIN(pu * S[i], pv * T[i])
			 - termUt * (1.0 - R_pow(U[i],powUt))
			 - termVt * (1.0 - R_pow(V[i],powVt)));

	  influ[i + j * (*n)] *= invsqrtn;
	  influ2[i + j * (*n)] *=  invsqrtn;
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
	  random[i] = norm_rand(); /*(unif_rand() < 0.5) ? -1.0 : 1.0 ;*/
	  mean += random[i];
	}
      mean /= *n;

      /* realization number l */
      s0[l] = 0.0;
      for (j = 0; j < *m; j++)
	{
	  process = 0.0;
	  for (i = 0; i < *n; i++)
	    if (*CFG)
	      process += (random[i] - mean) * influ[i + j * (*n)]
		+ random[i] * influ2[i + j * (*n)];
	    else
	      process += (random[i] - mean) * influ[i + j * (*n)];

	  s0[l] += process * process;
	}
      s0[l] /= *m;
    }

  PutRNGstate();


  Free(influ);
  Free(influ2);
  Free(random);

  Free(S);
  Free(T);
}

/***********************************************************************

 Extreme-Value Test based on An
 Test statistic

***********************************************************************/

void evtestA_stat(double *U, double *V, int *n, double *u, double *v, int *m,
		  int *CFG, double *stat, double *offset)
{
  int i, j;
  double s = 0.0, diff, cA0, cA1, Aj, t, loguv;
  double *S = Calloc(*n, double);
  double *T = Calloc(*n, double);

  for (i = 0; i < *n; i++)
    {
      S[i] = - log(U[i]);
      T[i] = - log(V[i]);
    }

  /* correction terms */
  if (*CFG)
    {
      cA0 = log_A_CFG(*n, S, T, 0.0);
      cA1 = log_A_CFG(*n, S, T, 1.0);
    }
  else
    {
      cA0 = inv_A_Pickands(*n, S, T, 0.0);
      cA1 = inv_A_Pickands(*n, S, T, 1.0);
    }

  for (j = 0; j < *m; j++)
    {
      loguv = log(u[j] *v[j]);
      t = log(v[j]) / loguv;

      if (*CFG)
	Aj = exp(log_A_CFG(*n, S, T, t)
		 - (1.0 - t) * cA0 - t * cA1);
      else
	Aj = 1.0 / (inv_A_Pickands(*n, S, T, t)
		    - (1.0 - t) * (cA0  - 1.0)
		    - t * (cA1 - 1.0));
      if (*offset < 0.0)
	diff = ecop(U, V, *n, u[j], v[j]) - exp(loguv * Aj);
      else
	diff = ecop(U, V, *n, u[j], v[j]) * (*n) / (*n+1)
	  + (*offset)/(*n+1) - exp(loguv * Aj);


      s += diff * diff;
    }

  *stat = s * (*n) / *m;

  Free(S);
  Free(T);
}


/***********************************************************************

 Extreme-Value Test based on An CFG - An Pickands
 Derivatives based on Cn

***********************************************************************/

void evTestAA_C(double *U, double *V, int *n, double *t, int *m,
		int *N, double *s0)
{
  double *influ = Calloc((*n) * (*m), double);
  double *random = Calloc(*n, double);

  double *S = Calloc(*n, double);
  double *T = Calloc(*n, double);
  double *Sp = Calloc(*n, double);
  double *Tp = Calloc(*n, double);
  double *Sm = Calloc(*n, double);
  double *Tm = Calloc(*n, double);

  double sumcfg, sump, Acfg, Ap, process, mean,
    invsqrtn = 1.0/sqrt(*n), minTkSi, minSkTi,
    lb = 1.0 / (*n + 1.0), ub = *n / (*n + 1.0),
    cA0cfg, cA1cfg, cA0p, cA1p;

  int i, j, k, l;

  /* temporary arrays */
  for (i = 0; i < *n; i++)
    {
      S[i] = - log(U[i]);
      T[i] = - log(V[i]);
      Sp[i] = - log(MIN(U[i] + invsqrtn, ub));
      Tp[i] = - log(MIN(V[i] + invsqrtn, ub));
      Sm[i] = - log(MAX(U[i] - invsqrtn, lb));
      Tm[i] = - log(MAX(V[i] - invsqrtn, lb));
    }

  /* correction terms */
  cA0cfg = log_A_CFG(*n, S, T, 0.0);
  cA1cfg = log_A_CFG(*n, S, T, 1.0);
  cA0p = inv_A_Pickands(*n, S, T, 0.0);
  cA1p = inv_A_Pickands(*n, S, T, 1.0);

  /* for each point of the grid */
  for (j = 0; j < *m; j++)
    {

      Acfg = exp(log_A_CFG(*n, S, T, t[j])
		 - (1.0 - t[j]) * cA0cfg - t[j] * cA1cfg);
      Ap = 1.0 / (inv_A_Pickands(*n, S, T, t[j])
		  - (1.0 - t[j]) * (cA0p  - 1.0)
		  - t[j] * (cA1p - 1.0));

      for (i = 0; i < *n; i++) /* for each pseudo-obs */
	{
	  sumcfg = 0.0;
	  sump = 0.0;
	  for (k = 0; k < *n; k++)
	    {
	      minTkSi = MIN(T[k] / t[j], S[i] / (1.0 - t[j]));
	      minSkTi = MIN(S[k] / (1.0 - t[j]), T[i] / t[j]);

	      sumcfg += - log(MIN(Sm[k] / (1.0 - t[j]), minTkSi))
		+ log(MIN(Sp[k] / (1.0 - t[j]), minTkSi))
		- log(MIN(Tm[k] / t[j], minSkTi))
		+ log(MIN(Tp[k] / t[j], minSkTi));

	      sump += MIN(Sm[k] / (1.0 - t[j]), minTkSi)
		- MIN(Sp[k] / (1.0 - t[j]), minTkSi)
		+ MIN(Tm[k] / t[j], minSkTi)
		- MIN(Tp[k] / t[j], minSkTi);
	    }
	  sumcfg *= invsqrtn / 2.0;
	  sump *= invsqrtn / 2.0;

	  /* influence term */
	  influ[i + j * (*n)] =
	    Acfg * (-log(MIN(S[i] / (1.0 - t[j]), T[i] / t[j])) - sumcfg)
	    + Ap * Ap * (MIN(S[i] / (1.0 - t[j]), T[i] / t[j]) - sump);

	  influ[i + j * (*n)] *=  invsqrtn;
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
	  random[i] = norm_rand(); /*(unif_rand() < 0.5) ? -1.0 : 1.0 ;*/
	  mean += random[i];
	}
      mean /= (double)(*n);

      /* realization number l */
      s0[l] = 0.0;
      for (j = 0; j < *m; j++)
	{
	  process = 0.0;
	  for (i = 0; i < *n; i++)
	    process += (random[i] - mean) * influ[i + j * (*n)];

	  s0[l] += process * process;
	}
      s0[l] /= (double)(*m);
    }

  PutRNGstate();


  Free(influ);
  Free(random);
  Free(S);
  Free(T);
  Free(Sp);
  Free(Tp);
  Free(Sm);
  Free(Tm);
}

/***********************************************************************

 Extreme-Value Test based on An CFG - An Pickands
 Derivatives based on An

***********************************************************************/

double ntgrand(double x, double termUcfg, double termVcfg,
	      double powUcfg, double powVcfg, double U, double V,
	      double t, double n)
{
  double x1t = R_pow(x,1-t), xt = R_pow(x,t);
  double indUcfg = (U <= x1t) - (int)(x1t * (n+1)) / n;
  double indVcfg = (V <= xt) - (int)(xt * (n+1)) / n;
  double xlogx = x * log(x), res = 0.0;
  if (indUcfg) res += termUcfg * R_pow(x,powUcfg) * indUcfg / xlogx;
  if (indVcfg) res += termVcfg * R_pow(x,powVcfg) * indVcfg / xlogx;
  return res;
}

void vec_ntgrand(double *x, int n, void *ex)
{
  int i;
  double *arg = ex;
  for (i = 0; i < n; i++)
    x[i] = ntgrand(x[i], arg[0], arg[1], arg[2], arg[3],
		    arg[4], arg[5], arg[6], arg[7]);
  return;
}

void evTestAA_derA(double *U, double *V, int *n, double *t, int *m,
		   int *N, double *s0)
{
  double *influ = Calloc((*n) * (*m), double);
  double *influ2 = Calloc((*n) * (*m), double);
  double *random = Calloc(*n, double);

  double *S = Calloc(*n, double);
  double *T = Calloc(*n, double);

  double Acfg, Ap, Acfgm, Acfgp, Apm, App, dAcfg, dAp,
    tj, tjp, tjm, termUcfg, termVcfg, termUp, termVp,
    powUcfg, powVcfg, powUp, powVp, process, mean,
    invsqrtn = 1.0/sqrt(*n), cA0cfg, cA1cfg, cA0p, cA1p;

  int i, j, l;

  /* for numerical integration: begin */
  double result, abserr;
  int last, neval, ier;
  int limit=100;
  double reltol=0.0001;
  double abstol=0.0001;
  int lenw = 4 * limit;
  int *iwork = Calloc(limit, int);
  double *work = Calloc(lenw, double);
  double ex[8], lower = 0.0, upper = 1.0;
  /* lower = 1.0 / (*n + 1), upper = *n / (*n + 1.0); */
  /* for numerical integration: end */

  /* temporary arrays */
  for (i = 0; i < *n; i++)
    {
      S[i] = - log(U[i]);
      T[i] = - log(V[i]);
    }

  /* correction terms */
  cA0cfg = log_A_CFG(*n, S, T, 0.0);
  cA1cfg = log_A_CFG(*n, S, T, 1.0);
  cA0p = inv_A_Pickands(*n, S, T, 0.0);
  cA1p = inv_A_Pickands(*n, S, T, 1.0);

  /* for each point of the grid */
  for (j = 0; j < *m; j++)
    {
      tj = t[j];
      tjp = tj + invsqrtn;
      tjm = tj - invsqrtn;

      Acfg = exp(log_A_CFG(*n, S, T, tj)
		 - (1.0 - tj) * cA0cfg - tj * cA1cfg);
      Acfgp = exp(log_A_CFG(*n, S, T, tjp)
		  - (1.0 - tjp) * cA0cfg - tjp * cA1cfg);
      Acfgm = exp(log_A_CFG(*n, S, T, tjm)
		  - (1.0 - tjm) * cA0cfg - tjm * cA1cfg);
      dAcfg = (Acfgp - Acfgm) / (2.0 * invsqrtn);
      termUcfg = Acfg - tj * dAcfg;
      termVcfg = Acfg + (1.0 - tj) * dAcfg;
      powUcfg = Acfg + tj - 1.0;
      powVcfg = Acfg - tj;


      Ap = 1.0 / (inv_A_Pickands(*n, S, T, tj)
		  - (1.0 - tj) * (cA0p  - 1.0)
		  - t[j] * (cA1p - 1.0));
      App = 1.0 / (inv_A_Pickands(*n, S, T, tjp)
		   - (1.0 - tjp) * (cA0p  - 1.0)
		   - tjp * (cA1p - 1.0));
      Apm = 1.0 / (inv_A_Pickands(*n, S, T, tjm)
		   - (1.0 - tjm) * (cA0p  - 1.0)
		   - tjm * (cA1p - 1.0));
      dAp = (App - Apm) / (2.0 * invsqrtn);
      termUp = (Ap - tj * dAp) / (Ap + tj - 1.0);
      termVp = (Ap + (1.0 - tj) * dAp) / (Ap - tj);
      powUp = (Ap + tj - 1.0) / (1.0 - tj);
      powVp = (Ap - tj) / tj;

      for (i = 0; i < *n; i++) /* for each pseudo-obs */
	{
	  ex[0] = termUcfg; ex[1] = termVcfg;
	  ex[2] = powUcfg; ex[3] = powVcfg;
	  ex[4] = U[i]; ex[5] = V[i];
	  ex[6] = t[j]; ex[7] = *n;

	  Rdqags(vec_ntgrand, (void *)ex, &lower, &upper,
		 &abstol,  &reltol,
		 &result,  &abserr,  &neval,  &ier,
		 &limit,  &lenw, &last,
		 iwork, work);

	  /*if (ier)
	    Rprintf("%lf %lf %d\n", result, abserr, ier); */

	  /* influence matrices */
	  influ[i + j * (*n)] = - Acfg * log(MIN(S[i] / (1.0 - t[j]), T[i] / t[j]))
	    + Ap * Ap * (MIN(S[i]/ (1.0 - t[j]), T[i] / t[j])
			 - termUp * (1.0 - R_pow(U[i],powUp))
			 - termVp * (1.0 - R_pow(V[i],powVp)));

	  influ2[i + j * (*n)] =  - Acfg * result;

	  influ[i + j * (*n)] *=  invsqrtn;
	  influ2[i + j * (*n)] *=  invsqrtn;
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
	  random[i] = norm_rand(); /*(unif_rand() < 0.5) ? -1.0 : 1.0 ;*/
	  mean += random[i];
	}
      mean /= (double)(*n);

      /* realization number l */
      s0[l] = 0.0;
      for (j = 0; j < *m; j++)
	{
	  process = 0.0;
	  for (i = 0; i < *n; i++)
	    process += (random[i] - mean) * influ[i + j * (*n)]
	      + random[i] * influ2[i + j * (*n)];

	  s0[l] += process * process;
	}
      s0[l] /= (double)(*m);
    }

  PutRNGstate();


  Free(influ);
  Free(influ2);
  Free(random);
  Free(S);
  Free(T);
  Free(iwork);
  Free(work);
}

/***********************************************************************

 Extreme-Value Test based on An CFG - An Pickands
 Test statistic

***********************************************************************/

void evTestAA_stat(double *S, double *T, int *n, double *t, int *m,
		   double *stat)
{
  int j;
  double s = 0.0, diff, cA0cfg, cA1cfg, cA0p, cA1p;

  /* correction terms */
  cA0cfg = log_A_CFG(*n, S, T, 0.0);
  cA1cfg = log_A_CFG(*n, S, T, 1.0);
  cA0p = inv_A_Pickands(*n, S, T, 0.0);
  cA1p = inv_A_Pickands(*n, S, T, 1.0);

  for (j = 0; j < *m; j++)
    {
      diff = exp(log_A_CFG(*n, S, T, t[j])
		 - (1.0 - t[j]) * cA0cfg - t[j] * cA1cfg)
	- 1.0/(inv_A_Pickands(*n, S, T, t[j])
	       - (1.0 - t[j]) * (cA0p  - 1.0)
	       - t[j] * (cA1p - 1.0));
      s += diff * diff;
    }

  *stat = s * (double)(*n) / (double)(*m);
}
