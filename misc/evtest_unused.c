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
