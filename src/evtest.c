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
 * @file   evtest.c
 * @author Ivan Kojadinovic
 * @date   2011
 *
 * @brief Tests of extreme-value dependence
 *
 */

#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <R_ext/Utils.h>
#include "copula.h"
#include "empcop.h"

/**
 * Test of extreme-value dependence based on the empirical copula
 *
 * @param U pseudo-observation
 * @param n sample size
 * @param p dimension
 * @param g grid
 * @param m grid size
 * @param N number of multiplier replications
 * @param tg array of powers (related to max-stability)
 * @param nt number of powers
 * @param s0 N replicates of the test statistic under the null
 * @param der2n if > 0, bins of size 2/sqrt(n), otherwise smaller near 0/1
 * @param o offset
 * @param stat value of the test statistic
 * @author Ivan Kojadinovic
 */
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
    invsqrtn = 1./sqrt(*n), tmpu, tmpv, diff, ecut, ecu;

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

	  ecterm = R_pow(multCn(U, *n, *p, ut, 1, 0, 0.), (1 - t)/t) / t;

	  /* derivatives */
	  for (k = 0; k < *p; k++)
	    {
	      /* bins = 2 invsqrtn */
	      if (*der2n)
		{
		  denom = 2. * invsqrtn;

		  /* u and v */
		  if (u[k] < invsqrtn)
		    {
		      tmpu = u[k]; tmpv = v[k];
		      u[k] = 2. * invsqrtn; v[k] = 0.;
		      der[k] = der_multCn(U, *n, *p, u, v, denom);
		      u[k] = tmpu; v[k] = tmpv;
		    }
		  else if (u[k] > 1. - invsqrtn)
		    {
		      tmpu = u[k]; tmpv = v[k];
		      u[k] = 1.; v[k] = 1. - 2. * invsqrtn;
		      der[k] = der_multCn(U, *n, *p, u, v, denom);
		      u[k] = tmpu; v[k] = tmpv;
		    }
		  else
		    {
		      u[k] += invsqrtn; v[k] -= invsqrtn;
		      der[k] = der_multCn(U, *n, *p, u, v, denom);
		      u[k] -= invsqrtn; v[k] += invsqrtn;
		    }

		  /* ut and vt */
		  if (ut[k] < invsqrtn)
		    {
		      tmpu = ut[k]; tmpv = vt[k];
		      ut[k] = 2. * invsqrtn; vt[k] = 0.;
		      dert[k] = der_multCn(U, *n, *p, ut, vt, denom);
		      ut[k] = tmpu; vt[k] = tmpv;
		    }
		  else if (ut[k] > 1. - invsqrtn)
		    {
		      tmpu = ut[k]; tmpv = vt[k];
		      ut[k] = 1.; vt[k] = 1. - 2. * invsqrtn;
		      dert[k] = der_multCn(U, *n, *p, ut, vt, denom);
		      ut[k] = tmpu; vt[k] = tmpv;
		    }
		  else
		    {
		      ut[k] += invsqrtn; vt[k] -= invsqrtn;
		      dert[k] = der_multCn(U, *n, *p, ut, vt, denom);
		      ut[k] -= invsqrtn; vt[k] += invsqrtn;
		    }
		}
	      else /* bins decreasing near 0 and 1 */
		{
		  u[k] += invsqrtn; v[k] -= invsqrtn;
		  denom = MIN(u[k], 1.) - MAX(v[k], 0.);
		  der[k] = der_multCn(U, *n, *p, u, v, denom);
		  u[k] -= invsqrtn; v[k] += invsqrtn;

		  ut[k] += invsqrtn; vt[k] -= invsqrtn;
		  denom = MIN(ut[k], 1.) - MAX(vt[k], 0.);
		  dert[k] = der_multCn(U, *n, *p, ut, vt, denom);
		  ut[k] -= invsqrtn; vt[k] += invsqrtn;
		}
	    }

	  /* for each pseudo-obs */
	  for (i = 0; i < *n; i++)
	    {
	      ind = 1.;
	      indt = 1.;
	      d = 0.;
	      dt = 0.;
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
      mean = 0.;
      for (i=0;i<*n;i++)
	{
	  random[i] = norm_rand(); /*(unif_rand() < 0.5) ? -1. : 1. ;*/
	  mean += random[i];
	}
      mean /= *n;

      /* realization number l */
      for (c = 0; c < *nt; c++)
	{
	  k = c + l * (*nt);
	  s0[k] = 0.;
	  for (j = 0; j < *m; j++)
	    {
	      process = 0.;
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
      stat[c] = 0.;
      t = tg[c];
      for (j = 0; j < *m; j++)
	{
	  for (k = 0; k < *p; k++)
	    {
	      u[k] = g[j + k * (*m)];
	      ut[k] =  R_pow(u[k], t);
	    }
	  ecut = multCn(U, *n, *p, ut, 1, 0, *o);
	  ecu = multCn(U, *n, *p, u, 1, 0, *o);
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

/**
 * Test of extreme-value dependence based on An
 * Derivatives based on Cn - see JMVA paper
 *
 * @param U pseudo-observations
 * @param V pseudo-observations
 * @param n sample size
 * @param u grid
 * @param v grid
 * @param m grid size
 * @param CFG if > 0, then An=CFG, otherwise Pickands
 * @param N number of multiplier replications
 * @param s0 N replications of the test statistic
 * @author Ivan Kojadinovic
 */
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
    mean, invsqrtn = 1. / sqrt(*n), minTkSi,  minSkTi,
    lb = 1. / (*n + 1.), ub = *n / (*n + 1.),
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
      cA0 = biv_logACFG(*n, S, T, 0.);
      cA1 = biv_logACFG(*n, S, T, 1.);
    }
  else
    {
      cA0 = biv_invAP(*n, S, T, 0.);
      cA1 = biv_invAP(*n, S, T, 1.);
    }

  /* for each point of the grid */
  for (j = 0; j < *m; j++)
    {
      loguv = log(u[j] * v[j]);

      pu = loguv / log(u[j]);
      pv = loguv / log(v[j]);

      d1 = der1bivCn(U, V, *n, u[j], v[j]);
      d2 = der2bivCn(U, V, *n, u[j], v[j]);

      t = 1. / pv;

      if (*CFG)
	{
	  A = exp(biv_logACFG(*n, S, T, t)
		  - (1. - t) * cA0 - t * cA1);
	  loguvA = loguv * A;
	  Aterm = exp(loguvA) * loguvA;
	}
      else
	{
	  A = 1. / (biv_invAP(*n, S, T, t)
		     - (1. - t) * (cA0  - 1.)
		     - t * (cA1 - 1.));
	  loguvA = loguv * A;
	  Aterm = exp(loguvA) * loguvA * A;
	}

      for (i = 0; i < *n; i++) /* for each pseudo-obs */
	{
	  sum = 0.;
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
	  sum *= invsqrtn / 2.;

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
      mean = 0.;
      for (i=0;i<*n;i++)
	{
	  random[i] = norm_rand(); /*(unif_rand() < 0.5) ? -1. : 1. ;*/
	  mean += random[i];
	}
      mean /= *n;

      /* realization number l */
      s0[l] = 0.;
      for (j = 0; j < *m; j++)
	{
	  process = 0.;
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

/* FIXME ?? -- as the integrand uses indicators and hence is discontinuous, the numerical integrator gives tons
   of 'ier=1' accuracy problem errors.   Rdqags() source (ier=1) mentions that one really should split the interval
   into piecewise continuous sub-intervals one self !!
 */

// The integrand [called from vec_intgr() from evtestA_derA()] :
double intgr(double x, double termUt, double termVt, double powUt, double powVt,
	     double U, double V, double t, double n)
{
  double
      x1t = R_pow(x,1-t),  indUt = (U <= x1t) - (int)(x1t * (n+1)) / n,
      xt  = R_pow(x, t),   indVt = (V <= xt ) - (int)(xt  * (n+1)) / n,
      res = 0.;
  if(indUt || indVt) {
      double xlogx = x * log(x);
      if (indUt) res += termUt * R_pow(x,powUt) * indUt / xlogx;
      if (indVt) res += termVt * R_pow(x,powVt) * indVt / xlogx;
  }
  return res;
}


// Vectorized intgr(), called from numerical integrator Rdqags(), from evtestA_derA
void vec_intgr(double *x, int n, void *ex)
{
  double *arg = (double *) ex;
  for (int i = 0; i < n; i++)
    x[i] = intgr(x[i],
		 arg[0], arg[1], arg[2], arg[3],
		 arg[4], arg[5], arg[6], arg[7]);
  return;
}

/**
 * Test of extreme-value dependence based on An
 * Derivatives based on An - see JMVA paper
 *
 * @param U pseudo-observations
 * @param V pseudo-observations
 * @param n sample size
 * @param u grid
 * @param v grid
 * @param m grid size
 * @param CFG if != 0, then An=CFG, otherwise Pickands
 * @param trace_lev: if > 0, print "progress" output to console
 * @param N number of multiplier replications
 * @param s0 ["Output"]: N replications of the test statistic
 * @author Ivan Kojadinovic
 */
void evtestA_derA(double *U, double *V, int *n,
		  double *u, double *v, int *m, int *CFG_tr__etc,
		  int *N, double *s0)
{
    int CFG       = CFG_tr__etc[0],
	trace_lev = CFG_tr__etc[1],
	report_err= CFG_tr__etc[2];
  double  *influ = Calloc((*n) * (*m), double), *influ2;
  influ2 = CFG ? Calloc((*n) * (*m), double) : NULL;
  double *random = Calloc(*n, double);

  /* for numerical integration [Rdqags() <==> R's integrate()]: begin */
  double result, abserr;
  int last, neval, ier;
  double reltol=0.0001;
  double abstol=0.0001;
  int limit=100;         int   *iwork = Calloc(limit, int);
  int lenw = 4 * limit;  double *work = Calloc(lenw,  double);
  double ex[8],
      lower = 0., // = 1. / (*n + 1.)
      upper = 1.; // = *n / (*n + 1.)
  int nErr = 0, errCnts[6] = {0,0, 0,0, 0,0};
  /* for numerical integration: end */

  double // temporary arrays  S[], T[] :
      *S = Calloc(*n, double),
      *T = Calloc(*n, double);
  for (int i = 0; i < *n; i++)
    {
      S[i] = - log(U[i]);
      T[i] = - log(V[i]);
    }

  /* correction terms */
  double cA0, cA1;
  if (CFG)
    {
      cA0 = biv_logACFG(*n, S, T, 0.);
      cA1 = biv_logACFG(*n, S, T, 1.);
    }
  else
    {
      cA0 = biv_invAP(*n, S, T, 0.);
      cA1 = biv_invAP(*n, S, T, 1.);
    }

  /* for each point of the grid */
  for (int j = 0; j < *m; j++) {
      double
	  uv = u[j] * v[j],
	  loguv = log(uv),
	  pu = loguv / log(u[j]),
	  pv = loguv / log(v[j]);

      if(trace_lev)
	  REprintf(trace_lev >= 2 ? "\n>-> j=%d: pv=%g ...\n" : "j=%d, pv=%11g.\n",
		   j, pv);
      double
	  t = 1. / pv,
	  invsqrtn = 1. / sqrt(*n),
	  tp = t + invsqrtn,
	  tm = t - invsqrtn;

      double A, loguvA, Aterm, Ap, Am, dAt, termUt, termVt, powUt, powVt;
      if (CFG)
	{
	  A = exp(biv_logACFG(*n, S, T, t)
		  - (1. - t) * cA0 - t * cA1);
	  loguvA = loguv * A;
	  Aterm = exp(loguvA) * loguvA;

	  Ap = exp(biv_logACFG(*n, S, T, tp)
		   - (1. - tp) * cA0 - tp * cA1);
	  Am = exp(biv_logACFG(*n, S, T, tm)
		   - (1. - tm) * cA0 - tm * cA1);
	  dAt = (Ap - Am) / (2. * invsqrtn);

	  termUt = A - t * dAt;
	  termVt = A + (1. - t) * dAt;
	  powUt = A + t - 1.;
	  powVt = A - t;
	}
      else
	{
	  A = 1. / (biv_invAP(*n, S, T, t)
		     - (1. - t) * (cA0  - 1.)
		     - t * (cA1 - 1.));
	  loguvA = loguv * A;
	  Aterm = exp(loguvA) * loguvA * A;

	  Ap = 1. / (biv_invAP(*n, S, T, tp)
		       - (1. - tp) * (cA0  - 1.)
		       - tp * (cA1 - 1.));
	  Am = 1. / (biv_invAP(*n, S, T, tm)
		       - (1. - tm) * (cA0  - 1.)
		       - tm * (cA1 - 1.));
	  dAt = (Ap - Am) / (2. * invsqrtn);

	  termUt = (A - t * dAt) / (A + t - 1.);
	  termVt = (A + (1. - t) * dAt) / (A - t);
	  powUt = (A + t - 1.) / (1. - t);
	  powVt = (A - t) / t;
	}

     R_CheckUserInterrupt();
      /*d1 = (A  -    t   * dAt) * R_pow(uv, A - 1 + t);
	d2 = (A + (1 - t) * dAt) * R_pow(uv, A - t);*/
      double
	  d1 = der1bivCn(U, V, *n, u[j], v[j]),
	  d2 = der2bivCn(U, V, *n, u[j], v[j]);

      for (int i = 0; i < *n; i++) { /* for each pseudo-obs */
	  if (CFG) {
	      if(i == 0) {
		  ex[0] = termUt; ex[1] = termVt;
		  ex[2] = powUt;  ex[3] = powVt;
		  ex[6] = t;      ex[7] = *n;
	      }
	      ex[4] = U[i];   ex[5] = V[i];

	      Rdqags(vec_intgr, (void *)ex, &lower, &upper,
		     &abstol,  &reltol,
		     &result,  &abserr,  &neval,  &ier,
		     &limit,  &lenw, &last,
		     iwork, work);

	      if(ier && report_err) {
		  // ier=1  happens a lot (99%) -- because integrand has jumps
		  /* if(ier > 1) */
		  /*     REprintf(" i=%d: Rdqags()=%lf, abserr=%lf, ier=%d *error*\n", i, */
		  /* 	       result, abserr, ier); */
		  nErr++;
		  errCnts[ier-1]++;
	      }
	      if(trace_lev >= 2)
		  REprintf(" i=%d: Rdqags()=%10g, abserr=%6g, neval=%d, last=%d\n", i,
			   result, abserr, neval, last);

	      influ[i + j * (*n)] = (U[i] <= u[j]) * (V[i] <= v[j])
		  - d1 * (U[i] <= u[j]) - d2 * (V[i] <= v[j])
		  + Aterm * log(MIN(pu * S[i], pv * T[i]));
	      influ2[i + j * (*n)] =  Aterm * result * invsqrtn;
	  }
	  else {
	      influ[i + j * (*n)] = (U[i] <= u[j]) * (V[i] <= v[j])
		  - d1 * (U[i] <= u[j]) - d2 * (V[i] <= v[j])
		  + Aterm * (MIN(pu * S[i], pv * T[i])
			     - termUt * (1. - R_pow(U[i],powUt))
			     - termVt * (1. - R_pow(V[i],powVt)));
	  }
	  influ [i + j * (*n)] *= invsqrtn;
      } // for i

    } // for j
  if(trace_lev)
      REprintf("\n---> completed influ[%d,%d] integration computations\n", *n, *m);
  if(nErr && report_err) {
      #define REPRINT6(IEX)				\
        for(int i=0; i < 6; i++) REprintf(" %3d", IEX);	\
        REprintf("\n")
      REprintf("Had %d integration errors (= %4.1f%%) with the following counts of each ier=1,..,6:\n",
	       nErr, 100. * nErr / ((double) *n * *m) );
      REPRINT6(   i+1    );
      REPRINT6(errCnts[i]);
  }
  GetRNGstate();

  /* generate N approximate realizations */
  if(trace_lev >= 2) REprintf("Random Realizations:\n");
  for (int l = 0; l < *N; l++) {
      /* generate n variates */
      double mean = 0.;
      for (int i=0; i < *n; i++) {
	  random[i] = norm_rand(); /*(unif_rand() < 0.5) ? -1. : 1. ;*/
	  mean += random[i];
      }
      mean /= *n;
      if(trace_lev >= 2) {
	  if(l % 50 == 1) REprintf(" l=%d\n", l); else REprintf(".");
      }
      /* realization number l */
      s0[l] = 0.;
      for (int j = 0; j < *m; j++) {
	  double process = 0.;
	  for (int i = 0; i < *n; i++)
	      if (CFG)
		  process += (random[i] - mean) * influ[i + j * (*n)]
		      + random[i] * influ2[i + j * (*n)];
	      else
		  process += (random[i] - mean) * influ[i + j * (*n)];

	  s0[l] += process * process;
      }
      s0[l] /= *m;
  }
  if(trace_lev)
      REprintf("\n---> Finished N=%d random realizations s0[]\n", *N);
  PutRNGstate();


  Free(influ);
  if(CFG) Free(influ2);
  Free(random);

  Free(S);
  Free(T);
}

/***********************************************************************

 Extreme-Value Test based on An
 Test statistic

***********************************************************************/
/**
 * Test statistic for the test of extreme-value
 * dependence based on An - see JMVA paper
 *
 * @param U pseudo-observations
 * @param V pseudo-observations
 * @param n sample size
 * @param u grid
 * @param v grid
 * @param m grid size
 * @param CFG if > 0, then An=CFG, otherwise Pickands
 * @param stat value of the test statistic
 * @param offset offset for the test statistic
 * @author Ivan Kojadinovic
 */
void evtestA_stat(double *U, double *V, int *n, double *u, double *v, int *m,
		  int *CFG, double *stat, double *offset)
{
  int i, j;
  double s = 0., diff, cA0, cA1, Aj, t, loguv;
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
      cA0 = biv_logACFG(*n, S, T, 0.);
      cA1 = biv_logACFG(*n, S, T, 1.);
    }
  else
    {
      cA0 = biv_invAP(*n, S, T, 0.);
      cA1 = biv_invAP(*n, S, T, 1.);
    }

  for (j = 0; j < *m; j++)
    {
      loguv = log(u[j] *v[j]);
      t = log(v[j]) / loguv;

      if (*CFG)
	Aj = exp(biv_logACFG(*n, S, T, t)
		 - (1. - t) * cA0 - t * cA1);
      else
	Aj = 1. / (biv_invAP(*n, S, T, t)
		    - (1. - t) * (cA0  - 1.)
		    - t * (cA1 - 1.));
      if (*offset < 0.)
	diff = bivCn(U, V, *n, u[j], v[j]) - exp(loguv * Aj);
      else
	diff = bivCn(U, V, *n, u[j], v[j]) * (*n) / (*n+1)
	  + (*offset)/(*n+1) - exp(loguv * Aj);

      s += diff * diff;
    }

  *stat = s * (*n) / *m;

  Free(S);
  Free(T);
}
