/*#################################################################################
##
##   R package Copula by Jun Yan and Ivan Kojadinovic Copyright (C) 2008, 2009
##
##   This file is part of the R package copula.
##
##   The R package copula is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package copula is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with the R package copula. If not, see <http://www.gnu.org/licenses/>.
##
#################################################################################*/


/***********************************************************************
  
 Testing wheter a copula belongs to the Extreme-Value class
 
***********************************************************************/

#include <R.h>
#include <Rmath.h>
#include "Anfun.h"

/***********************************************************************
  
 Empirical copula
 
***********************************************************************/

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
  return (res + o)/(n + 1.0);
}

/***********************************************************************
  
 Derivative of the empirical copula
 
***********************************************************************/

double derec(double *U, int n, int p, double *u, double *v, double o)
{
  return (ec(U, n, p, u, o) - ec(U, n, p, v, o)) *  sqrt(n) / 2.0; 
}

/***********************************************************************
  
 Extreme-Value Test
 
***********************************************************************/

void evtest(double *U, int *n, int *p, double *g, int *m, 
	    int *N, double *tg,  int *nt, double *o, double *s0)
{
  double *influ = Calloc((*n) * (*m) * (*nt), double);
  double *random = Calloc(*n, double);
 
  double *u = Calloc(*p, double);
  double *v = Calloc(*p, double);
  double *ut = Calloc(*p, double);
  double *vt = Calloc(*p, double);
  double *der = Calloc(*p, double);
  double *dert = Calloc(*p, double);
 
  double t, ecterm, process, mean, ind, indt, d, dt, invsqrtn = 1.0/sqrt(*n);

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
	  
	  ecterm = t *  R_pow(ec(U, *n, *p, u, *o), t - 1); 
	  
	  /* derivatives */
	  for (k = 0; k < *p; k++)
	    {
	      u[k] += invsqrtn;
	      v[k] -= invsqrtn;
	      der[k] = derec(U, *n, *p, u, v, *o);
	      u[k] -= invsqrtn;
	      v[k] += invsqrtn;
	      
	      ut[k] += invsqrtn;
	      vt[k] -= invsqrtn;
	      dert[k] = derec(U, *n, *p, ut, vt, *o);
	      ut[k] -= invsqrtn;
	      vt[k] += invsqrtn;
	      
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
		  dt += dert[k] * (U[i + k * (*n)] <= ut[k]);
		}
	      influ[i + j * (*n) + c * (*n) * (*m)] 
		= (indt - dt - ecterm * (ind - d)) * invsqrtn;
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
  
 Extreme-Value test statistic
 
***********************************************************************/
    
void evtest_stat(double *U, int *n, int *p, double *g, int *m,  
		 double *tg,  int *nt, double *o, double *stat)
{
  int j, k, c;
  double diff, t;
  double *u = Calloc(*p, double);
  double *ut = Calloc(*p, double);

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
	  diff = ec(U, *n, *p, ut, *o) - R_pow(ec(U, *n, *p, u, *o), t);
	  stat[c] += diff * diff;
	}
      stat[c] = stat[c] * (*n) / (*m);
    }

  Free(u);
  Free(ut);
}


/***********************************************************************
  
 Extreme-Value Test based on Anfun
 
***********************************************************************/

double ecop(double *U, double *V, int n, double u, double v, double o)
{
  int i;
  double res = 0.0;
  
  for (i = 0; i < n; i++)
      res += (U[i] <= u) * (V[i] <= v);
  return (res + o)/(n + 1.0);
}

double der1ec(double *U, double *V, int n, double u, double v, double o)
{
  double invsqrtn = 1.0 / sqrt(n);
  return (ecop(U, V, n, u + invsqrtn, v, o) - ecop(U, V, n, u - invsqrtn, v, o)) 
    / (2.0 * invsqrtn);
}

double der2ec(double *U, double *V, int n, double u, double v, double o)
{
  double invsqrtn = 1.0 / sqrt(n);
  return (ecop(U, V, n, u ,v + invsqrtn, o) - ecop(U, V, n, u ,v - invsqrtn, o)) 
     / (2.0 * invsqrtn);
}

void evtestA(double *U, double *V, int *n, double *u, double *v, 
	     int *m, double *o, int *CFG, int *N, double *s0)
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

      d1 = der1ec(U, V, *n, u[j], v[j], *o);
      d2 = der2ec(U, V, *n, u[j], v[j], *o);

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

    
void evtestA_stat(double *U, double *V, int *n, double *u, double *v, int *m,  
		  double *o, int *CFG, double *stat)
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
      
      diff = ecop(U, V, *n, u[j], v[j], *o) - exp(loguv * Aj);
      s += diff * diff;
    }
  
  *stat = s * (*n) / *m;

  Free(S);
  Free(T);
}
