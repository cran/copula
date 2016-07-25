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


#ifndef COPULA_GOF_H
#define COPULA_GOF_H

#include <R.h>

void cramer_vonMises(int *n, int *p, double *U, double *Ctheta,
		     double *stat);
void cramer_vonMises_grid(int *p, double *U, int *n, double *V, int *m,
			  double *Ctheta, double *stat);

void multiplier(int *p, double *u0, int *m, double *u, int *n, double *b,
		double *influ, double *denom, int *N, double *s0, int *verbose);

void cramer_vonMises_Pickands(int n, int m, double *S,
			      double *T, double *Atheta,
			      double *stat);

void cramer_vonMises_CFG(int n, int m, double *S,
			 double *T, double *Atheta,
			 double *stat);

void cramer_vonMises_Afun(int *n, int *m, double *S,
			  double *T, double *Atheta,
			  double *stat, int *CFG);

#endif
