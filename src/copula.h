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


#ifndef COPULA_DEFS_H
#define COPULA_DEFS_H

#include <R.h>

#include "An.h"
#include "empcop.h"
#include "gof.h"
#include "set_utils.h"
#include "indepTests.h"

// ./logseries.c: __FIXME__ also have rLog_vec_c()  from nacopula
void rlogseries_R      (int *n, double *alpha, int *val);
void rlogseries_R_ln1p (int *n, double *h,  double *val);


// ./fgm.c:
void validity_fgm(int *p, double *alpha, int *valid);
void rfgm(int *p, double *alpha, int *n, double *x);

// ./evtest.c:
void evtest(double *U, int *n, int *p, double *g, int *m,
	    int *N, double *tg, int *nt, double *s0, int *der2n,
	    double *o, double *stat);

void evtestA(double *U, double *V, int *n, double *u, double *v,
	     int *m, int *CFG, int *N, double *s0);

void evtestA_derA(double *U, double *V, int *n,
		  double *u, double *v, int *m, int *CFG_tr__etc, int *N, double *s0);

void evtestA_stat(double *U, double *V, int *n, double *u, double *v, int *m,
		  int *CFG, double *stat, double *offset);

// ./exchtest.c:
void evsymtest(double *U, double *V, int *n, double *t, int *m,
	       int *CFG, int *N, double *s0);

void evsymtest_derA(double *U, double *V, int *n, double *t, int *m,
		    int *CFG, int *N, double *s0);

void evsymtest_stat(double *S, double *T, int *n, double *t, int *m,
		    int *CFG, double *stat);

void exchtestCn(double *U, double *V, int *n, double *u, double *v,
		int *m, int *N, double *s0);

void exchtestCn_stat(double *U, double *V, int *n, double *u, double *v,
		     int *m, double *stat);

void radsymtestCn_stat(double *U, int *n, int *p, double *V, int *m,
		       double *stat);

// R_debye.c:
// "_C": nameclash - already have R level 'debye_1'
void debye_1_C(double *x, int *len, double *val, double *err, int *status);
void debye_2_C(double *x, int *len, double *val, double *err, int *status);
void debye_3(double *x, int *len, double *val, double *err, int *status);
void debye_4(double *x, int *len, double *val, double *err, int *status);


#endif
