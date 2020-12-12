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


#ifndef NACOPULA_DEFS_H
#define NACOPULA_DEFS_H

#include <R.h>
#include <Rinternals.h>
#include "copula_int.h"

SEXP sinc_c(SEXP x_);
SEXP A__c(SEXP x_, SEXP alpha, SEXP I_alpha);

SEXP rstable_c(SEXP n, SEXP alpha);
SEXP retstable_c(SEXP V0_, SEXP h, SEXP alpha, SEXP method);

SEXP rLog_vec_c(SEXP n_, SEXP p_, SEXP Ip_);
SEXP rSibuya_vec_c(SEXP n_, SEXP alpha_);
SEXP rF01Frank_vec_c(SEXP V0_, SEXP theta_0_, SEXP theta_1_, SEXP rej_, SEXP approx_);
SEXP rF01Joe_vec_c(SEXP V0_, SEXP alpha_, SEXP approx_);

SEXP gofT2stat_c(SEXP u1_, SEXP u2_);

/**
 * C API---for "us" but maybe also other R packages
 * "export" it via ../inst/include/
*/
SEXP polyn_eval(SEXP coef, SEXP x);

double sinc_MM(double x);
double A_(double x, double alpha);
double BdB0(double x, double alpha);

// retstable.c :
double rstable0(double alpha);
double rstable (double alpha);
void rstable_vec(double S[], const R_xlen_t n, const double alpha);
void retstable_MH(double *St, const double V0[], double h, double alpha, R_xlen_t n);
void retstable_LD(double *St, const double V0[], double h, double alpha, R_xlen_t n);

double rLog(double p, double Ip);
double rSibuya(double alpha, double gamma_1_a /**< == Gamma(1 - alpha) */);
double rSibuya_sum(R_xlen_t n, double alpha, double gamma_1_a /**< == Gamma(1 - alpha) */);
void rSibuya_vec(double* V, R_xlen_t n, double alpha);
double rF01Frank(double V0, double theta0, double theta1, double p0, double p1,
	double gamma_1_a, double rej, int approx);
void rF01Frank_vec(double *V01, const double *V0, R_xlen_t n, double theta_0,
	double theta_1, double rej, int approx);
double rF01Joe(double V0, double alpha, double gamma_1_a /**< == Gamma(1 - alpha) */,
	int approx);
void rF01Joe_vec(double* V01, const double *V0, R_xlen_t n, double alpha, double approx);

#endif
