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
 * @file   indepTests.h
 * @author Ivan Kojadinovic
 * @date   December 2007
 *
 * @brief  Multivariate and "vectorial" tests of independence and serial
 *         independence based on the empirical copula process
 *
 */

#ifndef INDEPTESTS_H
#define INDEPTESTS_H

#include "copula_int.h"

// indepTest_utils.c  --- temporary arrays -------------------------------------
void K_array(int n, int p, const double J[], double *K);
void L_array(int n, int p, const double K[], double *L);

double I_n(int n, int p, double *J, double *K, double *L);
double M_A_n(int n, int p, double *J, double *K, double *L, int A);

// different versions of temporary array J -----------------------
void J_m (int n, int p, const int b[], const double U[], const int R[], double *J);
void J_s (int n, int p,                const double U[],                double *J);

void J_sm(int n, int p, int q, const double U[], const int B[], double *J);
void J_u (int n, int p,        const double R[],                double *J);

// utility function: text progress bar 
void progressBar(int k, int N, int w);

// multIndepTest.c  --- Independence test among random vectors ------------------
void bootstrap_MA_I(int *n, int *N, int *p, int *b, double *U, int *m,
		    double *MA0, double *I0, int *subset, char **subset_char,
		    int *verbose);
void empirical_copula_test_rv(double *U, int *n, int *p, int *b, int *m, double *MA0,
			      double *I0, int *N, int *subset, double *MA, double *I,
			      double *pval, double *fisher, double *tippett, double *Ipval);

// serialIndepTest.c  --- Serial independence test of Genest and Remillard (2004) ---------
void simulate_empirical_copula_serial(int *n, int *N, int *p, int *m,
				      double *TA0, double *G0, int *subset,
				      char **subset_char, double *fisher0,
				      double *tippett0, int *verbose);
void empirical_copula_test_serial(double *U, int *n, int *p, int *m, double *TA0, double *G0,
				  int *N, int *subset, double *TA, double *G, double *pval,
				  double *fisher, double *tippett, double *globpval,
				  double *fisher0, double *tippett0);


// multSerialIndepTest.c  --- Multivariate serial independence test ------------------------
void bootstrap_serial(int *n, int *N, int *p, int *q, double *U, int *m,
		      double *MA0, double *I0, int *subset, char **subset_char,
		      int *verbose);
void empirical_copula_test_rv_serial(double *U, int *n, int *p, int *q, int *m, double *MA0,
				     double *I0, int *N, int *subset, double *MA, double *I,
				     double *pval, double *fisher, double *tippett, double *Ipval);

// ecIndepTest.c  --- Multivariate independence test of Genest and Remillard (2004) ------------
void simulate_empirical_copula(int *n, int *N, int *p, int *m, double *TA0,
			       double *G0, int *subset, char **subset_char,
			       double *fisher0, double *tippett0, int *verbose);
void empirical_copula_test(double *R, int *n, int *p, int *m, double *TA0, double *G0,
			   int *N, int *subset, double *TA, double *G, double *pval,
			   double *fisher, double *tippett, double *globpval,
			   double *fisher0, double *tippett0);

#endif


