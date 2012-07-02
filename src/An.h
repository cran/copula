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
 * @file   An.h
 * @author Ivan Kojadinovic
 * @date   2009-2012
 *
 * @brief Rank-based versions of the Pickands and CFG estimators
 *        of the Pickands dependence function. Bivariate and
 *        multivariate versions.
 *
 */


#ifndef ANFUN_H
#define ANFUN_H

#define MIN(x,y) ((x)<(y)?(x):(y))
#define MAX(x,y) ((x)>(y)?(x):(y))

double biv_invAP(int n, double *S, double *T, double t);
double biv_logACFG(int n, double *S, double *T, double t);

// called from R:
void biv_AP(int *n, double *S, double *T, double *t, int *m,
		int *corrected, double *A);

void biv_ACFG(int *n, double *S, double *T, double *t, int *m,
	   int *corrected, double *A);

void mult_A(double *U, int *n, int *d, double *w, int *m,
	    double *AP, double *ACFG, double *AHT);

#endif
