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


/*****************************************************************************

  Set function utilities

  Ivan Kojadinovic, May 2007

*****************************************************************************/

int card(int n);
double sum_binom(int n, int k);
void k_power_set(int *n, int *k, int *power_set);
void k_power_set_char(int *n, int *k, int *k_power_set, char **subset);
void natural2binary(int *n, double *sf, int *power_set, double *sf_out);

