/*#################################################################
##   Copula R package by Jun Yan Copyright (C) 2007
##
##   Set function utilities
##   Copyright (C) 2007 Ivan Kojadinovic <ivan.kojadinovic@univ-nantes.fr>
##
##   This program is free software; you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation; either version 2 of the License, or
##   (at your option) any later version.
##
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License along
##   with this program; if not, write to the Free Software Foundation, Inc.,
##   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
##
#################################################################*/

/*****************************************************************************

  Set function utilities

  Ivan Kojadinovic, May 2007

*****************************************************************************/

int card(int n);
double sum_binom(int n, int k);
void k_power_set(int *n, int *k, int *power_set);
void k_power_set_char(int *n, int *k, int *k_power_set, char **subset);
void natural2binary(int *n, double *sf, int *power_set, double *sf_out);

