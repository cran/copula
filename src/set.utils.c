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


/*****************************************************************************

  Set function utilities

  Ivan Kojadinovic, May 2007

*****************************************************************************/

#include <R.h>
#include <Rmath.h>

#include "set.utils.h"

/*****************************************************************************

  Counts the number of bits equal to 1

*****************************************************************************/

int card(int n)
{
  int i;
  for(i=0; n; n >>= 1)
    i += n & 1;
  return(i);
}

/*****************************************************************************

  sum_{i=0}^k choose(n,i)

*****************************************************************************/

double sum_binom(int n, int k)
{
  int i;
  double s = 1.0;
  for (i=1;i<=k;i++)
    s += choose(n,i);
  return s;
}

/*****************************************************************************

  Generation of the first k + 1 levels of the power set of X
  in the "natural" order. Recursive function.

*****************************************************************************/

void k_power_set_rec(int n, int k, int last, int *power_set, int *b)
{
  int i, istart;

  /* look for the leftmost 1 in b and start to fill blank cases
     with 1 left from this position */
  istart = n;

  if(*b > 0)
    while(!(*b & 1<<(istart-1)))
      istart--;
  else
    istart = 0;

  for(i=istart; i<n; i++) {

    last++;
    *(power_set+last) = *b + (1<<i);
  }

  if(last != (int)sum_binom(n,k) - 1)
    k_power_set_rec(n, k, last, power_set, b+1);
}

/*****************************************************************************

  Generation of the first k + 1 levels of the power set of X
  in the "natural" order. Wrapps the previous function.

*****************************************************************************/

void k_power_set(int *n, int *k, int *power_set)
{
  power_set[0] = 0;
  k_power_set_rec(*n, *k, 0, power_set, power_set);
}

/*****************************************************************************

  Conversion binary -> list of elements
  b : binary code to be converted
  x : result in an array (ordered liste of elements)

*****************************************************************************/

void binary2subset(int n, int b, int *x)
{
  int i;
  for(i=0; i<n; i++)
    if(b & 1 << i) {

      *x = i;
      x++;
    }
}

/*****************************************************************************

  Converts the k power set of X in the "natural" order to char**.

*****************************************************************************/

#define SET_MAX 4

void k_power_set_char(int *n, int *sb, int *k_power_set, char **subset)
{
  int i, j;
  int x[32];
  char string[255];

  sprintf(subset[0],"{}");

  for(i=1; i<*sb; i++) {

    for(j=0; j<*n; j++)
      x[j]=0;

    binary2subset(*n,k_power_set[i],x);

    subset[i] = (char *) R_alloc(SET_MAX * (*n), sizeof(char));

    sprintf(subset[i],"{%d",x[0]+1);

    for(j=1; j<card(k_power_set[i]); j++) {

      sprintf(string,",%d", x[j]+1);
      strcat(subset[i],string);
    }

    strcat(subset[i],"}");
  }
}

/*****************************************************************************

  Writing a set function given in "natural" order in the binary order
  Pre: power_set contains the power_set in "natural" order

*****************************************************************************/

void natural2binary(int *n, double *sf, int *power_set, double *sf_out) {

  int i;
  for(i=0; i<(1<<*n); i++)
    sf_out[power_set[i]] = sf[i];
}

