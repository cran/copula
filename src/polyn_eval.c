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


#include "nacopula.h"

/**
 * Polynomial evaluation via Horner scheme
 *
 * @param coef coefficients
 * @param x evaluation point(s)
 * @return values
 * @author Marius Hofert
 */
SEXP polyn_eval(SEXP coef, SEXP x)
{
 SEXP result;
 R_xlen_t n = XLENGTH(x);
 int m = LENGTH(coef);
 // deal with integer or numeric -- NULL cannot (yet?) be coerced
 if(isNull(x)) { result = allocVector(REALSXP, 0); return result; }
 if(!isNull(coef))
     coef = isReal(coef) ? Rf_duplicate(coef) : coerceVector(coef, REALSXP);
 PROTECT(coef);
 PROTECT(x = isReal(x) ? Rf_duplicate(x) : coerceVector(x, REALSXP));
 PROTECT(result = Rf_duplicate(x));
 double *cf = REAL(coef), *xx = REAL(x), *res = REAL(result);
 for(R_xlen_t i = 0; i < n; i++) {
   double r, xi = xx[i];
   if(m == 0) {
      r = 0.;
   } else {
      int j = m-1;
      r = cf[j];
      for (j--; j >= 0; j--)
          r = cf[j] + r * xi;
   }
   res[i] = r;
 }
 UNPROTECT(3);
 return result;
}
