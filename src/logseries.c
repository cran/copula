#include <math.h>
#include "R.h"
#include "Rmath.h"

/* Kemp A.W. 1981, Applied Statistics 30(3) pp249--253 */


/* Algorithm LS */
int rlogseries_LS ( double alpha ) {
  double t = - alpha / log(1 - alpha);
  double u = runif(0.0, 1.0), p = t;
  int x = 1;
  while (u > p) {
    u = u - p;
    x = x + 1;
    p = p * alpha * (x - 1) / x;
  }
  return (x);
}


/* Algorithm LK */
int rlogseries_LK ( double alpha ) {
  double h = log(1 - alpha);
  double u2 = runif(0.0, 1.0), x = 1.0, u1, q;
  if ( u2 > alpha ) return (int) x;
  u1 = runif(0.0, 1.0);
  q = 1 - exp(u1 * h);
  if ( u2 < q * q ) return (int) (1 + log(u2) / log(q));
  if ( u2 > q ) return 1;
  return 2;
}


/* R wrapper */
void rlogseries_R ( int *n,  double *alpha, int *val) {
  int i;
  double thres = 0.95;

  GetRNGstate();
  for (i = 0; i < *n; i++) {
    if ( alpha[i] < thres ) val[i] = rlogseries_LS ( alpha[i] );
    else val[i] = rlogseries_LK ( alpha[i] );
  }
  PutRNGstate();
}
