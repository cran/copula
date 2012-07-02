#include <R.h>

#ifndef __MYGSL_H__
#define __MYGSL_H__

/* from gsl_version.h */
#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

enum {
  GSL_SUCCESS  = 0,
  GSL_FAILURE  = -1,
  GSL_CONTINUE = -2,  /* iteration has not converged */
  GSL_EDOM     = 1,   /* input domain error, e.g sqrt(-1) */
  GSL_ERANGE   = 2,   /* output range error, e.g. exp(1e100) */
  GSL_EFAULT   = 3,   /* invalid pointer */
  GSL_EINVAL   = 4,   /* invalid argument supplied by user */
  GSL_EFAILED  = 5,   /* generic failure */
  GSL_EFACTOR  = 6,   /* factorization failed */
  GSL_ESANITY  = 7,   /* sanity check failed - shouldn't happen */
  GSL_ENOMEM   = 8,   /* malloc failed */
  GSL_EBADFUNC = 9,   /* problem with user-supplied function */
  GSL_ERUNAWAY = 10,  /* iterative process is out of control */
  GSL_EMAXITER = 11,  /* exceeded max number of iterations */
  GSL_EZERODIV = 12,  /* tried to divide by zero */
  GSL_EBADTOL  = 13,  /* user specified an invalid tolerance */
  GSL_ETOL     = 14,  /* failed to reach the specified tolerance */
  GSL_EUNDRFLW = 15,  /* underflow */
  GSL_EOVRFLW  = 16,  /* overflow  */
  GSL_ELOSS    = 17,  /* loss of accuracy */
  GSL_EROUND   = 18,  /* failed because of roundoff error */
  GSL_EBADLEN  = 19,  /* matrix, vector lengths are not conformant */
  GSL_ENOTSQR  = 20,  /* matrix not square */
  GSL_ESING    = 21,  /* apparent singularity detected */
  GSL_EDIVERGE = 22,  /* integral or series is divergent */
  GSL_EUNSUP   = 23,  /* requested feature is not supported by the hardware */
  GSL_EUNIMPL  = 24,  /* requested feature not (yet) implemented */
  GSL_ECACHE   = 25,  /* cache limit exceeded */
  GSL_ETABLE   = 26,  /* table limit exceeded */
  GSL_ENOPROG  = 27,  /* iteration is not making progress towards solution */
  GSL_ENOPROGJ = 28,  /* jacobian evaluations are not improving the solution */
  GSL_ETOLF    = 29,  /* cannot reach the specified tolerance in F */
  GSL_ETOLX    = 30,  /* cannot reach the specified tolerance in X */
  GSL_ETOLG    = 31,  /* cannot reach the specified tolerance in gradient */
  GSL_EOF      = 32   /* end of file */
} ;


/* define gsl_sf_result; from specfunc/gsl_sf_result.h */
struct gsl_sf_result_struct {
  double val;
  double err;
};
typedef struct gsl_sf_result_struct gsl_sf_result;

/* check for underflow; from specfunc/check.h */
#define CHECK_UNDERFLOW(r) if (fabs((r)->val) < GSL_DBL_MIN) GSL_ERROR("underflow", GSL_EUNDRFLW);

/* my own definition of GSL_ERROR */
/* GSL_ERROR: call the error handler, and return the error code */
#define GSL_ERROR(reason, gsl_errno) \
       do { \
       Rprintf("ERROR: %d\n", gsl_errno); \
       return gsl_errno ; \
       } while (0)

/* from specfunc/error.h */
#define DOMAIN_ERROR(result) do { (result)->val = GSL_NAN; (result)->err = GSL_NAN; GSL_ERROR ("domain error", GSL_EDOM); } while(0)

/* from gsl_nan.h and gsl_div.c infnan.c */
#ifdef INFINITY
# define GSL_POSINF INFINITY
# define GSL_NEGINF (-INFINITY)
#elif defined(HUGE_VAL)
# define GSL_POSINF HUGE_VAL
# define GSL_NEGINF (-HUGE_VAL)
#else
# define GSL_POSINF (gsl_posinf())
# define GSL_NEGINF (gsl_neginf())
#endif

#ifdef NAN
# define GSL_NAN NAN
#elif defined(INFINITY)
# define GSL_NAN (INFINITY/INFINITY)
#else
# define GSL_NAN (gsl_nan())
#endif

double gsl_fdiv (const double x, const double y);

double gsl_nan (void);


/* from gsl_sf_debye.h */
int     gsl_sf_debye_1_e(const double x, gsl_sf_result * result);
int     gsl_sf_debye_2_e(const double x, gsl_sf_result * result);
int     gsl_sf_debye_3_e(const double x, gsl_sf_result * result);
int     gsl_sf_debye_4_e(const double x, gsl_sf_result * result);

#include "gsl_machine.h"

__END_DECLS


#endif /* __MYGSL_H__ */
