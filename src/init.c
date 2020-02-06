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


#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "copula.h"
#include "nacopula.h"

#define CDEF(name)  {#name, (DL_FUNC) &name, sizeof(name ## _t)/sizeof(name ## _t[0]), name ##_t}

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}


// ./An.c ///////////////////////////////////////////////////////////////////

static R_NativePrimitiveArgType biv_ACFG_t[7] = {
    INTSXP, REALSXP, REALSXP, REALSXP, /* m: */ INTSXP, INTSXP, REALSXP
};
static R_NativePrimitiveArgType biv_AP_t[7] = {
    INTSXP, REALSXP, REALSXP, REALSXP, /* m: */ INTSXP, INTSXP, REALSXP
};

static R_NativePrimitiveArgType mult_A_t[8] = {
    REALSXP, INTSXP, INTSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP
};

// ./empcop.c ///////////////////////////////////////////////////////////////////
static R_NativePrimitiveArgType Cn_C_t[8] = {
	REALSXP, INTSXP, INTSXP, REALSXP, INTSXP, REALSXP, REALSXP, INTSXP
};

// ./fgm.c /////////////////////////////////////////////////////////////////////

static R_NativePrimitiveArgType validity_fgm_t[3] = { INTSXP, REALSXP, INTSXP };
static R_NativePrimitiveArgType rfgm_t[4] = { INTSXP, REALSXP, INTSXP, REALSXP };


// ./set_utils.c ///////////////////////////////////////////////////////////////

static R_NativePrimitiveArgType k_power_set_t[3] = { INTSXP, INTSXP, INTSXP };
static R_NativePrimitiveArgType k_power_set_char_t[4] = {
    INTSXP, INTSXP, INTSXP, STRSXP };


// ./gof.c & gof.h /////////////////////////////////////////////////////////////

static R_NativePrimitiveArgType cramer_vonMises_t[5] = {
    INTSXP, INTSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType cramer_vonMises_grid_t[7] = {
    INTSXP, REALSXP, INTSXP, REALSXP, INTSXP, REALSXP, REALSXP };
static R_NativePrimitiveArgType multiplier_t[11] = {
    INTSXP, REALSXP, INTSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP,
    INTSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType cramer_vonMises_Pickands_t[6] = {
    INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP };
static R_NativePrimitiveArgType cramer_vonMises_CFG_t[6] = {
    INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP };
// MM: FIXME *_Afun is superfluous
static R_NativePrimitiveArgType cramer_vonMises_Afun_t[7] = {
    INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP };


// ./logseries.c ///////////////////////////////////////////////////////////////

static R_NativePrimitiveArgType rlogseries_R_t     [3] = { INTSXP, REALSXP, INTSXP };
static R_NativePrimitiveArgType rlogseries_R_ln1p_t[3] = { INTSXP, REALSXP, REALSXP };


// ./evtest.c //////////////////////////////////////////////////////////////////

static R_NativePrimitiveArgType evtest_t[] = {
    REALSXP, INTSXP, INTSXP, REALSXP, INTSXP,
    INTSXP, REALSXP, INTSXP, REALSXP, INTSXP,
    REALSXP, REALSXP
};
static R_NativePrimitiveArgType evtestA_t[] = {REALSXP, REALSXP, INTSXP, REALSXP, REALSXP,
	     INTSXP, INTSXP, INTSXP, REALSXP};

static R_NativePrimitiveArgType evtestA_derA_t[] = {
    REALSXP, REALSXP, INTSXP,
    REALSXP, REALSXP, INTSXP,
    INTSXP, INTSXP, REALSXP
};
static R_NativePrimitiveArgType evtestA_stat_t[] = {
    REALSXP, REALSXP, INTSXP,
    REALSXP, REALSXP, INTSXP,
    INTSXP, REALSXP, REALSXP
};

// ./exchtest.c ////////////////////////////////////////////////////////////////

static R_NativePrimitiveArgType evsymtest_t[] = {REALSXP, REALSXP, INTSXP, REALSXP, INTSXP,
	       INTSXP, INTSXP, REALSXP};

static R_NativePrimitiveArgType evsymtest_derA_t[] = {REALSXP, REALSXP, INTSXP, REALSXP, INTSXP,
		    INTSXP, INTSXP, REALSXP};

static R_NativePrimitiveArgType evsymtest_stat_t[] = {REALSXP, REALSXP, INTSXP, REALSXP, INTSXP,
		    INTSXP, REALSXP};

static R_NativePrimitiveArgType exchtestCn_t[] = {REALSXP, REALSXP, INTSXP, REALSXP, REALSXP,
		INTSXP, INTSXP, REALSXP};

static R_NativePrimitiveArgType exchtestCn_stat_t[] = {REALSXP, REALSXP, INTSXP, REALSXP, REALSXP,
		     INTSXP, REALSXP};

static R_NativePrimitiveArgType radsymtestCn_stat_t[] = {REALSXP, INTSXP, INTSXP, REALSXP, INTSXP,
							 REALSXP};

// ./multIndepTest.c ////////////////////////////////////////////////////////////

static R_NativePrimitiveArgType bootstrap_MA_I_t[] = {INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP,
	       REALSXP, REALSXP, INTSXP, STRSXP,
	       INTSXP};
static R_NativePrimitiveArgType empirical_copula_test_rv_t[] = {REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP,
			      REALSXP, INTSXP, INTSXP, REALSXP, REALSXP,
			      REALSXP, REALSXP, REALSXP, REALSXP};


// ./serialIndepTest.c ////////////////////////////////////////////////////////////

static R_NativePrimitiveArgType simulate_empirical_copula_serial_t[] = {INTSXP, INTSXP, INTSXP, INTSXP,
				      REALSXP, REALSXP, INTSXP,
				      STRSXP, REALSXP,
				      REALSXP, INTSXP};
static R_NativePrimitiveArgType empirical_copula_test_serial_t[] = {REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP,
				  INTSXP, INTSXP, REALSXP, REALSXP, REALSXP,
				  REALSXP, REALSXP, REALSXP,
				  REALSXP, REALSXP};


// ./multSerialIndepTest.c ///////////////////////////////////////////////////////////

static R_NativePrimitiveArgType bootstrap_serial_t[] = {INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP,
		      REALSXP, REALSXP, INTSXP, STRSXP,
		      INTSXP};
static R_NativePrimitiveArgType empirical_copula_test_rv_serial_t[] = {REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP,
				     REALSXP, INTSXP, INTSXP, REALSXP, REALSXP,
				     REALSXP, REALSXP, REALSXP, REALSXP};


// ./ecIndepTest.c ////////////////////////////////////////////////////////////

static R_NativePrimitiveArgType simulate_empirical_copula_t[] = {INTSXP, INTSXP, INTSXP, INTSXP, REALSXP,
			       REALSXP, INTSXP, STRSXP,
			       REALSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType empirical_copula_test_t[] = {REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP,
			   INTSXP, INTSXP, REALSXP, REALSXP, REALSXP,
			   REALSXP, REALSXP, REALSXP,
			   REALSXP, REALSXP};


// ./R_debye.c /////////////////////////////////////////////////////////////////

/* no longer needed; now using debye_1 from package gsl
static R_NativePrimitiveArgType debye_1_C_t[5] = { REALSXP, INTSXP, REALSXP,REALSXP, INTSXP };
static R_NativePrimitiveArgType debye_2_C_t[5] = { REALSXP, INTSXP, REALSXP,REALSXP, INTSXP };
static R_NativePrimitiveArgType debye_3_t[5] = { REALSXP, INTSXP, REALSXP,REALSXP, INTSXP };
static R_NativePrimitiveArgType debye_4_t[5] = { REALSXP, INTSXP, REALSXP,REALSXP, INTSXP };
*/

static const R_CMethodDef CEntries[]  = {
    CDEF(biv_ACFG),
    CDEF(biv_AP),
    CDEF(mult_A),
    CDEF(Cn_C),

    CDEF(validity_fgm),
    CDEF(rfgm),
    CDEF(k_power_set),
    CDEF(k_power_set_char),

    CDEF(cramer_vonMises),
    CDEF(cramer_vonMises_grid),
    CDEF(multiplier),
    CDEF(cramer_vonMises_Pickands),
    CDEF(cramer_vonMises_CFG),
    CDEF(cramer_vonMises_Afun),

    CDEF(rlogseries_R),
    CDEF(rlogseries_R_ln1p),

    CDEF(evtest),
    CDEF(evtestA),
    CDEF(evtestA_derA),
    CDEF(evtestA_stat),
    CDEF(evsymtest),
    CDEF(evsymtest_derA),
    CDEF(evsymtest_stat),
    CDEF(exchtestCn),
    CDEF(exchtestCn_stat),
    CDEF(radsymtestCn_stat),

    CDEF(bootstrap_MA_I),
    CDEF(empirical_copula_test_rv),
    CDEF(simulate_empirical_copula_serial),
    CDEF(empirical_copula_test_serial),
    CDEF(bootstrap_serial),
    CDEF(empirical_copula_test_rv_serial),
    CDEF(simulate_empirical_copula),
    CDEF(empirical_copula_test),
    /*
    CDEF(debye_1_C),
    CDEF(debye_2_C),
    CDEF(debye_3),
    CDEF(debye_4),
    */
    {NULL, NULL, 0}
};

static R_CallMethodDef CallEntries[] = {
    CALLDEF(sinc_c, 1),
    CALLDEF(A__c, 3),
    CALLDEF(polyn_eval, 2),

    CALLDEF(rstable_c, 2),
    CALLDEF(retstable_c, 4),

    CALLDEF(rLog_vec_c, 3),
    CALLDEF(rSibuya_vec_c, 2),

    CALLDEF(rF01Frank_vec_c, 5),
    CALLDEF(rF01Joe_vec_c, 3),

    CALLDEF(gofT2stat_c, 2),

    {NULL, NULL, 0}
};

/**
 * register routines
 * @param dll pointer
 * @return none
 * @author Martin Maechler
 */
void
#ifdef HAVE_VISIBILITY_ATTRIBUTE
__attribute__ ((visibility ("default")))
#endif
R_init_copula(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
