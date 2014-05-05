/**
 * @file   copula_int.h
 * @author Martin Maechler
 * @date   March 2014
 *
 * @brief R support for internationalized messages
 *        which might be translated { -> gettext() }
 */

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("copula", String)
#else
#define _(String) (String)
#define dngettext(pkg, String, StringP, N) (N > 1 ? StringP : String)
#endif
