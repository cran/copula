## Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

require(copula)
## fCopulae: don't do on CRAN, and really "can not" suggest fCopulae
tryfCop <- TRUE # for interactive convenience
## when run as BATCH:
tryfCop <- nzchar(Sys.getenv("R_copula_check_fCop")) ||
    identical("true", unname(Sys.getenv("R_MM_PKG_CHECKING")))

if(tryfCop) { ## will only "work" if not "--as-cran"
    .r <- require
    tryfCop <- suppressWarnings(.r(fCopulae, quietly=TRUE)) }
tryfCop

## From system.file("test-tools-1.R", package="Matrix"):
## Make sure errors are signaled
assertError <- function(expr) {
    d.expr <- deparse(substitute(expr))
    t.res <- tryCatch(expr, error = function(e) e)
    if(!inherits(t.res, "error"))
	stop(d.expr, "\n\t did not give an error", call. = FALSE)
    invisible(t.res)
}
numTailIndexLower <- function(copula, u) {
  ## u is a vector approaching 0
  pCopula(cbind(u, u, deparse.level = 0), copula) / u
}

numTailIndexUpper <- function(copula, u) {
  # u is a vector approaching 1
  (1 - 2 * u + pCopula(cbind(u, u, deparse.level = 0), copula)) / (1 - u)
}

(u.0 <- sort(outer(c(1,2,5), 10^-(1:5)), decreasing=TRUE)[-(1:2)])
## 0.1, 0.05, 0.02, 0.01, ..... 1e-5
u.1 <- 1 - u.0

### Upper Tail Dependence ---------------------------

# R/Copula:
gumbC3  <- gumbelCopula(param= 3, dim = 2)
gumbC20 <- gumbelCopula(param=20, dim = 2)
gumbC40 <- gumbelCopula(param=40, dim = 2)

ut20 <- numTailIndexUpper(gumbC20, u.1)
(ut40 <- numTailIndexUpper(gumbC40, u.1))

stopifnot(
 all.equal(tailIndex(gumbC20)[["upper"]],
           numTailIndexUpper(gumbC20, 1 - 1e-7), tol=1e-8)
 ,
 all.equal(tailIndex(gumbC40)[["upper"]],
           numTailIndexUpper(gumbC40, 1 - 1e-7), tol=1e-8)
)

if(tryfCop) { ## Rmetrics
    C <- parchmCopula(u.1,u.1, alpha=40, type = "4", alternative = TRUE)
    stopifnot(all.equal(ut40,     (1-2*u.1+C)/(1-u.1),
                        check.attr=FALSE, tol= 1e-14))
}


### Lower Tail Dependence-------------------------

S <- cbind(u.0,u.0)
# R/Copula:
## C  <- pCopula(dim = 2, copula = gumbelCopula(param=20), S)
## (C1  <- C/u.0)
(lt20 <- numTailIndexLower(gumbC20, u.0))

if(tryfCop) { ## Rmetrics
    C <-  parchmCopula(S, alpha=20, type = "4", alternative = FALSE)
    stopifnot(all.equal(lt20, C/u.0, check.attr=FALSE, tol= 1e-14))
}

signif(numTailIndexLower(gumbC3, 10^-(5*(1:40))),   3)#--> 0
## but for large theta, the convergence (to 0) is *MUCH* slower:
signif(numTailIndexLower(gumbC20, 10^-(5*(1:40))),  3)


###-------------- Frank --------------------------
Frank2 <- frankCopula(param=2, dim = 2)
tailIndex(Frank2) # 0 0

## Upper and lower tail dependence
(tl <- numTailIndexLower(Frank2, u.0))
stopifnot(all.equal(tl, numTailIndexUpper(Frank2, u.1), tol=1e-10))

stopifnot(
  (tu1 <- numTailIndexUpper(Frank2, .99999)) < .00003
,
  all.equal(tu1, numTailIndexLower(Frank2, .00001), tol=1e-6)
,
  (tu2 <- numTailIndexUpper(Frank2, 1-1e-6)) < 3e-6
,
  all.equal(tu2, numTailIndexLower(Frank2, 1e-6), tol= 1e-4)
)



###-------------- Elliptic --------------------------

u2 <- cbind(u.0,u.1)

(t.7.3 <- tCopula(0.7, df=3, dim = 2))
(t.9.2 <- tCopula(0.9, df=2, dim = 2))

t.frac <- tCopula(0.9, df=2.5, dim = 2)
## fractional df  currently (must) *fail* for pCopula
assertError(pCopula(cbind(u.0,u.1), t.frac))

stopifnot( {
    ft <- dCopula(u2, t.frac)
    all.equal(ft, dCopula(u2[,2:1], t.frac), tol=1e-15)
 },
 !is.unsorted(ft)
 ,
 all.equal(tailIndex(t.7.3)[["upper"]],
           numTailIndexUpper(t.7.3, 1 - 1e-8), tol=1e-5)
 ,
 all.equal(tailIndex(t.9.2)[["upper"]],
           numTailIndexUpper(t.9.2, 1 - 1e-8), tol=1e-7)
 ,
 all.equal(tailIndex(t.7.3)[["lower"]],
           numTailIndexLower(t.7.3, 1e-8), tol=1e-5)
 ,
 all.equal(tailIndex(t.9.2)[["lower"]],
           numTailIndexLower(t.9.2, 1e-8), tol=1e-7)
)

(ut. <- numTailIndexUpper(t.7.3, u.1))

if(tryfCop && .r(fCopulae)) { ## Rmetrics
    p.fC <- pellipticalCopula(u = u.1, v = u.1, rho = 0.7, param = c(nu=3))
    p. <- pCopula(u = cbind(u.1, u.1), t.7.3)
    ## they are really not "so equal"
    stopifnot(
    all.equal(p.fC, p., check.attr=FALSE, tol= 0.002)
    )
}



cat('Time elapsed: ', proc.time(),'\n') # for ''statistical reasons''
