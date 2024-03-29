## ----message=FALSE------------------------------------------------------------
require(copula)
doExtras <- FALSE

## -----------------------------------------------------------------------------
##' @title [m]inus Log-Likelihood for Archimedean Copulas ("fast version")
##' @param theta parameter (length 1 for our current families)
##' @param acop Archimedean copula (of class "acopula")
##' @param u data matrix n x d
##' @param n.MC if > 0 MC is applied with sample size equal to n.MC; otherwise,
##'        the exact formula is used
##' @param ... potential further arguments, passed to <acop> @dacopula()
##' @return negative log-likelihood
##' @author Martin Maechler (Marius originally)
mLogL <- function(theta, acop, u, n.MC=0, ...) { # -(log-likelihood)
    -sum(acop@dacopula(u, theta, n.MC=n.MC, log=TRUE, ...))
}

## -----------------------------------------------------------------------------
##' @title Plotting the Negative Log-Likelihood for Archimedean Copulas
##' @param cop an outer_nacopula (currently with no children)
##' @param u n x d  data matrix
##' @param xlim x-range for curve() plotting
##' @param main title for curve()
##' @param XtrArgs a list of further arguments for mLogL()
##' @param ... further arguments for curve()
##' @return invisible()
##' @author Martin Maechler
curveLogL <- function(cop, u, xlim, main, XtrArgs=list(), ...) {
    unam <- deparse(substitute(u))
    stopifnot(is(cop, "outer_nacopula"), is.list(XtrArgs),
              (d <- ncol(u)) >= 2, d == dim(cop),
              length(cop@childCops) == 0# not yet *nested* A.copulas
              )
    acop <- cop@copula
    th. <- acop@theta # the true theta
    acop <- setTheta(acop, NA) # so it's clear, the true theta is not used below
    if(missing(main)) {
        tau. <- cop@copula@tau(th.)
        main <- substitute("Neg. Log Lik."~ -italic(l)(theta, UU) ~ TXT ~~
			   FUN(theta['*'] == Th) %=>% tau['*'] == Tau,
			   list(UU = unam,
				TXT= sprintf("; n=%d, d=%d;  A.cop",
					     nrow(u), d),
				FUN = acop@name,
				Th = format(th.,digits=3),
				Tau = format(tau., digits=3)))
    }
    r <- curve(do.call(Vectorize(mLogL, "theta"), c(list(x, acop, u), XtrArgs)),
               xlim=xlim, main=main,
               xlab = expression(theta),
               ylab = substitute(- log(L(theta, u, ~~ COP)), list(COP=acop@name)),
               ...)
    if(is.finite(th.))
        axis(1, at = th., labels=expression(theta["*"]),
             lwd=2, col="dark gray", tck = -1/30)
    else warning("non-finite cop@copula@theta = ", th.)
    axis(1, at = initOpt(acop@name),
         labels = FALSE, lwd = 2, col = 2, tck = 1/20)
    invisible(r)
}

## -----------------------------------------------------------------------------
op <- options("copula:verboseUsingRmpfr"=TRUE) # see when "Rmpfr" methods are chosen automatically

## -----------------------------------------------------------------------------
n <- 200
d <- 100
tau <- 0.2
theta <- copJoe@iTau(tau)
cop <- onacopulaL("Joe", list(theta,1:d))
theta

## ----fig.align="center", fig.width=7.5, fig.height=6--------------------------
set.seed(1)
U1 <- rnacopula(n,cop)
enacopula(U1, cop, "mle") # 1.432885 --  fine
th4 <- 1 + (1:4)/4
mL.tr <- c(-3558.5, -3734.4, -3299.5, -2505.)
mLt1 <- sapply(th4, function(th) mLogL(th, cop@copula, U1, method="log.poly")) # default
mLt2 <- sapply(th4, function(th) mLogL(th, cop@copula, U1, method="log1p"))
mLt3 <- sapply(th4, function(th) mLogL(th, cop@copula, U1, method="poly"))
stopifnot(all.equal(mLt1, mL.tr, tolerance=5e-5),
          all.equal(mLt2, mL.tr, tolerance=5e-5),
          all.equal(mLt3, mL.tr, tolerance=5e-5))
system.time(r1l  <- curveLogL(cop, U1, c(1, 2.5), X=list(method="log.poly")))
mtext("all three polyJ() methods on top of each other")
system.time({
    r1J <- curveLogL(cop, U1, c(1, 2.5), X=list(method="poly"),
                     add=TRUE, col=adjustcolor("red", .4))
    r1m  <- curveLogL(cop, U1, c(1, 2.5), X=list(method="log1p"),
                      add=TRUE, col=adjustcolor("blue",.5))
})

## ----fig.align="center", fig.width=7.5, fig.height=6--------------------------
U2 <- rnacopula(n,cop)
summary(dCopula(U2, cop)) # => density for the *correct* parameter looks okay
## hmm: max = 5.5e177
if(doExtras)
    system.time(r2 <- curveLogL(cop, U2, c(1, 2.5)))
stopifnot(all.equal(enacopula(U2, cop, "mle"), 1.43992755, tolerance=1e-5),
          all.equal(mLogL(1.8, cop@copula, U2), -4070.1953,tolerance=1e-5)) # (was -Inf)

## ----fig.align="center", fig.width=7.5, fig.height=6--------------------------
U3 <- rnacopula(n,cop)
(th. <- enacopula(U3, cop, "mle")) # 1.4495
system.time(r3 <- curveLogL(cop, U3, c(1, 2.5)))
axis(1, at = th., label = quote(hat(theta)))

## ----fig.align="center", fig.width=7.5, fig.height=6--------------------------
U4 <- rnacopula(n,cop)
enacopula(U4, cop, "mle") # 1.4519  (prev. was 2.351 : "completely wrong")
summary(dCopula(U4, cop)) # ok (had one Inf)
if(doExtras)
    system.time(r4 <- curveLogL(cop, U4, c(1, 2.5)))
mLogL(2.2351, cop@copula, U4)
mLogL(1.5,    cop@copula, U4)
mLogL(1.2,    cop@copula, U4)
if(doExtras) # each curve takes almost 2 sec
    system.time({
        curveLogL(cop, U4, c(1, 1.01))
        curveLogL(cop, U4, c(1, 1.0001))
        curveLogL(cop, U4, c(1, 1.000001))
    })
## --> limit goes *VERY* steeply up to  0
## --> theta 1.164 is about the boundary:
stopifnot(identical(setTheta(cop, 1.164), onacopula(cop@copula, C(1.164, 1:100))),
	  all.equal(600.59577,
		    cop@copula@dacopula(U4[118,,drop=FALSE],
					theta=1.164, log = TRUE), tolerance=1e-5)) # was "Inf"

## -----------------------------------------------------------------------------
n <- 200
d <- 150
tau <- 0.3
(theta <- copJoe@iTau(tau))
cop <- onacopulaL("Joe",list(theta,1:d))

## ----fig.align="center", fig.width=7.5, fig.height=6--------------------------
set.seed(47)
U. <- rnacopula(n,cop)
enacopula(U., cop, "mle") # 1.784578
system.time(r. <- curveLogL(cop, U., c(1.1, 3)))
## => still looks very good

## -----------------------------------------------------------------------------
d <- 180
tau <- 0.4
(theta <- copJoe@iTau(tau))
cop <- onacopulaL("Joe",list(theta,1:d))

## ----fig.align="center", fig.width=7.5, fig.height=6--------------------------
U. <- rnacopula(n,cop)
enacopula(U., cop, "mle") # 2.217582
if(doExtras)
system.time(r. <- curveLogL(cop, U., c(1.1, 4)))
## => still looks very good

## -----------------------------------------------------------------------------
n <- 200
d <- 50 # smaller 'd' -- so as to not need 'Rmpfr' here
tau <- 0.2
(theta <- copGumbel@iTau(tau))
cop <- onacopulaL("Gumbel",list(theta,1:d))

## ----fig.align="center", fig.width=7.5, fig.height=6--------------------------
set.seed(1)
U1 <- rnacopula(n,cop)
if(doExtras) {
    U2 <- rnacopula(n,cop)
    U3 <- rnacopula(n,cop)
}
enacopula(U1, cop, "mle") # 1.227659 (was 1.241927)
##--> Plots with "many" likelihood evaluations
system.time(r1 <- curveLogL(cop, U1, c(1, 2.1)))
if(doExtras) system.time({
    mtext("and two other generated samples")
    r2 <- curveLogL(cop, U2, c(1, 2.1), add=TRUE)
    r3 <- curveLogL(cop, U3, c(1, 2.1), add=TRUE)
})

## -----------------------------------------------------------------------------
d <- 150
tau <- 0.6
(theta <- copGumbel@iTau(tau))
cG.5 <- onacopulaL("Gumbel",list(theta,1:d))

## -----------------------------------------------------------------------------
set.seed(17)
U4 <- rnacopula(n,cG.5)
U5 <- rnacopula(n,cG.5)
U6 <- rnacopula(n,cG.5)
if(doExtras) { ## "Rmpfr" is used {2012-06-21}: -- therefore about 18 seconds!
 tol <- if(interactive()) 1e-12 else 1e-8
 print(system.time(
 ee. <- c(enacopula(U4, cG.5, "mle", tol=tol),
          enacopula(U5, cG.5, "mle", tol=tol),
          enacopula(U6, cG.5, "mle", tol=tol))))
dput(ee.)# in case the following fails
## tol=1e-12 Linux nb-mm3 3.2.0-25-generic x86_64 (2012-06-23):
##   c(2.47567251789004, 2.48424484287686, 2.50410767129408)
##   c(2.475672518,      2.484244763,      2.504107671),
stopifnot(all.equal(ee., c(2.475672518, 2.484244763, 2.504107671),
		    tolerance= max(1e-7, 16*tol)))
}
## --> Plots with "many" likelihood evaluations
th. <- seq(1, 3, by= 1/4)
if(doExtras) # "default2012" (polyG default) partly uses Rmpfr here:
system.time(r4   <- sapply(th., mLogL, acop=cG.5@copula, u=U4))## 25.6 sec
## whereas this (polyG method) is very fast {and still ok}:
system.time(r4.p <- sapply(th., mLogL, acop=cG.5@copula, u=U4, method="pois"))
r4. <- c(0, -18375.33, -21948.033, -24294.995, -25775.502,
         -26562.609, -26772.767, -26490.809, -25781.224)
stopifnot(!doExtras ||
          all.equal(r4,   r4., tolerance = 8e-8),
          all.equal(r4.p, r4., tolerance = 8e-8))
## --> use fast method here as well:
system.time(r5.p <- sapply(th., mLogL, acop=cG.5@copula, u=U5, method="pois"))
system.time(r6.p <- sapply(th., mLogL, acop=cG.5@copula, u=U6, method="pois"))
if(doExtras) {
    if(FALSE) # for speed analysis, etc
        debug(copula:::polyG)
    mLogL(1.65, cG.5@copula, U4) # -23472.96
}
dd <- dCopula(U4, setTheta(cG.5, 1.64), log = TRUE,
              method = if(doExtras)"default" else "pois")
summary(dd)
stopifnot(!is.na(dd), # no NaN's anymore
	  40 < range(dd), range(dd) < 710)

## -----------------------------------------------------------------------------
n <- 64
d <- 5
tau <- 0.8
(theta <- copFrank@iTau(tau))
cop <- onacopulaL("Frank", list(theta,1:d))

## ----fig.align="center", fig.width=7.5, fig.height=6--------------------------
set.seed(11) # these seeds give no problems: 101, 41, 21
U. <- rnacopula(n,cop)
cop@copula <- setTheta(cop@copula, NA) # forget the true theta
system.time(f.ML <- emle(U., cop)); f.ML # --> fine: theta = 18.033, Log-lik = 314.01
if(doExtras)
    system.time(f.mlMC <- emle(U., cop, n.MC = 1e4)) # with MC
stopifnot(all.equal(unname(coef(f.ML)), 18.03331, tolerance= 1e-6),
	  all.equal(f.ML@min, -314.0143, tolerance=1e-6),
	  !doExtras || ## Simulate MLE (= SMLE) is "extra" random,  hmm...
	  all.equal(unname(coef(f.mlMC)), 17.8, tolerance= 0.01)
	  ##		   64-bit ubuntu: 17.817523
	  ##		 ? 64-bit Mac:	  17.741
	 )

cop@copula <- setTheta(cop@copula, theta)
r. <- curveLogL(cop, U., c(1, 200)) # => now looks fine
tail(as.data.frame(r.), 15)
stopifnot( is.finite( r.$y ),
	  ## and is convex (everywhere):
	  diff(r.$y, d=2) > 0)
options(op) # revert to previous state

## ----echo=FALSE---------------------------------------------------------------
print(sessionInfo(), locale=FALSE)

