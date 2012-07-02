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

## NB: run from ../tests/gof-ex.R --> keep "robust" (and quick) !
##              ~~~~~~~~~~~~~~~~~  (MM: see ../tests/gof-ex.Rout.save )

require(copula)

##' @title [m]inus Log Likelihood for Archimedean copula "fast version"
##' @param theta parameter (length 1 for our current families)
##' @param acop Archimedean copula (of class "acopula")
##' @param u data matrix n x d
##' @param n.MC if > 0 MC is applied with sample size equal to n.MC; otherwise,
##'        the exact formula is used
##' @param ... potential further arguments, passed to <acop> @dacopula()
##' @return
##' @author Martin Maechler (Marius originally)
mLogL <- function(theta, acop, u, n.MC=0, ...) { # -(log-likelihood)
    -sum(acop@dacopula(u, theta, n.MC=n.MC, log=TRUE, ...))
}

##' @title
##' @param cop an outer_nacopula (currently with no children)
##' @param u n x d  data matrix
##' @param xlim x-range for curve() plotting
##' @param main title for curve()
##' @param XtrArgs a list of further arguments for mLogL()
##' @param ... further arguments for curve()
##' @return
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

(doExtras <- interactive() || nzchar(Sys.getenv("R_copula_check_extra")))
## Want to see when "Rmpfr" methods are chosen automatically:
options("copula:verboseUsingRmpfr" = TRUE)


### "Joe", tau = 0.2 ###########################################################

n <- 200
d <- 100
tau <- 0.2
(theta <- copJoe@tauInv(tau))# 1.44381
(cop <- onacopulaL("Joe",list(theta,1:d)))

set.seed(1)
U1 <- rnacopula(n,cop)
enacopula(U1, cop, "mle") # 1.432885 --  fine
(th4 <- 1 + (1:4)/4)
mL.tr <- c(-3558.5, -3734.4, -3299.5, -2505.)
mLt1 <- sapply(th4, function(th) mLogL(th, cop@copula, U1, method="log.poly")) # default
mLt2 <- sapply(th4, function(th) mLogL(th, cop@copula, U1, method="log1p"))
mLt3 <- sapply(th4, function(th) mLogL(th, cop@copula, U1, method="poly"))
stopifnot(all.equal(mLt1, mL.tr, tol=5e-5),
          all.equal(mLt2, mL.tr, tol=5e-5),
          all.equal(mLt3, mL.tr, tol=5e-5))

system.time(r1l  <- curveLogL(cop, U1, c(1, 2.5), X=list(method="log.poly")))
if(doExtras) {
mtext("all three polyJ() methods on top of each other")
system.time(r1J  <- curveLogL(cop, U1, c(1, 2.5), X=list(method="poly"),
                              add=TRUE, col=adjustcolor("red", .4)))
system.time(r1m  <- curveLogL(cop, U1, c(1, 2.5), X=list(method="log1p"),
                              add=TRUE, col=adjustcolor("blue",.5)))
}

U2 <- rnacopula(n,cop)
## the density for the *correct* parameter looks okay
summary(dnacopula(cop, U2))
## hmm:  max = 5.5e177
if(doExtras)
system.time(r2 <- curveLogL(cop, U2, c(1, 2.5)))
stopifnot(all.equal(enacopula(U2, cop, "mle"), 1.43991422),
          all.equal(mLogL(1.8, cop@copula, U2), -4070.14762))# (was -Inf)

U3 <- rnacopula(n,cop)
enacopula(U3, cop, "mle") # 1.4495
if(doExtras)
system.time(r3 <- curveLogL(cop, U3, c(1, 2.5)))

U4 <- rnacopula(n,cop)
enacopula(U4, cop, "mle") # 1.4519  was 2.351..  "completely wrong"
summary(dnacopula(cop, U4)) # ok (had one Inf)
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
##--> limit goes *VERY* steeply  up to  0

##--> theta 1.164 is about the boundary:
stopifnot(all.equal(600.59261959,
  cop@copula@dacopula(U4[118,], theta=1.164, log = TRUE)))## was "Inf"


### "Joe", harder cases: d = 150, tau = 0.3 ####################################

n <- 200
d <- 150
tau <- 0.3
(theta <- copJoe@tauInv(tau))# 1.772
(cop <- onacopulaL("Joe",list(theta,1:d)))
set.seed(47)
U. <- rnacopula(n,cop)
enacopula(U., cop, "mle") # 1.784578
system.time(r. <- curveLogL(cop, U., c(1.1, 3)))
## still looks very good

### "Joe", even harder: d = 180, tau = 0.4 #####################################

d <- 180
tau <- 0.4
(theta <- copJoe@tauInv(tau))# 2.219
(cop <- onacopulaL("Joe",list(theta,1:d)))
U. <- rnacopula(n,cop)
enacopula(U., cop, "mle") # 2.217582
if(doExtras)
system.time(r. <- curveLogL(cop, U., c(1.1, 4)))
## still looks very good


### Similar for Gumbel ########################################################

n <- 200
d <- 50 # smaller 'd' -- so as to not need 'Rmpfr' here
tau <- 0.2
(theta <- copGumbel@tauInv(tau))# 1.25
(cop <- onacopulaL("Gumbel",list(theta,1:d)))

set.seed(1)
U1 <- rnacopula(n,cop)
U2 <- rnacopula(n,cop)
U3 <- rnacopula(n,cop)

enacopula(U1, cop, "mle") # 1.227659 (was 1.241927)
##--> Plots with "many" likelihood evaluations
system.time(r1 <- curveLogL(cop, U1, c(1, 2.1)))
mtext("and two other generated samples")
system.time(r2 <- curveLogL(cop, U2, c(1, 2.1), add=TRUE))
system.time(r3 <- curveLogL(cop, U3, c(1, 2.1), add=TRUE))

### "Gumbel", harder: d = 150, tau = 0.6 #######################################

d <- 150
tau <- 0.6
(theta <- copGumbel@tauInv(tau))# 2.5
cG.5 <- onacopulaL("Gumbel",list(theta,1:d))

set.seed(17)
U4 <- rnacopula(n,cG.5)
U5 <- rnacopula(n,cG.5)
U6 <- rnacopula(n,cG.5)

if(doExtras) { ## "Rmpfr" is used {2012-06-21}: -- therefore about 18 seconds!
tol <- if(interactive()) 1e-12 else 1e-8
system.time(
 ee. <- c(enacopula(U4, cG.5, "mle", tol=tol),
          enacopula(U5, cG.5, "mle", tol=tol),
          enacopula(U6, cG.5, "mle", tol=tol)))
dput(ee.)# in case the following fails
## tol=1e-12 Linux nb-mm3 3.2.0-25-generic x86_64 (2012-06-23):
##   c(2.47567251789004, 2.48424484287686, 2.50410767129408)
##   c(2.475672518,      2.484244763,      2.504107671),
stopifnot(all.equal(ee., c(2.475672518, 2.484244763, 2.504107671),
		    tol= max(1e-7, 16*tol)))
}

##--> Plots with "many" likelihood evaluations
(th. <- seq(1, 3, by= 1/4))
if(doExtras)## "default2012" (polyG default) partly uses Rmpfr here:
system.time(r4   <- sapply(th., mLogL, acop=cG.5@copula, u=U4))## 25.6 sec
## whereas this (polyG method) is very fast {and still ok}:
system.time(r4.p <- sapply(th., mLogL, acop=cG.5@copula, u=U4, method="pois"))
r4. <- c(0, -18375.33, -21948.033, -24294.995, -25775.502,
         -26562.609, -26772.767, -26490.809, -25781.224)
stopifnot(!doExtras ||
          all.equal(r4,   r4., tol = 8e-8),
          all.equal(r4.p, r4., tol = 8e-8))
##--> use fast method here as well:
system.time(r5.p <- sapply(th., mLogL, acop=cG.5@copula, u=U5, method="pois"))
system.time(r6.p <- sapply(th., mLogL, acop=cG.5@copula, u=U6, method="pois"))

if(FALSE) ## for speed analysis, etc
    debug(copula:::polyG)
mLogL(1.65, cG.5@copula, U4)# -23472.96

dd <- cG.5@copula@dacopula(U4, 1.64, log = TRUE,
			   method = if(doExtras)"default" else "pois")
summary(dd)
stopifnot(!is.na(dd), ## no NaN's anymore
	  40 < range(dd), range(dd) < 710)


### "Frank", a case we found "hard" ############################################

n <- 64
d <- 5
tau <- 0.8
(theta <- copFrank@tauInv(tau))# 18.192
(cop <- onacopulaL("Frank",list(theta,1:d)))
set.seed(11) ## these seeds give no problem: 101, 41, 21
U. <- rnacopula(n,cop)
## forget the true theta:
cop@copula <- setTheta(cop@copula, NA)
system.time(f.ML <- emle(U., cop)); f.ML ## --> fine: theta = 18.033, Log-lik = 314.01
if(doExtras)## with MC :
    system.time(f.mlMC <- emle(U., cop, n.MC = 1e4))## takes a while
    ## (7.5 sec on nb-mm3 2010)
stopifnot(
          all.equal(unname(coef(f.ML)), 18.03331, tol= 1e-6)
	  ,
	  all.equal(f.ML@min, -314.0143, tol=1e-6)
	  ,
	  !doExtras || ## Simulate MLE (= SMLE) is "extra" random,  hmm...
	  all.equal(unname(coef(f.mlMC)), 17.8, tol= 0.01)
	  ##		   64-bit ubuntu: 17.817523
	  ##		 ? 64-bit Mac:	  17.741
	 )

cop@copula <- setTheta(cop@copula, theta)# for the plot:
r. <- curveLogL(cop, U., c(1, 200))
## now looks fine
tail(as.data.frame(r.), 15)
stopifnot( is.finite( r.$y ),
	  ## and is convex (everywhere):
	  diff(r.$y, d=2) > 0)
