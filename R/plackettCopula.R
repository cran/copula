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


plackettCopula <- function(param = NA_real_) {
    ## expressions for cdf and pdf: -- FIXME still used ??
    cdfE <- quote({ I <- 1 + (alpha - 1) * (u1 + u2)
      0.5 / (alpha - 1) * sqrt(I - (I^2 - 4 * alpha * (alpha - 1) * u1 * u2))
    })
    pdfE <- quote(((1 + (alpha - 1) * (u1 + u2))^2 - 4 * alpha * (alpha - 1) * u1 * u2)^(- 3/2) *
                  alpha * (1 + (alpha - 1) * (u1 + u2 - 2 * u1 * u2)))
    dim <- 2L
    new("plackettCopula",
        dimension = dim,
        parameters = param[1],
        exprdist = expression(cdf = cdfE, pdf = pdfE),## FIXME: still used ??
        param.names = "alpha",
        param.lowbnd = 0,
        param.upbnd = Inf,
        fullname = "<deprecated slot>")# "Plackett copula family"
}

pplackettCopula <- function(u, copula) {
  ## dim == 2
  u1 <- u[,1]
  u2 <- u[,2]
  if(is.na(alpha <- copula@parameters[1])) stop("parameter is NA")
  if(alpha == 1) return(u1 * u2) # independence copula
  ## Joe (1997, p.141)
  eta <- alpha - 1
  ## FIXME: better approximation for alpha = 1+eta ~= 1 (<=>  eta ~= 0)
  Ieuu <- 1 + eta * (u1 + u2)
  0.5 / eta * (Ieuu - sqrt(Ieuu^2 - 4 * alpha * eta * u1 * u2))
}

## a version for numerical exploration: 'alpha' instead copula:
.dplackettCopula <- function(u, alpha, log = FALSE) {
  ## dim == 2
  stopifnot(identical(ncol(u), 2L), length(alpha) == 1, !is.na(alpha))
  if(alpha == 1) ## independence Copula
    return(rep(if(log) 0 else 1, nrow(u)))
  u1 <- u[,1]
  u2 <- u[,2]
  eta <- alpha - 1
  ## Joe (1997, p.141)
  ## ((1 + eta * (u1 + u2))^2 - 4 * alpha * eta * u1*u2)^(- 3/2) * alpha * (1 + eta * (u1 + u2 - 2 * u1*u2))
  u1.u2 <- u1 + u2 ; tu12 <- 2*u1*u2
  ## T1 := (1 + eta * u1.u2)^2 - 2 * alpha * eta * tu12 # = 1 + eta*[ u1.u2*(2 + eta*u1.u2) - 2*alpha*tu12 ]
  t1 <- if(log) eta*(u1.u2*(2 + eta*u1.u2) - 2*alpha*tu12) # = T1-1
        else (1 + eta * u1.u2)^2 - 2 * alpha * eta * tu12  # = T1
  t2 <- eta * (u1.u2 - tu12)
  if(log)
    -1.5*log1p(t1) + log(alpha) + log1p(t2)  # = log[ T1^(-3/2) * alpha * (1 + t2) ]
  else
    t1^-1.5 * alpha * (1 + t2)
}

dplackettCopula <- function(u, copula, log = FALSE, ...) {
  ## dim == 2
  if(is.na(alpha <- copula@parameters[1])) stop("parameter is NA")
  .dplackettCopula(u, alpha=alpha, log=log)
}



rplackettCopula <- function(n, copula) {
  u1 <- runif(n)
  u2 <- runif(n)
  psi <- min(1e50, copula@parameters) # "truncate", so even 'Inf' works
  if(is.na(psi)) stop("parameter is NA")
  ## Johnson (1987, p.193)
  a <- u2 * (1 - u2)
  A <- psi + a * (psi - 1)^2
  B <- 2 * a * (u1 * psi^2 + 1 - u1) + psi * (1 - 2 * a)
  D <- sqrt(psi * (psi + 4 * a * u1 * (1 - u1) * (1 - psi)^2))
  v <- (B - (1 - 2 * u2) * D) / 2 / A
  cbind(u1, v, deparse.level=0L)
}



## .plackettTau$assoMeasFun <<- in sysdata.rda <<-- ../inst/docs/tauRho/getSysdataImage.R
## NB: this is not very accurate, e.g., for alpha=1, see ../tests/moments.R
plackettTauFun <- function(alpha) {
  ## ss <- .plackettTau$ss # == 0.5
  ## forwardTransf <- .plackettTau$trFuns$forwardTransf # == x^(ss)
  theta <- sqrt(alpha) # == forwardTransf(alpha, ss)
  valFun <- .plackettTau$assoMeasFun$valFun
  ## ifelse(theta <= 1, valFun(theta), -valFun(1/theta)), more efficiently:
  idx <- theta <= 1
  val <- alpha
  val[ idx] <-   valFun( theta[idx] )
  val[!idx] <- - valFun(1 / theta[!idx])
  val
}
## tauPlackettCopula <- function(copula) plackettTauFun(copula@parameters[1])

plackettdTau <- function(alpha) {
  ss <- .plackettTau$ss # == 0.5
  ## forwardTransf <- .plackettTau$trFuns$forwardTransf # == x^(ss)
  forwardDer <- .plackettTau$trFuns$forwardDer
  valFun <- .plackettTau$assoMeasFun$valFun
  theta <- sqrt(alpha) # == forwardTransf(alpha, ss)
  ## c(ifelse(alpha <= 1, valFun(theta, 1) * forwardDer(alpha, ss),
  ##         valFun(1/theta, 1) * forwardDer(alpha, ss) / theta^2)) :
  idx <- theta <= 1
  val <- alpha
  val[ idx] <- valFun(   theta[idx], 1) * forwardDer(alpha[ idx], ss)
  val[!idx] <- valFun(1/theta[!idx], 1) * forwardDer(alpha[!idx], ss) / theta[!idx]^2
  val
}
## dTauPlackettCopula <- function(copula)   plackettdTau(copula@parameters[[1L]])

iTauPlackettCopula <- function(copula, tau) {
  plackettTauInvLt1 <- approxfun(x = .plackettTau$assoMeasFun$fm$ysmth,
                                 y = .plackettTau$assoMeasFun$fm$x)
  ## theta <- ifelse(tau <= 0, plackettTauInvLt1(tau), 1 / plackettTauInvLt1(-tau)) :
  nT <- tau <= 0
  theta <- tau
  theta[ nT] <-     plackettTauInvLt1( tau[ nT])
  theta[!nT] <- 1 / plackettTauInvLt1(-tau[!nT])
  ## ss <- .plackettTau$ss
  theta^2 ## == .plackettTau$trFuns$backwardTransf(theta, ss)
}

.dRhoPlackett <- function(alpha) {
  ## FIXME: de l'Hopital  for alpha ~= 1
  (2 * (2 - 2 * alpha + (1 + alpha) * log(alpha))) / (alpha - 1)^3
}
dRhoPlackettCopula <- function(copula) .dRhoPlackett(copula@parameters)

.rhoPlackett <- function(alpha) { # (*not* vectorized - FIXME ..)
  if (alpha == 0) -1
  else if(alpha == 1) 0
  ## FIXME: de l'Hopital  for alpha ~= 1
  else (alpha + 1) / (alpha - 1) - 2 * alpha * log(alpha) / (alpha - 1)^2
}
rhoPlackettCopula <- function(copula) .rhoPlackett(copula@parameters)


setMethod("pCopula", signature("matrix", "plackettCopula"), pplackettCopula)
setMethod("dCopula", signature("matrix", "plackettCopula"), dplackettCopula)
## pCopula() and dCopula() *generic* already deal with non-matrix case!

setMethod("rCopula", signature("numeric", "plackettCopula"), rplackettCopula)

setMethod("tau", signature("plackettCopula"),
          function(copula) plackettTauFun(copula@parameters))

setMethod("rho", signature("plackettCopula"), rhoPlackettCopula)

setMethod("iTau", signature("plackettCopula"), iTauPlackettCopula)
setMethod("iRho", signature("plackettCopula"), iRhoCopula)

setMethod("dTau", signature("plackettCopula"),
          function(copula) plackettdTau(copula@parameters))

setMethod("dRho", signature("plackettCopula"), dRhoPlackettCopula)

setMethod("lambda", signature("plackettCopula"), function(copula) c(lower=0, upper=0))
