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


##' @title A (Pickands dependence) function of a bivariate tevCopula
##' @param copula a bivariate tevCopula object
##' @param w points at which A is to be evaluated
##' @return A(w) for the tevCopula
##' @author Jun Yan
ATev <- function(copula, w) {
  rho <- copula@parameters[1]
  nu <- getdf(copula) # defined in tCopula.R
  wnu <- (w / (1 - w))^(1 / nu)
  x <- (wnu - rho) / sqrt(1 - rho^2) * sqrt(nu + 1)
  y <- (1 / wnu - rho) / sqrt(1 - rho^2) * sqrt(nu + 1)
  A <- w * pt(x, nu + 1) + (1 - w) * pt(y, nu + 1)
  ifelse(w == 0 | w == 1, 1, A)
}

##' @title First and second derivatives of A of a bivariate tevCopula
##' @param copula a bivariate tevCopula object
##' @param w points at which A is to be evaluated
##' @return data frame with two coplumns: A'(w) and A''(w)
dAduTev <- function(copula, w) {
  rho <- copula@parameters[1]
  nu <- getdf(copula) # defined in tCopula.R
  ## prepare dx, dy
  wnu <- (w / (1 - w))^(1 / nu)
  x <- (wnu - rho) / sqrt(1 - rho^2) * sqrt(nu + 1)
  y <- (1 / wnu - rho) / sqrt(1 - rho^2) * sqrt(nu + 1)
  dx <- eval(deriv(~ (w / (1 - w))^(1 / nu) / sqrt(1 - rho^2) * sqrt(nu + 1), "w", hessian=TRUE), list(w = w))
  dxdw <- c(attr(dx, "gradient"))
  d2xdw2 <- c(attr(dx, "hessian"))
  dy <- eval(deriv(~ (w / (1 - w))^( - 1 / nu) / sqrt(1 - rho^2) * sqrt(nu + 1), "w", hessian=TRUE), list(w = w))
  dydw <- c(attr(dy, "gradient"))
  d2ydw2 <- c(attr(dy, "hessian"))
  ## prepare ddtx, ddty, derivative of the t-(nu) density
  ## a <-  gamma(0.5 * (nu + 1)) / gamma(0.5 * nu) / sqrt(nu * pi)
  dens <- ~ gamma(0.5 * (nu + 1)) / gamma(0.5 * nu) / sqrt(nu * pi)  * (1 + u^2 / nu)^(-0.5 * (nu + 1))
  ddens <- deriv(dens, "u")
  ddtx <- c(attr(eval(ddens, list(u = x, nu = nu + 1)), "gradient"))
  ddty <- c(attr(eval(ddens, list(u = y, nu = nu + 1)), "gradient"))
  ## now collect the results
  der1 <- pt(x, nu + 1) + w * dt(x, nu + 1) * dxdw - pt(y, nu + 1) + (1 - w) * dt(y, nu + 1) * dydw
  der2 <- dt(x, nu + 1) * dxdw +
    dt(x, nu + 1) * dxdw + w * ddtx * dxdw^2 + w * dt(x, nu + 1) * d2xdw2 +
      (- dt(y, nu + 1) * dydw ) +
        (- dt(y, nu + 1) * dydw) + (1 - w) * ddty * dydw^2 + (1 - w) * dt(y, nu + 1) * d2ydw2
  data.frame(der1 = der1, der2 = der2)
}

##' @title Constructor of a bivariate tevCopula
##' @param param dispersion parameter of the tevCopula
##' @param df numeric scalar, degrees of freedom
##' @param df.fixed logical, whether df is fixed
##' @return a bivariate tevCopula
##' FIXME: we should allow multivariate tevCopula???
tevCopula <- function(param = NA_real_, df = 4, df.fixed = FALSE) {
  dim <- 2L
  pdim <- length(param)
  parameters <- c(param, df)
  param.names <- c(paste("rho", 1:pdim, sep="."), "df")
  param.lowbnd <- c(rep(-1, pdim), 1e-6)
  param.upbnd <- c(rep(1, pdim), Inf)
  ## similar to tCopula
  ## attr(parameters, "fixed") <-
  ##     c(if(is.null(fixed <- attr(param, "fixed"))) ## IK: parameters -> param
  ##           rep(FALSE, pdim)
  ##       else fixed,
  ##       df.fixed)

  ## IK: new version
  attr(parameters, "fixed") <- c(!isFreeP(param), df.fixed)

  new("tevCopula",
             dimension = dim,
             parameters = parameters,
             df.fixed = df.fixed,
             param.names = param.names,
             param.lowbnd = param.lowbnd,
             param.upbnd = param.upbnd,
             fullname = "<deprecated slot>")# paste("t-EV copula family", if(df.fixed) paste("df fixed at", df))
}

##' @title Evaluation of distribution function of a tevCopula at u
ptevCopula <- function(u, copula) {
  ## dim <- copula@dimension
  ## for (i in 1:dim) assign(paste0("u", i), u[,i])
  ## dim = 2
  u1 <- u[,1]; u2 <- u[,2]
  p <- (r <- uu <- u1 * u2) > 0
  p <- p & (nna <- !is.na(p)) # p: positive uu
  logu <- log(uu[p])
  r[p ] <- exp(logu * ATev(copula, log(u2[p]) / logu))
  r[!p & nna] <- 0 # when one u is zero
  r
}

##' @title Evaluation of density function of a tevCopula at u
dtevCopula <- function(u, copula, log=FALSE, ...) {
  ## dim = 2
  u1 <- u[,1]; u2 <- u[,2]
  C <- ptevCopula(u, copula)
  logu <- log(u1 * u2)
  w <- log(u2) / logu
  wexpr <- ~ log(u2) / log(u1 * u2)
  dw <- eval(deriv(wexpr, c("u1", "u2"), hessian=TRUE), list(u1 = u1, u2 = u2))
  dwdu1 <- c(attr(dw, "gradient")[,"u1"])
  dwdu2 <- c(attr(dw, "gradient")[,"u2"])
  d2wdu1du2 <- c(attr(dw, "hessian")[,"u1","u2"])
  A <- ATev(copula, w)
  Ader <- dAduTev(copula, w)
  Ader1 <- Ader$der1; Ader2 <- Ader$der2
  ## dCdu1 <- C * (1 / u1 * A + logu * Ader1 * dwdu1)
  dCdu2 <- C * (1 / u2 * A + logu * Ader1 * dwdu2)
  ## pdf = d2Cdu1du2
  pdf <- dCdu2 * (1 / u1 * A + logu * Ader1 * dwdu1) +
    C * (1 / u1 * Ader1 * dwdu2 + 1 / u2 * Ader1 * dwdu1 +
         logu * Ader2 * dwdu2 * dwdu1 + logu * Ader1 * d2wdu1du2)

  ## FIXME: improve log-case
  if(log) log(pdf) else pdf
}

##' @title 1st and 2nd derivatives of C wrt to arguments of a bivariate symmetric evCopula
##' TODO:
dCduSymEvCopula <- function(copula, u, ...) {
  mat <- matrix(NA, nrow(u), 2)
  pcop <- pCopula(u, copula)
  loguv <- log(u[,1]) + log(u[,2])
  w <- log(u[,2]) / loguv
  a <- A(copula, w)
  aDer <- dAdu(copula, w)$der1
  mat[,1] <- pcop * (a / u[,1] - loguv * aDer * log(u[,2]) / (loguv)^2 / u[,1])
  mat[,2] <- pcop * (a / u[,2] + loguv * aDer * log(u[,1]) / (loguv)^2 / u[,2])
  mat
}

setMethod("dCdu", signature("tevCopula"), dCduSymEvCopula)


## This block is copied from ../../copulaUtils/assoc/ ##########################

## tau

tevTauFun <- function(alpha) {
  ss <- .tevTau$ss
  forwardTransf <- .tevTau$trFuns$forwardTransf
  valFun <- .tevTau$assoMeasFun$valFun
  theta <- forwardTransf(alpha, ss)

  valFun(theta)
}

tauTevCopula <- function(copula) {
  alpha <- copula@parameters[1]
  tevTauFun(alpha)
}

iTauTevCopula <- function(copula, tau) {
    if (any(neg <- tau < 0)) {
        warning("For the t-ev copula, tau must be >= 0. Replacing negative values by 0.")
        tau[neg] <- 0
    }
    tevTauInv <- approxfun(x = .tevTau$assoMeasFun$fm$ysmth,
                           y = .tevTau$assoMeasFun$fm$x, rule=2)

    ss <- .tevTau$ss
    theta <- tevTauInv(tau)
    .tevTau$trFuns$backwardTransf(theta, ss)
}

tevdTau <- function(alpha) {
  ss <- .tevTau$ss
  forwardTransf <- .tevTau$trFuns$forwardTransf
  forwardDer <- .tevTau$trFuns$forwardDer
  valFun <- .tevTau$assoMeasFun$valFun
  theta <- forwardTransf(alpha, ss)

  valFun(theta, 1) * forwardDer(alpha, ss)
}

dTauTevCopula <- function(copula) {
  alpha <- copula@parameters[1]
  tevdTau(alpha)
}

## rho

tevRhoFun <- function(alpha) {
  ss <- .tevRho$ss
  forwardTransf <- .tevRho$trFuns$forwardTransf
  valFun <- .tevRho$assoMeasFun$valFun
  theta <- forwardTransf(alpha, ss)

  valFun(theta)
}

rhoTevCopula <- function(copula) {
  alpha <- copula@parameters[1]
  tevRhoFun(alpha)
}

iRhoTevCopula <- function(copula, rho) {
    if (any(neg <- rho < 0)) {
        warning("For the t-ev copula, rho must be >= 0. Replacing negative values by 0.")
        rho[neg] <- 0
     }

    tevRhoInv <- approxfun(x = .tevRho$assoMeasFun$fm$ysmth,
                           y = .tevRho$assoMeasFun$fm$x, rule = 2)

    ss <- .tevRho$ss
    theta <- tevRhoInv(rho)
    .tevRho$trFuns$backwardTransf(theta, ss)
}

tevdRho <- function(alpha) {
  ss <- .tevRho$ss
  forwardTransf <- .tevRho$trFuns$forwardTransf
  forwardDer <- .tevRho$trFuns$forwardDer
  valFun <- .tevRho$assoMeasFun$valFun
  theta <- forwardTransf(alpha, ss)

  valFun(theta, 1) * forwardDer(alpha, ss)
}

dRhoTevCopula <- function(copula) {
  alpha <- copula@parameters[1]
  tevdRho(alpha)
}

################################################################################

setMethod("pCopula", signature("matrix", "tevCopula"), ptevCopula)
setMethod("dCopula", signature("matrix", "tevCopula"), dtevCopula)


setMethod("A", signature("tevCopula"), ATev)
setMethod("dAdu", signature("tevCopula"), dAduTev)

setMethod("tau", signature("tevCopula"), tauTevCopula)
setMethod("rho", signature("tevCopula"), rhoTevCopula)

setMethod("iTau", signature("tevCopula"), iTauTevCopula)
setMethod("iRho", signature("tevCopula"), iRhoTevCopula)

setMethod("dTau", signature("tevCopula"), dTauTevCopula)
setMethod("dRho", signature("tevCopula"), dRhoTevCopula)
