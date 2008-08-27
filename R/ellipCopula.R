getSigma <- function(copula) {
  dim <- copula@dimension
  rho <- copula@getRho(copula)
  sigma <- diag(dim)
  if (copula@dispstr == "ex") {
    sigma[lower.tri(sigma)] <- rho[1]
    sigma[upper.tri(sigma)] <- rho[1]
  }
  else if (copula@dispstr == "ar1") {
    for (i in 1:dim)  for (j in 1:dim)  sigma[i,j] <- rho ^ abs(i - j)
  }
  else if (copula@dispstr == "un") {
    sigma[lower.tri(sigma)] <- rho
    sigma[upper.tri(sigma)] <- t(sigma)[upper.tri(sigma)]
  }
  else if (copula@dispstr == "toep") {
    for (i in 1:dim) for (j in 1:dim)
      if (i != j) sigma[i,j] <- rho[abs(i - j)]
  }
  sigma
}

ellipCopula <- function(family, param, dim = 2, dispstr = "ex", df = 5, ...) {
  familiesImplemented <- c("normal", "t")
  fam <- pmatch(family, familiesImplemented, -1)
  if (fam == -1)
    stop(paste("Valid family names are", familiesImplemented))
  copula <- switch(fam,
                   normalCopula(param, dim = dim, dispstr = dispstr),
                   tCopula(param, dim = dim, dispstr = dispstr, df = df)
                   )
  copula
}

calibKendallsTauEllipCopula <- function(copula, tau) {
  sin((tau * pi) / 2)
}

calibSpearmansRhoEllipCopula <- function(copula, rho) {
  sin (pi * rho / 6) * 2
}

tauDerEllipCopula <- function(copula)
  {
    return( 2 / (pi * sqrt(1 - copula@parameters^2)) )
  }

rhoDerEllipCopula <- function(copula)
  {
    return( 6 / (pi * sqrt(4 - copula@parameters^2)) )
  }

setMethod("calibKendallsTau", signature("ellipCopula"), calibKendallsTauEllipCopula)
setMethod("calibSpearmansRho", signature("ellipCopula"), calibSpearmansRhoEllipCopula)

setMethod("tauDer", signature("ellipCopula"), tauDerEllipCopula)
setMethod("rhoDer", signature("ellipCopula"), rhoDerEllipCopula)
