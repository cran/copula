evCopula <- function(family, param, dim = 2, ...) {
  familiesImplemented <- c("galambos", "gumbel", "huslerReiss")
  fam <- pmatch(family, familiesImplemented, -1)
  if (fam == -1)
    stop(paste("Valid family names are", familiesImplemented))
  copula <- switch(fam,
                   galambosCopula(param),
                   gumbelCopula(param),
                   huslerReissCopula(param)
                   )
  copula
}


tailIndexEvCopula <- function(copula) {
  lower <- 0
  upper <- 2 - 2 * Afun(copula, 0.5)
  c(lower=lower, upper=upper)
}


hdensity <- function(copula, z) {
  A <- Afun(copula, z)
  ders <- AfunDer(copula, z)
  1 + (1 - 2 * z) * ders$der1 / A + z * (1 - z) * (ders$der2 * A - ders$der1^2) / A^2
}

revCopula <- function(copula, n) {
  ## Caperaa, Fougeres, and Genest (2000, Journal of Multivariate Analysis)
  ## This algorithm has low efficiency for high dependence, in which case,
  ## hdensity is pretty much centered around 0.5.
  ## In particular, it generates peculiar numbers for galambosCopula with
  ## alpha >= 30. Don't how to solve yet.
  
  M <- hdensity(copula, 0.5) ## maximum obtained at 0.5 for symmetric copula
  z <- rep(NA, n)
  ndone <- 0
  while (TRUE) {
    ucand <- runif(n)
    accept <- runif(n) <= hdensity(copula, ucand) / M
    accept <- accept & (!is.na(accept)) ## hdensity can be NA at some ucand
    ngood <- sum(accept)
    if (ngood == 0) next
    ngood <- min(ngood, n - ndone)
    z[ndone + 1:ngood] <- ucand[accept][1:ngood]
    ndone <- ndone + ngood
    if (ndone == n) break    
  }
  ders <- AfunDer(copula, z)
  pz <- z * (1 - z) * ders$der2 / hdensity(copula, z) / Afun(copula, z)
  w <- rep(NA, n)
  mix1 <- runif(n) <= pz
  nmix1 <- sum(mix1)
  if (any(mix1)) w[mix1] <- runif(nmix1)
  if (any(!mix1)) w[!mix1] <- runif(n - nmix1) * runif(n - nmix1)
  ## CFG (2000, p.39)
  cbind(exp(z * log(w)/Afun(copula, z)),
        exp((1 - z) * log(w)/Afun(copula, z)))
}


#### These one-dimensional numerical integration is quite accurate.
#### They are much better than two-dimensional integration function adapt.

kendallsTauEvCopula <- function(copula) {
  integrand <- function(x) x * (1 - x) / Afun(copula, x) * AfunDer(copula, x)$der2
  integrate(integrand, 0, 1)$value
}

spearmansRhoEvCopula <- function(copula) {
  integrand <- function(x) 1 / (Afun(copula, x) + 1)^2
  12 * integrate(integrand, 0, 1)$value - 3
}

setMethod("tailIndex", signature("evCopula"), tailIndexEvCopula)
setMethod("rcopula", signature("evCopula"), revCopula)
setMethod("kendallsTau", signature("evCopula"), kendallsTauEvCopula)
setMethod("spearmansRho", signature("evCopula"), spearmansRhoEvCopula)
