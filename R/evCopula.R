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

setMethod("tailIndex", signature("evCopula"), tailIndexEvCopula)
