archmCopula <- function(family, param, dim = 2, ...) {
  familiesImplemented <- c("clayton", "frank", "gumbel", "amh")
  fam <- pmatch(family, familiesImplemented, -1)
  if (fam == -1)
    stop(paste("Valid family names are", familiesImplemented))
  copula <- switch(fam,
                   claytonCopula(param, dim = dim),
                   frankCopula(param, dim = dim),
                   gumbelCopula(param, dim = dim),
                   amhCopula(param, dim = 2)
                   )
  copula
}



kendallsTauArchmCopula <- function(copula) {
  integrand <- function(x) genFun(copula, x) / genFunDer1(copula, x)
  1 + 4 * integrate(integrand, 0, 1)$value
}

setMethod("kendallsTau", signature("archmCopula"), kendallsTauArchmCopula)
