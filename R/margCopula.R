##' @title Marginal copula of a given copula
##' @param copula an input copula
##' @param keep a logical vector indicating which margins to keep
##' @return marginal copula
##' @note Currently only for ellipCopula and archmCopula
##' @author Jun Yan

setGeneric("margCopula", function(copula, keep) {
    stopifnot(length(keep) == copula@dimension, sum(keep) >= 2,
              is(copula, "tCopula") || is(copula, "normalCopula") || is(copula, "archmCopula"))
    standardGeneric("margCopula")
})

margNormalCopula <- function(copula, keep) {
    dim <- sum(keep)
    if (copula@dispstr == "ex")
        normalCopula(getTheta(copula), dim = dim, dispstr = "ex")
    else { # ar1, toep, and un all become un
        sigma <- getSigma(copula)[keep, keep]
        param <- sigma[lower.tri(sigma)]
        normalCopula(param, dim = dim, dispstr = "un")
    }
}

margTCopula <- function(copula, keep) {
    dim <- sum(keep)
    if (copula@dispstr == "ex")
        tCopula(copula@getRho(copula), dim = dim, dispstr = "ex",
                df = getdf(copula), df.fixed = copula@df.fixed)
    else { # ar1, toep, and un all become un
        sigma <- getSigma(copula)[keep, keep]
        param <- sigma[lower.tri(sigma)]
        tCopula(param, dim = dim, dispstr = "un",
                df = getdf(copula), df.fixed = copula@df.fixed)
    }
}

setMethod("margCopula", signature("normalCopula", "logical"), margNormalCopula)
setMethod("margCopula", signature("tCopula", "logical"), margTCopula)

setMethod("margCopula", signature("archmCopula", "logical"),
          function(copula, keep) {
              dim <- sum(keep)
              copula@dimension <- dim
              copula
          })
