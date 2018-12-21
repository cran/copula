##' @title Marginal Copula of a Given Copula
##' @param copula input copula
##' @param keep logical vector indicating which margins to keep
##' @return marginal copula (corresponding to the components in 'keep')
##' @note Currently only for ellipCopula and archmCopula
##' @author Jun Yan
##' @note Consider a 5-dim AR1 structured correlation matrix based on rho = 0.8
##'       and use keep = c(FALSE, TRUE, FALSE, TRUE, TRUE) to see that
##'       the new correlation matrix has off-diagonal entries rho^2, rho^3, rho
##'       => neither AR1- nor Toeplitz-structured anymore
setGeneric("margCopula", function(copula, keep) {
    stopifnot(length(keep) == copula@dimension, sum(keep) >= 2,
              is(copula, "tCopula") || is(copula, "normalCopula") || is(copula, "archmCopula"))
    standardGeneric("margCopula")
})

## Normal copula
margNormalCopula <- function(copula, keep) {
    dim <- sum(keep) # dimension of the new, marginal copula
    if (copula@dispstr == "ex")
        normalCopula(getTheta(copula), dim = dim, dispstr = "ex")
    else { # ar1, toep, and un all become un
        normalCopula(P2p(getSigma(copula)[keep, keep]),
                     dim = dim, dispstr = "un")
    }
}

## t copula
margTCopula <- function(copula, keep) {
    dim <- sum(keep) # dimension of the new, marginal copula
    if (copula@dispstr == "ex")
        tCopula(copula@getRho(copula), dim = dim, dispstr = "ex",
                df = getdf(copula), df.fixed = copula@df.fixed)
    else { # ar1, toep, and un all become un
        ## Note: (At least) if the new copula is bivariate or 'keep' starts at 1
        ##       and is 'connected', then the structure is actually 'AR1'
        ##       This could be checked via...
        ##       d <- dim(copula) # old copula dimension
        ##       is.AR1 <- (dim == 2) || (keep[1] && all(diff(seq_len(d)[keep]) == 1)) # still 'AR1' structure?
        ##       param <- if(is.AR1) copula@getRho(copula) else P2p(getSigma(copula)[keep, keep])
        ##       if(!is.AR1) dispstr <- "toep"
        ##       ... but makes things just more complicated
        ##       (not always returning an object of the same type)
        tCopula(P2p(getSigma(copula)[keep, keep]),
                dim = dim, dispstr = "un",
                df = getdf(copula), df.fixed = copula@df.fixed)
    }
}

setMethod("margCopula", signature("normalCopula", "logical"), margNormalCopula)
setMethod("margCopula", signature("tCopula",      "logical"), margTCopula)

## Archimedean copulas
setMethod("margCopula", signature("archmCopula", "logical"),
          function(copula, keep) {
              dim <- sum(keep)
              copula@dimension <- dim
              copula
          })
