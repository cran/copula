\name{coeffG}
\title{Coefficients of Polynomial used for Gumbel Copula}
\alias{coeffG}
\description{
  Compute the coefficients \eqn{a_{d,k}(\theta)}{a[d,k](\theta)} involved in the
  generator (psi) derivatives and the copula density of Gumbel copulas.

  For non-small dimensions \eqn{d}, these are numerically challenging to
  compute accurately.
}
\usage{
coeffG(d, alpha,
       method = c("sort", "horner", "direct", "dsumSibuya",
                  paste("dsSib", eval(formals(dsumSibuya)$method), sep = ".")),
       log = FALSE, verbose = FALSE)
}
\arguments{
  \item{d}{number of coefficients, (the copula dimension), d >= 1.}
  \item{alpha}{parameter \eqn{1/\theta} in \eqn{(0,1]}; you may use
    \code{\link[Rmpfr]{mpfr}(alph, precBits = <n_prec>)}
    for higher precision methods (\code{"Rmpfr*"}) from package
    \CRANpkg{Rmpfr}.}
  \item{method}{a \code{\link{character}} string, one of
    \describe{
      \item{\code{"sort"}:}{compute coefficients via \eqn{exp(log())}
	pulling out the maximum, and sort.}
      \item{\code{"horner"}:}{uses polynomial evaluation, our internal
	\code{polynEval()}.}
      \item{\code{"direct"}:}{brute force approach.}
      \item{\code{"dsSib.<FOO>"}:}{uses \code{\link{dsumSibuya}(..., method= "<FOO>")}.}
    }
  }
  \item{log}{logical determining if the logarithm (\code{\link{log}}) is
    to be returned.}
  \item{verbose}{logical indicating if some information should be shown,
    currently for \code{method == "sort"} only.}
}
\value{
  a numeric vector of length \code{d}, of values
  \deqn{% latex
    a_k(\theta, d) = (-1)^{d-k}\sum_{j=k}^d \alpha^j * s(d,j) * S(j,k),
    k \in \{1,\ldots,d\}.}{% non-latex; can use greek letters:
    a_k(\theta, d) = (-1)^(d-k) Sum(j=k..d; \alpha^j * s(d,j) * S(j,k)),
    k in 1..d.
  }
}
\author{
  Marius Hofert and Martin Maechler
}
\note{There are still known numerical problems (with non-"Rmpfr" methods; and
  those are slow), e.g., for d=100,
  alpha=0.8 and  \eqn{sign(s(n,k)) = (-1)^{n-k}}{sign(s(n,k)) = (-1)^(n-k)}.

  As a consequence, the \code{method}s and its defaults may change in
  the future, and so the exact implementation of \code{coeffG()} is
  still considered somewhat experimental.
}
\examples{
a.k  <- coeffG(16, 0.55)
plot(a.k, xlab = quote(k), ylab = quote(a[k]),
     main = "coeffG(16, 0.55)", log = "y", type = "o", col = 2)
a.kH <- coeffG(16, 0.55, method = "horner")
stopifnot(all.equal(a.k, a.kH, tol = 1e-11))# 1.10e-13 (64-bit Lnx, nb-mm4)
%% maybe more from ../demo/G_ak.R
}
\keyword{arith}
