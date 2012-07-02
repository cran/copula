\name{log1mexp}
\alias{log1mexp}
\title{Compute  f(a) = \eqn{\mathrm{log}}{log}(1 - \eqn{\mathrm{exp}}{exp}(-a))
  numerically optimally}
\description{
  Compute f(a) = log(1 - exp(-a)) stably
}
\usage{
  log1mexp(a, cutoff = log(2))
}
\arguments{
  \item{a}{numeric vector of positive values}
  \item{cutoff}{positive number; \code{log(2)} is \dQuote{optimal},
    %% see below, TODO
    but the exact value is unimportant, and anything in
    \eqn{[0.5, 1]} is fine.}
}
\value{
  f(a) == log(1 - exp(-a)) == log1p(-exp(-a)) ==
  log(-expm1(-a))
}
\author{
  Martin Maechler, May 2002.
}
\references{
  Martin Maechler (2012)
  \emph{}
  \url{}
%  see also <<log1mexp>> in ../inst/doc/Frank-Rmpfr.Rnw
%   which has {log(-expm1(-theta))} ... underflows for large
%   theta, and hence should be replaced by the generally
%   accurate %% as shown by Maechler (2011) --- FIXME:
%   finally right that small paper!! %%
%   ========================= ~/R/MM/NUMERICS/log1-exp.R
}
\examples{
a <- 2^seq(-58,10, length = 256)
fExpr <- expression(
          log(1 - exp(-a)),
          log(-expm1(-a)),
          log1p(-exp(-a)),
          log1mexp(a))
names(fExpr) <- c("DEF", "expm1", "log1p", "F")
str(fa <- do.call(cbind, as.list(fExpr)))
matplot(a, fa, type = "l", log = "x")# don't see the big difference
legend("topleft", fExpr, col=1:4, lty=1:4, bty="n")

head(fa)# expm1() works here
tail(fa)# log1p() works here

matplot(a, -fa, type = "l", log = "xy")# don't see the big difference
legend("left", fExpr, col=1:4, lty=1:4, bty="n")
}
\keyword{math}
