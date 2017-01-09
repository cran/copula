\name{log1mexp}
\alias{log1pexp}
\alias{log1mexp}
\title{Compute  f(a) = \eqn{\mathrm{log}}{log}(1 +/- \eqn{\mathrm{exp}}{exp}(-a))
  Numerically Optimally}
\description{
  Compute f(a) = log(1 - exp(-a)), respectively
  g(x) = log(1 + exp(x))  quickly numerically accurately.
}
\usage{
log1mexp(a, cutoff = log(2))
log1pexp(x, c0 = -37, c1 = 18, c2 = 33.3)
}
\arguments{
  \item{a}{numeric vector of positive values}
  \item{x}{numeric vector}
  \item{cutoff}{positive number; \code{log(2)} is \dQuote{optimal},
    %% see below, TODO
    but the exact value is unimportant, and anything in
    \eqn{[0.5, 1]} is fine.}
  \item{c0, c1, c2}{cutoffs for \code{log1pexp}; see below.}
}
\value{
  f(a) == log(1 - exp(-a)) == log1p(-exp(-a)) ==
  log(-expm1(-a))

  or

  g(x) == log(1 + exp(x)) == log1p(exp(x))

  computed accurately and quickly
}
%\author{Martin Maechler, May 2002; \code{log1pexp()} in 2012}
\references{%  ~/R/Pkgs/Rmpfr/vignettes/log1mexp-note.Rnw
  Martin \enc{Mächler}{Maechler} (2012).
  Accurately Computing \eqn{\log(1-\exp(-|a|))};
  \url{https://CRAN.R-project.org/package=Rmpfr/vignettes/log1mexp-note.pdf}.
  %  see also <<log1mexp>> in ../inst/doc/Frank-Rmpfr.Rnw
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
head(fa)# expm1() works here
tail(fa)# log1p() works here

## graphically:
lwd <- 1.5*(5:2); col <- adjustcolor(1:4, 0.4)
op <- par(mfcol=c(1,2), mgp = c(1.25, .6, 0), mar = .1+c(3,2,1,1))
  matplot(a, fa, type = "l", log = "x", col=col, lwd=lwd)
  legend("topleft", fExpr, col=col, lwd=lwd, lty=1:4, bty="n")
  # expm1() & log1mexp() work here

  matplot(a, -fa, type = "l", log = "xy", col=col, lwd=lwd)
  legend("left", paste("-",fExpr), col=col, lwd=lwd, lty=1:4, bty="n")
  # log1p() & log1mexp() work here
par(op)

curve(log1pexp, -10, 10, asp=1)
abline(0,1, h=0,v=0, lty=3, col="gray")

## Cutoff c1 for log1pexp() -- not often "needed":
curve(log1p(exp(x)) - log1pexp(x), 16, 20, n=2049)
## need for *some* cutoff:
x <- seq(700, 720, by=2)
cbind(x, log1p(exp(x)), log1pexp(x))

## Cutoff c2 for log1pexp():
curve((x+exp(-x)) - x, 20, 40, n=1025)
curve((x+exp(-x)) - x, 33.1, 33.5, n=1025)
}
\keyword{math}
