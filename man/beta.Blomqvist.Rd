\name{beta.Blomqvist}
\title{Sample and Population Version of Blomqvist's Beta for Archimedean Copulas}
\alias{beta.}
\alias{betan}
\alias{beta.hat}%deprecated
\description{
  Compute the population (\code{beta.()}) and sample
  (\code{betan()}) version of Blomqvist's beta for an Archimedean
  copula.

  See the reference below for definitions and formulas.
}
\usage{
beta.(cop, theta, d, scaling=FALSE)
betan(u, scaling=FALSE)
}
\arguments{
  \item{cop}{an Archimedean copula (of dimension \eqn{d}) to be
    estimated.}
  \item{theta}{copula parameter.}
  \item{d}{dimension.}
  \item{scaling}{logical, if true, the factors 2^(d-1)/(2^(d-1)-1) and
    2^(1-d) in Blomqvist's beta are omitted.}
  \item{u}{For \code{betan}: 	(\eqn{n\times d}{n x d})-matrix of
      d-dimensional observations distributed according to the copula.}
}
\value{
  \describe{
    \item{\code{beta.}:}{a number, being the population version of
      Blomqvist's beta for the corresponding
      \ifelse{latex}{Archi-medean}{Archimedean} copula;}
    \item{\code{betan}:}{a number, being the sample version of
      Blomqvist's beta for the given data.}
  }
}
\references{
  Schmid and Schmidt (2007),
  Nonparametric inference on multivariate versions of Blomqvist's beta
  and related measures of tail dependence,
  \emph{Metrika} \bold{66}, 323--354.
}
\author{Marius Hofert}
\seealso{
  \code{\linkS4class{acopula}}
}
\examples{
beta.(copGumbel, 2.5, d = 5)

d.set <- c(2:6, 8, 10, 15, 20, 30)
cols <- adjustcolor(colorRampPalette(c("red", "orange", "blue"),
                                     space = "Lab")(length(d.set)), 0.8)
## AMH:
for(i in seq_along(d.set))
   curve(Vectorize(beta.,"theta")(copAMH, x, d = d.set[i]), 0, .999999,
         main = "Blomqvist's beta(.) for  AMH",
         xlab = expression(theta), ylab = expression(beta(theta, AMH)),
         add=(i > 1), lwd=2, col=cols[i])
mtext("NB:  d=2 and d=3 are the same")
legend("topleft", paste("d =",d.set), bty="n", lwd=2, col=cols)

## Gumbel:
for(i in seq_along(d.set))
   curve(Vectorize(beta.,"theta")(copGumbel, x, d = d.set[i]), 1, 10,
         main = "Blomqvist's beta(.) for  Gumbel",
         xlab = expression(theta), ylab = expression(beta(theta, Gumbel)),
         add=(i > 1), lwd=2, col=cols[i])
legend("bottomright", paste("d =",d.set), bty="n", lwd=2, col=cols)

## Clayton:
for(i in seq_along(d.set))
   curve(Vectorize(beta.,"theta")(copClayton, x, d = d.set[i]), 1e-5, 10,
         main = "Blomqvist's beta(.) for  Clayton",
         xlab = expression(theta), ylab = expression(beta(theta, Gumbel)),
         add=(i > 1), lwd=2, col=cols[i])
legend("bottomright", paste("d =",d.set), bty="n", lwd=2, col=cols)

## Joe:
for(i in seq_along(d.set))
   curve(Vectorize(beta.,"theta")(copJoe, x, d = d.set[i]), 1, 10,
         main = "Blomqvist's beta(.) for  Joe",
         xlab = expression(theta), ylab = expression(beta(theta, Gumbel)),
         add=(i > 1), lwd=2, col=cols[i])
legend("bottomright", paste("d =",d.set), bty="n", lwd=2, col=cols)

## Frank:
for(i in seq_along(d.set))
   curve(Vectorize(beta.,"theta")(copFrank, x, d = d.set[i]), 1e-5, 50,
         main = "Blomqvist's beta(.) for  Frank",
         xlab = expression(theta), ylab = expression(beta(theta, Gumbel)),
         add=(i > 1), lwd=2, col=cols[i])
legend("bottomright", paste("d =",d.set), bty="n", lwd=2, col=cols)

## Shows the numeric problems:
curve(Vectorize(beta.,"theta")(copFrank, x, d = 29), 35, 42, col="violet")
}
\keyword{multivariate}
\keyword{distribution}
