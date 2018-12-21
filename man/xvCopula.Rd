\name{xvCopula}
\title{Model (copula) selection based on \code{k}-fold cross-validation}
\alias{xvCopula}
\description{
  Computes the leave-one-out cross-validation criterion (or a
  \code{k}-fold version of it) for the hypothesized parametric copula
  family using, by default, maximum pseudo-likelihood estimation.

  The leave-one-out criterion is a crossvalidated log likelihood. It is
  denoted by \eqn{\widehat{xv}_n}{xv^n} in Gr\enc{ø}{o}nneberg and Hjort
  (2014) and defined in equation (42) therein. When computed for several
  parametric copula families, it is thus meaningful to select the family
  maximizing the criterion.

  For \eqn{k < n}, \eqn{n} the sample size, the \code{k}-fold version is
  an approximation of the leave-one-out criterion that uses \eqn{k}
  randomly chosen (almost) equally sized data blocks instead of \eqn{n}.
  When \eqn{n} is large, \eqn{k}-fold cross-validation is considerably
  faster (if \eqn{k} is \dQuote{small} compared to \eqn{n}).
}
\usage{
xvCopula(copula, x, k = NULL, verbose = interactive(),
         ties.method = eval(formals(rank)$ties.method), \dots)
}
\arguments{
  \item{copula}{object of class \code{"\linkS4class{copula}"}
    representing the hypothesized copula family.}
  \item{x}{a data matrix that will be transformed to pseudo-observations.}
  \item{k}{the number of data blocks; if \code{k = NULL}, \code{nrow(x)}
    blocks are considered (which corresponds to leave-one-out
    cross-validation).}
  \item{verbose}{a logical indicating if progress of the cross validation
    should be displayed via \code{\link[utils]{txtProgressBar}}.}
  \item{ties.method}{string specifying how ranks should be computed if
     there are ties in any of the coordinate samples of \code{x} and
     fitting is based on maximum pseudo-likelihood; passed
     to \code{\link{pobs}}.}
  \item{\dots}{additional arguments passed to \code{\link{fitCopula}()}.}
}
%\details{}
\value{
  A real number equal to the cross-validation criterion
  multiplied by the sample size.% \eqn{n}
}
%\author{Ivan Kojadinovic and Martin Maechler}
\references{
  Gr\enc{ø}{o}nneberg, S., and Hjort, N.L. (2014)
  The copula information criteria.
  \emph{Scandinavian Journal of Statistics} \bold{41}, 436--459.
}
\note{Note that \code{k}-fold cross-validation with \code{k < n}
  shuffles the lines of \code{x} prior to forming the blocks.  The
  result thus depends on the value of the random seed.

  The default estimation method is maximum pseudo-likelihood estimation
  but this can be changed if necessary along with all the other
  arguments of \code{\link{fitCopula}()}.
}
\seealso{
  \code{\link{fitCopula}()} for the underlying estimation procedure and
  \code{\link{gofCopula}()} for goodness-of-fit tests.
}
\examples{
\dontshow{set.seed(47)}% reproducibility
## A two-dimensional data example ----------------------------------
x <- rCopula(200, claytonCopula(3))

\donttest{% slow
## Model (copula) selection -- takes time: each fits 200 copulas to 199 obs.
xvCopula(gumbelCopula(), x)
xvCopula(frankCopula(), x)
xvCopula(joeCopula(), x)
xvCopula(claytonCopula(), x)
xvCopula(normalCopula(), x)
xvCopula(tCopula(), x)
xvCopula(plackettCopula(), x)
}% dont test

## The same with 5-fold cross-validation [to save time ...]
set.seed(1) # k-fold is random (for k < n) !
xvCopula(gumbelCopula(),  x, k=5)
xvCopula(frankCopula(),   x, k=5)
xvCopula(joeCopula(),     x, k=5)
xvCopula(claytonCopula(), x, k=5)
xvCopula(normalCopula(),  x, k=5)
xvCopula(tCopula(),       x, k=5)
xvCopula(plackettCopula(),x, k=5)
}
\keyword{models}
\keyword{multivariate}
