\name{corKendall}
\alias{corKendall}
\title{(Fast) Computation of Pairwise Kendall's Taus}
\description{
  For a data matrix \code{x}, compute the Kendall's tau
  \dQuote{correlation} matrix, i.e., all pairwise Kendall's taus
  between the columns of \code{x}.

  By default and when \code{x} has no missing values
  (\code{\link{NA}}s), the fast \eqn{O(n log(n))} algorithm of
  \code{\link[pcaPP]{cor.fk}()} is used.
}
\usage{
corKendall(x, checkNA = TRUE,
           use = if(checkNA && anyNA(x)) "pairwise" else "everything")
}
\arguments{
  \item{x}{data, a n x p matrix (or less efficiently a
    data.frame), or a numeric vector which is treated as n x 1 matrix.}
  \item{checkNA}{logical indicating if \code{x} should be checked for
    \code{\link{NA}}s and in the case of NA's \emph{and} when \code{use}
    is not specified (\code{\link{missing}}), \code{cor(*, use =
      "pairwise")} should be used.  Note that \code{corKendall(x,
      checkNA = FALSE)} will produce an error when \code{x} has NA's.}
  \item{use}{a string to determine the treatment of \code{\link{NA}}s in
    \code{x}, see \code{\link{cor}}; its default determined via
    \code{checkNA}.  When this differs from \code{"everything"}, \R's
    \code{\link{cor}} is used; otherwise \CRANpkg{pcaPP}'s
    \code{\link[pcaPP]{cor.fk}()} which cannot deal with \code{\link{NA}}s.}
}
% \author{Martin Maechler}
\value{
  The \eqn{p \times p}{p x p} matrix \eqn{K} of pairwise Kendall's taus, with
  \code{K[i,j] := tau(x[,i], x[,j])}.
}
\seealso{
  \code{\link[pcaPP]{cor.fk}()} from \CRANpkg{pcaPP} (used by default
  when there are no missing values (\code{NA}s) in \code{x}).

  \code{\link{etau}()} or \code{\link{fitCopula}(*, method = "itau")}
  make use of \code{corKendall()}.
}
\examples{
## If there are no NA's, corKendall() is faster than cor(*, "kendall")
## and gives the same :

system.time(C1 <- cor(swiss, method="kendall"))
system.time(C2 <- corKendall(swiss))
stopifnot(all.equal(C1, C2,  tol = 1e-5))

## In the case of missing values (NA), corKendall() reverts to
## cor(*, "kendall", use = "pairwise") {no longer very fast} :

swM <- swiss # shorter names and three missings:
colnames(swM) <- abbreviate(colnames(swiss), min=6)
swM[1,2] <- swM[7,3] <- swM[25,5] <- NA
(C3 <- corKendall(swM)) # now automatically uses the same as
stopifnot(identical(C3, cor(swM, method="kendall", use="pairwise")))
## and is quite close to the non-missing "truth":
stopifnot(all.equal(unname(C3), unname(C2), tol = 0.06)) # rel.diff.= 0.055

try(corKendall(swM, checkNA=FALSE)) # --> Error
## the error is really from  pcaPP::cor.fk(swM)
}
