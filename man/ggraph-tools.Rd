\name{ggraph-tools}
\alias{pairwiseCcop}
\alias{pairwiseIndepTest}
\alias{pviTest}
\alias{gpviTest}
\title{Computations for Graphical GOF Test via Pairwise Rosenblatt Transforms}
\description{
  Tools for computing a graphical goodness-of-fit (GOF) test based on pairwise
  Rosenblatt transformed data.
  \describe{
    \item{\code{pairwiseCcop()} }{computes a \eqn{(n,d,d)}-\code{\link{array}}
      which contains pairwise Rosenblatt-transformed data.}
    \item{\code{pairwiseIndepTest()} }{takes such an array as input and
      computes a \eqn{(d,d)}-\code{\link{matrix}} of test results from
      pairwise tests of independence (as by \code{\link{indepTest}()}).}
    \item{\code{pviTest()} }{can be used to extract the matrix of p-values from
      the return matrix of \code{pairwiseIndepTest()}.}
    \item{\code{gpviTest()} }{takes such a matrix of p-values and computes a
      global p-value with the method provided.}
  }
}
\usage{
pairwiseCcop(u, copula, ...)
pairwiseIndepTest(cu.u, N=256,
        iTest = indepTestSim(n, p=2, m=2, N=N, verbose = idT.verbose, ...),
        verbose=TRUE, idT.verbose = verbose, ...)

 pviTest(piTest)
gpviTest(pvalues, method=p.adjust.methods, globalFun=min)
}
\arguments{
  \item{u}{\eqn{(n,d)}-\code{\link{matrix}} of copula data.}
  \item{copula}{copula object used for the Rosenblatt transform
    (\eqn{H_0}{H[0]} copula).}
  \item{\dots}{additional arguments passed to the internal function
    which computes the conditional copulas (for \code{pairwiseCcop()}).
    Can be used to pass, for example, the degrees of freedom
    parameter \code{df} for t-copulas.

    For \code{pairwiseIndepTest()},
    \dots are passed to \code{\link{indepTestSim}()}.}

  \item{cu.u}{\eqn{(n,d,d)}-\code{\link{array}} as returned by
    \code{pairwiseCcop()}.}
  \item{N}{argument of \code{\link{indepTestSim}()}.}
  \item{iTest}{the result of (a version of) \code{\link{indepTestSim}()};
    as it does \emph{not} depend on the data, and is costly to compute,
    it can be computed separately and passed here.}
  \item{verbose}{\code{\link{integer}} (or \code{\link{logical}})
    indicating if and how much progress should be printed during the
    computation of the tests for independence.}
  \item{idT.verbose}{logical, passed as \code{verbose} argument to
    \code{\link{indepTestSim}()}.}
  \item{piTest}{\eqn{(d,d)}-\code{\link{matrix}} of \code{indepTest}
    objects as returned by \code{pairwiseIndepTest()}.}

  \item{pvalues}{\eqn{(d,d)}-\code{\link{matrix}} of p-values.}
  \item{method}{\code{\link{character}} vector of adjustment methods for
    p-values; see \code{\link{p.adjust.methods}} for more details.}
  \item{globalFun}{\code{\link{function}} determining how to compute a
    global p-value from a matrix of pairwise adjusted p-values.}
}
\note{
  If \code{u} are distributed according to or \dQuote{perfectly} sampled
  from a copula, p-values on GOF tests for that copula should be uniformly
  distributed in \eqn{[0,1]}.
}
\value{
  \describe{
    \item{pairwiseCcop}{\eqn{(n,d,d)}-\code{\link{array}} \code{cu.u}
      with \code{cu.u[i,j]} containing \eqn{C(u_i\,|\,u_j)}{C(u[,i]|u[,j])}
      for \eqn{i\neq j}{i!=j} and \eqn{u_i}{u[,i]} for \eqn{i=j}.}
    \item{pairwiseIndepTest}{\eqn{(d,d)}-\code{\link{matrix}} of lists
      with test results as returned by \code{\link{indepTest}()}.  The
      test results correspond to pairwise tests of independence as
      conducted by \code{\link{indepTest}()}.}
    \item{pviTest}{\eqn{(d,d)}-\code{\link{matrix}} of p-values.}
    \item{gpviTest}{global p-values for the specified methods.}
  }
}
%\author{Marius Hofert and Martin Maechler}
\references{
  Hofert and \enc{Mächler}{Maechler} (2014),
  see \code{\link{pairsRosenblatt}}.% -> ./pairsRosenblatt.Rd
}
\seealso{
  \code{\link{pairsRosenblatt}} % -> ./pairsRosenblatt.Rd
  for where these tools are used, including
  \code{\link{demo}(gof_graph)} for examples.
}
\examples{
## demo(gof_graph)
% FIXME add some examples using these *tools* ?
}
\keyword{multivariate}
\keyword{htest}
