\name{generator}
\docType{methods}
%----------------------Archimedean Copulas------------------
% iPsi
\alias{iPsi}
\alias{iPsi-methods}
\alias{iPsi,claytonCopula-method}
\alias{iPsi,frankCopula-method}
\alias{iPsi,gumbelCopula-method}
\alias{iPsi,amhCopula-method}
\alias{iPsi,joeCopula-method}
% psi
\alias{genInv}%deprecated
\alias{psi}
\alias{psi-methods}
\alias{psi,claytonCopula-method}
\alias{psi,frankCopula-method}
\alias{psi,gumbelCopula-method}
\alias{psi,amhCopula-method}
\alias{psi,joeCopula-method}
% iPsi
\alias{genFun}%deprecated
\alias{iPsi}
\alias{iPsi-methods}
\alias{iPsi,claytonCopula-method}
\alias{iPsi,frankCopula-method}
\alias{iPsi,gumbelCopula-method}
\alias{iPsi,amhCopula-method}
\alias{iPsi,joeCopula-method}
% diPsi
\alias{genFunDer1}%deprecated
\alias{genFunDer2}%deprecated
\alias{diPsi}
\alias{diPsi-methods}
\alias{diPsi,claytonCopula-method}
\alias{diPsi,frankCopula-method}
\alias{diPsi,gumbelCopula-method}
\alias{diPsi,amhCopula-method}
\alias{diPsi,joeCopula-method}
%----------------------Extreme Value Copulas------------------
% A
\alias{Afun}%deprecated
\alias{A}
\alias{A-methods}
\alias{A,galambosCopula-method}
\alias{A,gumbelCopula-method}
\alias{A,huslerReissCopula-method}
\alias{A,tawnCopula-method}
\alias{A,tevCopula-method}
\alias{A,asymCopula-method}
\alias{A,indepCopula-method}
% dAdu
\alias{AfunDer}%deprecated
\alias{dAdu}
\alias{dAdu-methods}
\alias{dAdu,galambosCopula-method}
\alias{dAdu,gumbelCopula-method}
\alias{dAdu,huslerReissCopula-method}
\alias{dAdu,tawnCopula-method}
\alias{dAdu,tevCopula-method}
% main
\title{Generator Functions for Archimedean and Extreme-Value Copulas}
\description{
  Methods to evaluate the generator function, the inverse generator
  function, and derivatives of the inverse of the generator function for
  Archimedean copulas. For extreme-value copulas, the \dQuote{Pickands
    dependence function} plays the role of a generator function.
}
\usage{
psi(copula, s)
%% TODO?  dPsi(copula, s)
iPsi(copula, u, \dots)
diPsi(copula, u, degree=1, log=FALSE, \dots)

A(copula, w)
dAdu(copula, w)
}
\arguments{
  \item{copula}{an object of class \code{"\linkS4class{copula}"}. }
  \item{u, s, w}{numerical vector at which these functions are to be
    evaluated. }
  \item{\dots}{further arguments for specific families.}
  \item{degree}{the degree of the derivative (defaults to 1).}
  \item{log}{logical indicating if the \code{\link{log}} of the
    \emph{absolute} derivative is desired.  Note that the derivatives
    of \eqn{psi} alternate in sign.}
}
\details{
  \code{psi()} and \code{iPsi()} are, respectively, the generator
  function \eqn{\psi()} and its inverse \eqn{\psi^{(-1)}}{\psi^(-1)} for
  an Archimedean copula, see \code{\link{pnacopula}} for definition and
  more details.

  \code{diPsi()} computes (currently only the first two) derivatives of
  \code{iPsi()} (\eqn{= \psi^{(-1)}}{= \psi^(-1)}).

  \code{A()}, the \dQuote{Pickands dependence function}, can be seen as the
  generator function of an extreme-value copula. For instance, in the
  bivariate case, we have the following result (see, e.g., Gudendorf and
  Segers 2009):

  A bivariate copula \eqn{C} is an extreme-value copula if and only if
  \deqn{C(u,v) = (uv)^{A(\log(v) / \log(uv))}, \qquad (u,v) \in (0,1]^2
  \setminus \{(1,1)\},}{ C(u, v) = (uv)^A(log(v) / log(uv)),   (u,v) in
  (0,1]^2 w/o {(1,1)},} where \eqn{A: [0, 1] \to [1/2, 1]}{A: [0,1] ->
  [1/2, 1]} is convex and satisfies \eqn{\max(t, 1-t) \le A(t) \le
  1}{max(t,1-t) <= A(t) <= 1} for all \eqn{t \in [0, 1]}{t in [0,1]}.
  \cr
  In the \eqn{d}-variate case, there is a similar characterization,
  except that this time, the Pickands dependence function \eqn{A} is
  defined on the \eqn{d}-dimensional unit simplex.

  \code{dAdu()} returns a data.frame containing the 1st and 2nd
  derivative of \code{A()}.
}
\references{
  Gudendorf, G. and Segers, J. (2009)
  {Nonparametric estimation of an extreme-value copula in arbitrary dimensions},
  {ArXiv}, \url{http://arxiv.org/abs/0911.1015}.
}
\seealso{
  Nonparametric estimators for \eqn{A()} are available, see
  \code{\link{An}}.
}
\examples{
## List the available methods (and their definitions):
showMethods("A")
showMethods("iPsi", incl=TRUE)
}
\keyword{multivariate}
