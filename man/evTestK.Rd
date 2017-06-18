\name{evTestK}
\alias{evTestK}
\title{Bivariate Test of Extreme-Value Dependence Based on Kendall's Distribution}
\description{
  Test of extreme-value dependence based on the  bivariate probability integral
  transformation. The test statistic is defined in Ben Ghorbal,
  G. \enc{Nešlehová}{Neslehova}, and Genest (2009).
}
\usage{
evTestK(x, method = c("fsample","asymptotic","jackknife"), ties = NA, N = 100)
}
\arguments{
  \item{x}{ a data matrix. }
  \item{method}{ specifies the variance estimation method; can be either
    \code{"fsample"} (finite-sample, the default), \code{"asymptotic"}
    or \code{"jackknife"}.}
  \item{ties}{ logical; if \code{TRUE}, the original test is adapted to
    take the presence of ties in the coordinate samples of \code{x} into
    account; the default value of \code{NA} indicates that the
    presence/absence of ties will be checked for automatically.}
  \item{N}{ number of samples to be used to estimate a bias term if
    \code{ties = TRUE}.}
}
\details{
  The code for this test was generously provided by Johanna
  G. \enc{Nešlehová}{Neslehova}.  More details are available in Appendix
  B of Ben Ghorbal, G. \enc{Nešlehová}{Neslehova} and Genest (2009).
}
\value{
  An object of \code{\link{class}} \code{htest} which is a list,
  some of the components of which are
  \item{statistic}{ value of the test statistic. }
  \item{p.value}{ corresponding p-value. }
}
%\author{Ivan Kojadinovic and Jun Yan}
\references{
  Ghorbal, M. B., Genest, C., and G. \enc{Nešlehová}{Neslehova}, J. (2009)
  On the test of Ghoudi, Khoudraji, and Rivest for extreme-value dependence.
  \emph{The Canadian Journal of Statistics} \bold{37}, 1--9.

  Kojadinovic, I. (2017). Some copula inference procedures adapted to
  the presence of ties. \emph{Computational Statistics and Data
  Analysis} \bold{112}, 24--41, \url{http://arxiv.org/abs/1609.05519}.
}
\seealso{ \code{\link{evTestC}}, \code{\link{evTestA}}, \code{\link{evCopula}},
  \code{\link{gofEVCopula}}, \code{\link{An}}.
}
\examples{
set.seed(321)
## Do the data come from an extreme-value copula?
evTestK(Ug <- rCopula(200, gumbelCopula(3))) # not significant => yes, EV
    dim(Uc <- rCopula(200, claytonCopula(3)))
## Clayton:                       tests are highly significant => no, not EV
(K1 <- evTestK(Uc))
(K2 <- evTestK(Uc, method = "asymptotic"))
\donttest{
system.time(print(K3 <- evTestK(Uc, method = "jackknife")))
## slower: 1.06 sec (2015 intel core i7)
}
\dontshow{
Xtra <- copula:::doExtras()
stopifnot(all.equal(K1$p.value, 3.9471315e-9, tolerance = .001),
          all.equal(K2$p.value, 2.1277296e-9, tolerance = .001),
          if(Xtra) { set.seed(1)
            print(system.time(K3 <- evTestK(Uc, method = "jackknife")))
            all.equal(K3$p.value, 0.0 , tolerance = .001) # ??
          } else TRUE
         )
system.time(T1 <- evTestK(Ug))
system.time(T2 <- evTestK(Ug, method = "asymptotic"))
stopifnot(all.equal(T1$p.value, 0.1363025, tolerance = .001),
          all.equal(T2$p.value, 0.1217252, tolerance = .001),
          if(Xtra) { set.seed(2)
            print(system.time(T3 <- evTestK(Ug, method = "jackknife")))# 3.6 sec
            all.equal(T3$p.value, 0.1507736 , tolerance = .001) #?
          } else TRUE
         )
}%\dont show
}
\keyword{htest}
\keyword{multivariate}
