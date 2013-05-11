\name{evTestK}
\alias{evTestK}
\title{Bivariate Test of Extreme-Value Dependence Based on Kendall's Process}
\description{
  Test of extreme-value dependence based on the  bivariate probability integral
  transformation. The test statistic is defined in Ben Ghorbal,
  G. \enc{Nešlehová}{Neslehova}, and Genest (2009).
}
\usage{
evTestK(x, method = c("fsample","asymptotic","jackknife"))
}
\arguments{
  \item{x}{ a data matrix. }
  \item{method}{ specifies the variance estimation method; can be either
    \code{"fsample"} (finite-sample, the default), \code{"asymptotic"}
    or \code{"jackknife"}.}
}
\details{
  The code for this test was generously provided by Johanna G. \enc{Nešlehová}{Neslehova}.
  More details are available in Appendix B of Ben Ghorbal, G. \enc{Nešlehová}{Neslehova} and
  Genest (2009).
}
\value{
  Returns a list whose attributes are:
  \item{statistic}{ value of the test statistic. }
  \item{p.value}{ corresponding p-value. }
}
\references{
  Ghorbal, M. B., Genest, C., and G. \enc{Nešlehová}{Neslehova}, J. (2009)
  On the test of Ghoudi, Khoudraji, and Rivest for extreme-value dependence.
  \emph{The Canadian Journal of Statistics} \bold{37}, 1--9.
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
\dontrun{
(K3 <- evTestK(Uc, method = "jackknife"))# slow !
}
\dontshow{
Xtra <- copula:::doExtras()
stopifnot(all.equal(K1$p.value, 3.9471315e-9, tol = .001),
          all.equal(K2$p.value, 2.1277296e-9, tol = .001),
          if(Xtra) { set.seed(1)
            print(system.time(K3 <- evTestK(Uc, method = "jackknife")))
            all.equal(K3$p.value, 0.0 , tol = .001) # ??
          } else TRUE
         )
system.time(T1 <- evTestK(Ug))
system.time(T2 <- evTestK(Ug, method = "asymptotic"))
stopifnot(all.equal(T1$p.value, 0.1363025, tol = .001),
          all.equal(T2$p.value, 0.1217252, tol = .001),
          if(Xtra) { set.seed(2)
            print(system.time(T3 <- evTestK(Ug, method = "jackknife")))# 3.6 sec
            all.equal(T3$p.value, 0.0 , tol = .001) #?
          } else TRUE
         )
}%\dont show
}
\keyword{htest}
\keyword{multivariate}
