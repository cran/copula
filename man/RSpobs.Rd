\name{RSpobs}
\alias{RSpobs}
%% still too long: ..
\title{Pseudo-Observations of Radial and Uniform Part of Elliptical and Archimedean Copulas}
\description{
  Given a matrix of iid multivariate data from a meta-elliptical or
  meta-Archimedean model, \code{RSpobs()} computes pseudo-observations
  of the radial part \eqn{R} and the vector \eqn{\bm{S}}{S} which
  follows a uniform distribution on the unit sphere (for elliptical
  copulas) or the unit simplex (for Archimedean copulas). These
  quantities can be used for (graphical) goodness-of-fit tests, for example.
}
\usage{
RSpobs(x, do.pobs = TRUE, method = c("ellip", "archm"), ...)
}
\arguments{
  \item{x}{an \eqn{(n, d)}-\code{\link{matrix}} of data; if
    \code{do.pobs=FALSE}, the rows of \code{x} are assumed to lie in the
    \eqn{d}-dimensional unit hypercube (if they do not, this leads to an
    error).}
  \item{do.pobs}{\code{\link{logical}} indicating whether
    \code{\link{pobs}()} is applied to \code{x} for transforming the
    data to the \eqn{d}-dimensional unit hypercube.}
  \item{method}{\code{\link{character}} string indicating the assumed
    underlying model, being meta-elliptical if \code{method="ellip"} (in
    which case \code{S} should be approximately uniform on the
    \eqn{d}-dimensional unit sphere) or
    meta-Archimedean if \code{method="archm"} (in which case \code{S}
    should be approximately uniform on the \eqn{d}-dimensional unit simplex).}
  \item{\dots}{additional arguments passed to the implemented
    methods. These can be
    \describe{
      \item{\code{method="ellip"}}{\code{qQg()}
	(the quantile \code{\link{function}} of the (assumed)
	distribution function \eqn{G_g}{G[g]} as given in Genest,
	Hofert, G. \enc{Nešlehová}{Neslehova}
	(2014)); if provided, \code{qQg()} is used in the transformation for
	obtaining pseudo-observations of \eqn{R} and \eqn{\bm{S}}{S}
	(see the code for more details).}
      \item{\code{method="archm"}}{\code{iPsi()} (the assumed underlying
	generator inverse); if provided, \code{iPsi()} is used in the
	transformation for obtaining pseudo-observations of \eqn{R} and
	\eqn{\bm{S}}{S} (see the code for more details).}
    }
  }
}
\value{
  A \code{\link{list}} with components \code{R} (an
  \eqn{n}-\code{\link{vector}} containing the pseudo-observations
  of the radial part) and \code{S} (an \eqn{(n,
    d)}-\code{\link{matrix}} containing the pseudo-observations of the
  uniform distribution (on the unit sphere/simplex)).
}
\details{
  The construction of the pseudo-obersvations of the radial part and the
  uniform distribution on the unit sphere/simplex is described in
  Genest, Hofert, G. \enc{Nešlehová}{Neslehova} (2014).
}
%\author{Marius Hofert}
\references{
  Genest, C., Hofert, M., G. \enc{Nešlehová}{Neslehova}, J., (2014).
  Is the dependence Archimedean, elliptical, or what?
  \emph{To be submitted}.
}
\seealso{\code{\link{pobs}()} for computing the \dQuote{classical}
  pseudo-observations.}
\examples{
set.seed(100)
n <- 250 # sample size
d <- 5 # dimension
nu <- 3 # degrees of freedom

## Build a mean vector and a dispersion matrix,
## and generate multivariate t_nu data:
mu <- rev(seq_len(d)) # d, d-1, .., 1
L <- diag(d) # identity in dim d
L[lower.tri(L)] <- 1:(d*(d-1)/2)/d # Cholesky factor (diagonal > 0)
Sigma <- crossprod(L) # pos.-def. dispersion matrix (*not* covariance of X)
X <- rep(mu, each=n) + mvtnorm::rmvt(n, sigma=Sigma, df=nu) # multiv. t_nu data
## note: this is *wrong*: mvtnorm::rmvt(n, mean=mu, sigma=Sigma, df=nu)

## compute pseudo-observations of the radial part and uniform distribution
## once for 1a), once for 1b) below
RS.t    <- RSpobs(X, method="ellip", qQg=function(p) qt(p, df=nu)) # 'correct'
RS.norm <- RSpobs(X, method="ellip", qQg=qnorm) # for testing 'wrong' distribution
stopifnot(length(RS.norm$R) == n, length(RS.t$R) == n,
          dim(RS.norm$S) == c(n,d), dim(RS.t$S) == c(n,d))

## 1) Graphically testing the radial part

## 1a) Q-Q plot of R against the correct quantiles
qqplot2(RS.t$R, qF=function(p) sqrt(d*qf(p, df1=d, df2=nu)),
        main.args=list(text=as.expression(substitute(bold(italic(
                       F[list(d.,nu.)](r^2/d.))~~"Q-Q Plot"),
                       list(d.=d, nu.=nu))),
                       side=3, cex=1.3, line=1.1, xpd=NA))

## 1b) Q-Q plot of R against the quantiles of F_R for a multivariate normal
##     distribution
tit <- as.expression(substitute(bold(italic(chi[d.]) ~~ "Q-Q Plot"), list(d.=d)))
qqplot2(RS.norm$R, qF=function(p) sqrt(qchisq(p, df=d)),
        main.args=list(text=tit, side=3, cex=1.3, line=1.1, xpd=NA))

## 2) Graphically testing the angular distribution

## auxiliary function
qqp <- function(k, Bmat) {
    tit <- as.expression(substitute(plain("Beta")(s1,s2) ~~ bold("Q-Q Plot"),
                                    list(s1 = k/2, s2 = (ncol(Bmat)+1-k)/2)))
    qqplot2(Bmat[,k],
            qF=function(p) qbeta(p, shape1=k/2, shape2=(ncol(Bmat)+1-k)/2),
            main.args=list(text=tit, side=3, cex=1.3, line=1.1, xpd=NA))
}

## 2a) Q-Q plot of the 'correct' angular distribution
##     (Bmat[,k] should follow a Beta(k/2, (d-k)/2) distribution)
Bmat.t <- gofBTstat(RS.t$S)
qqp(1, Bmat=Bmat.t) # k=1
qqp(3, Bmat=Bmat.t) # k=3

## 2b) Q-Q plot of the 'wrong' angular distribution
Bmat.norm <- gofBTstat(RS.norm$S)
qqp(1, Bmat=Bmat.norm) # k=1
qqp(3, Bmat=Bmat.norm) # k=3

## 3) Graphically check independence between radial part and B_1 and B_3

## 'correct' distributions (multivariate t)
plot(pobs(cbind(RS.t$R, Bmat.t[,1])), # k = 1
          xlab=quote(italic(R)), ylab=quote(italic(B)[1]),
          main=quote(bold("Rank plot between"~~italic(R)~~"and"~~italic(B)[1])))
plot(pobs(cbind(RS.t$R, Bmat.t[,3])), # k = 3
	  xlab=quote(italic(R)), ylab=quote(italic(B)[3]),
          main=quote(bold("Rank plot between"~~italic(R)~~"and"~~italic(B)[3])))

## 'wrong' distributions (multivariate normal)
plot(pobs(cbind(RS.norm$R, Bmat.norm[,1])), # k = 1
          xlab=quote(italic(R)), ylab=quote(italic(B)[1]),
          main=quote(bold("Rank plot between"~~italic(R)~~"and"~~italic(B)[1])))
plot(pobs(cbind(RS.norm$R, Bmat.norm[,3])), # k = 3
	  xlab=quote(italic(R)), ylab=quote(italic(B)[3]),
          main=quote(bold("Rank plot between"~~italic(R)~~"and"~~italic(B)[3])))
}
\keyword{utilities}
