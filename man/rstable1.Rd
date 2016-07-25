\name{rstable1}
\title{Random numbers from (Skew) Stable Distributions}
\alias{rstable1}
\alias{rstable}
\description{
  Generate random numbers of the stable distribution
  \deqn{S(\alpha, \beta, \gamma, \delta; k)} with characteristic
  exponent \eqn{\alpha\in(0,2]}{alpha in (0,2]},
  skewness \eqn{\beta\in[-1,1]}{beta in [-1,1]},
  scale    \eqn{\gamma\in[0,\infty)}{gamma in [0,Inf)}, and
  location \eqn{\delta\in\mathbf{R}}{delta in IR}; see Nolan (2010) for
  the parameterization \eqn{k\in\{0,1\}}{k in {0,1}}.  The case
  \eqn{\gamma=0} is understood as the unit jump at \eqn{\delta}.
}
\usage{
rstable1(n, alpha, beta, gamma = 1, delta = 0, pm = 1)
}
\arguments{
  \item{n}{an \code{\link{integer}}, the number of observations to generate.}
    \item{alpha}{characteristic exponent \eqn{\alpha\in(0,2]}{alpha in (0,2]}.}
    \item{beta}{skewness \eqn{\beta\in[-1,1]}{beta in [-1,1]}.}
    \item{gamma}{scale \eqn{\gamma\in[0,\infty)}{gamma in [0,Inf)}.}
    \item{delta}{location \eqn{\delta\in\mathbf{R}}{delta in IR}.}
    \item{pm}{0 or 1, denoting which parametrization (as by Nolan) is used.}
}
\value{
  A \code{\link{numeric}} vector of length \code{n} containing the
  generated random variates.
}
\details{
  We use the approach of John Nolan for generating random variates of
  stable distributions.  The function \code{rstable1} provides two basic
  parametrizations, by default,

  \code{pm = 1}, the so called \dQuote{S}, \dQuote{S1}, or \dQuote{1}
  parameterization.  This is the parameterization used by Samorodnitsky and
  Taqqu (1994), and is a slight modification of Zolotarev's (A)
  parameterization.  It is the form with the most simple form of the
  characteristic function; see Nolan (2010, p. 8).

  \code{pm = 0} is the \dQuote{S0} parameterization: based on the (M)
  representation of Zolotarev for an alpha stable distribution with
  skewness beta.  Unlike the Zolotarev (M) parameterization, gamma and
  delta are straightforward scale and shift parameters.  This
  representation is continuous in all 4 parameters.
}
%% \author{
%%   Diethelm Wuertz wrote \code{\link[stabledist]{rstable}} (now in
%%   package \pkg{stabledist}) for Rmetrics;
%%   Martin Maechler vectorized it (also in \code{alpha},\dots), fixed it
%%   for \eqn{\alpha=1,\beta\ne 0}{alpha = 1, beta != 0} and sped it up.
%% }
\seealso{\code{\link[stabledist]{rstable}} which also allows the
  2-parametrization and provides further functionality for
  stable distributions.
}
\references{
  Chambers, J. M., Mallows, C. L., and Stuck, B. W. (1976),
  \emph{A Method for Simulating Stable Random Variables},
  J.  Amer.  Statist.  Assoc. \bold{71}, 340--344.

  Nolan, J. P. (2012),
  \emph{Stable Distributions---Models for Heavy Tailed Data},
  Birkhaeuser, in progress.

  Samoridnitsky, G. and Taqqu, M. S. (1994),
  \emph{Stable Non-Gaussian Random Processes, Stochastic Models
    with Infinite Variance},
  Chapman and Hall, New York.

% Nolan, J.P. (1999),
% \emph{Stable Distributions},
% Preprint, University Washington DC, 30 pages.

% Weron, A., Weron R. (1999),
%     \emph{Computer Simulation of Levy alpha-Stable Variables and
%     Processes},
%     Preprint Technical Univeristy of Wroclaw, 13 pages.
}
\examples{
   # Generate and plot a series of stable random variates
   set.seed(1953)
   r <- rstable1(n = 1000, alpha = 1.9, beta = 0.3)
   plot(r, type = "l", main = "stable: alpha=1.9 beta=0.3",
        col = "steelblue"); grid()

   hist(r, "Scott", prob = TRUE, ylim = c(0,0.3),
        main = "Stable S(1.9, 0.3; 1)")
   lines(density(r), col="red2", lwd = 2)
}
\keyword{distribution}
