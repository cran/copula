\name{pairsRosenblatt}
\alias{pairsRosenblatt}
\alias{pairsColList}
\title{Plots for Graphical GOF Test via Pairwise Rosenblatt Transforms}
\description{
  \code{pairsColList()} creates a \code{\link{list}} containing
  information about colors for a given matrix of (approximate aka
  \dQuote{pseudo}) p-values.  These colors are used in
  \code{pairsRosenblatt()} for visualizing a graphical
  goodness-of-fit test based on pairwise Rosenblatt transformed data.
}
\usage{
pairsRosenblatt(cu.u, pvalueMat=pviTest(pairwiseIndepTest(cu.u)),
                method = c("scatter", "QQchisq", "QQgamma",
                "PPchisq", "PPgamma", "none"),
                g1, g2, col = "B&W.contrast",
                colList = pairsColList(pvalueMat, col=col),
                main=NULL,
                sub = gpviString(pvalueMat, name = "pp-values"),
		panel = NULL, do.qqline = TRUE,
                keyOpt = list(title="pp-value", rug.at=pvalueMat), ...)

pairsColList(P, pdiv = c(1e-04, 0.001, 0.01, 0.05, 0.1, 0.5),
             signif.P = 0.05, pmin0 = 1e-05, bucketCols = NULL,
             fgColMat = NULL, bgColMat = NULL, col = "B&W.contrast",
             BWcutoff = 170,
             bg.col = c("ETHCL", "zurich", "zurich.by.fog", "baby",
                        "heat", "greenish"),
             bg.ncol.gap = floor(length(pdiv)/3),
             bg.col.bottom = NULL, bg.col.top = NULL, ...)
}
\arguments{
  \item{cu.u}{\eqn{(n,d,d)}-\code{\link{array}} of pairwise
    Rosenblatt-transformed observations as returned by
    \code{\link{pairwiseCcop}()}.}
  \item{pvalueMat}{\eqn{(d,d)}-\code{\link{matrix}} of p-values (or pp-values).}
  \item{method}{\code{\link{character}} indicating the plot method
    to be used. Currently possible are:
    \describe{
      \item{\code{"scatter"}}{a simple scatter plot.}
      \item{\code{"QQchisq"}}{a Q-Q plot after a map to the
        \eqn{\chi^2}{chi-square} distribution.}
      \item{\code{"QQgamma"}}{a Q-Q plot after a map to the
        gamma distribution.}
      \item{\code{"PPchisq"}}{a P-P plot after a map to the
        \eqn{\chi^2}{chi-square} distribution.}
      \item{\code{"PPgamma"}}{a P-P plot after a map to the
        gamma distribution.}
      \item{\code{"none"}}{no points are plotted.}
    }
    Note: These methods merely just set \code{g1} and \code{g2}
    correctly; see the code for more details.
  }
  \item{g1}{\code{\link{function}} from
    \eqn{[0,1]^n\to[0,1]^n}{[0,1]^n -> [0,1]^n}
    applied to "x" for plotting in one panel.}
  \item{g2}{\code{\link{function}} from
    \eqn{[0,1]^{n\times 2}\to[0,1]^n}{[0,1]^(n x 2) -> [0,1]^n}
    applied to "y" for plotting in one panel.}
  \item{colList}{\code{\link{list}} of colors and information as
    returned by \code{pairsColList()}.}
  \item{main}{title.}
  \item{sub}{sub-title with a smart default containing a global (p)p-value.}
  \item{panel}{a \code{panel} function as for \code{\link{pairs}}, or,
    by default, \code{NULL}, where the panel is set as \code{\link{points}}
    or \dQuote{\code{points} \eqn{+} \code{qqline}} if the method is
    "QQ...." and \code{do.qqline} is true.}
  \item{do.qqline}{if \code{method = "QQ...."}, specify if the plot
    panels should also draw a \code{\link{qqline}()}.}
  \item{keyOpt}{argument passed to \code{\link{.pairsCond}()} for
    options for the key.}
  \item{\dots}{additional arguments passed to \code{\link{.pairsCond}()}
    (for \code{pairsRosenblatt()})
    and to \code{heat_hcl()} (for \code{pairsColList};
    used to generate the color palette), see Details.}
  \item{P}{\eqn{d \times d}{d * d} \code{\link{matrix}} of p-values.}
  \item{pdiv}{numeric vector of strictly increasing p-values in
    (0,1) that determine the \dQuote{buckets} for the background
    colors of \code{\link{.pairsCond}()} which creates the
    pairs-like goodness-of-fit plot.}
  \item{signif.P}{significance level (must be an element of \code{pdiv}).}
  \item{pmin0}{a \code{\link{numeric}} indicating the lower endpoint of
    the p-value buckets if \code{pmin} is zero.  If set to 0, the lowest
    value of the p-value buckets will also be 0.

    Note that \code{pmin0} should be in (0, \code{min(pdiv)}) when using
    \code{pairsColList()} for \code{\link{.pairsCond}()}.}
  \item{bucketCols}{\code{\link{vector}} of length as \code{pdiv} containing the
    colors for the buckets.  If not specified, either
    \code{bg.col.bottom} and \code{bg.col.top} are used (if provided) or
    \code{bg.col} (if provided).}
  \item{fgColMat}{\eqn{(d,d)}-\code{\link{matrix}} with foreground colors (the
    default will be black if the background color is bright
    and white if it is dark; see also \code{BWcutoff}).}
  \item{bgColMat}{\eqn{(d,d)}-\code{\link{matrix}} of background colors;
    do not change this unless you know what you are doing.}
  \item{col}{foreground color (defaults to "B&W.contrast" which switches
    black/white according to \code{BWcutoff}), passed to
    \code{\link{.pairsCond}()}.  If \code{colList} is not specified,
    this color is used to construct the points' color.}
  \item{BWcutoff}{number in (0, 255) for switching
    foreground color if \code{col="B&W.contrast"}.}
  \item{bg.col}{color scheme for the background colors.}
  \item{bg.ncol.gap}{number of colors left out as "gap" for
    color buckets below/above \code{signif.P} (to make
    significance/non-significance more visible).}
  \item{bg.col.bottom}{\code{\link{vector}} of length 3 containing a HCL
    color specification. If \code{bg.col.bottom} is provided and
    \code{bucketCols} is not, \code{bg.col.bottom} is used as the color
    for the bucket of smallest p-values.}
  \item{bg.col.top}{\code{\link{vector}} of length 3 containing a HCL
    color specification. If \code{bg.col.top} is provided and
    \code{bucketCols} is not, \code{bg.col.top} is used as the color for
    the bucket of largest p-values.}
}
\value{
\describe{
  \item{\code{pairsRosenblatt}: }{invisibly, the result of \code{\link{.pairsCond}()}.}
  \item{\code{pairsColList}: }{a named \code{\link{list}} with components
    \describe{
      \item{\code{fgColMat}}{\code{\link{matrix}} of foreground colors.}
      \item{\code{bgColMat}}{\code{\link{matrix}} of background colors
	(corresponding to \code{P}).}
	\item{\code{bucketCols}}{\code{\link{vector}} containing the colors
        corresponding to \code{pvalueBuckets} as described above.}
      \item{\code{pvalueBuckets}}{\code{\link{vector}} containing the
	endpoints of the p-value buckets.}
    }
  }
}
}
\details{
  Extra arguments of \code{pairsRosenblatt()} are passed to
  \code{\link{.pairsCond}()}, these notably may include \code{key}, true
  by default, which draws a color key for the colors used as panel
  background encoding (pseudo) p-values.

  \code{pairsColList()} is basically an auxiliary function to specify
  the colors used in the graphical goodness-of-fit test as
  conducted by \code{pairsRosenblatt()}. The latter is described in
  detail in Hofert and \enc{Mächler}{Maechler} (2013).
  See also \code{\link{demo}(gof_graph)}.
}
\author{Marius Hofert, Martin \enc{Mächler}{Maechler}.}
\references{
  Hofert, M. and \enc{Mächler}{Maechler}, M. (2013)
  A graphical goodness-of-fit test for dependence models in higher dimensions;
  submitted.
}
\seealso{
  \code{\link{pairwiseCcop}} for the tools behind the scenes.
  \code{\link{demo}(gof_graph)} for examples.
}
\examples{
## 2-dim example {d = 2} ===============
##
## "t" Copula with 22. degrees of freedom; and (pairwise) tau = 0.5
nu <- 2.2 # degrees of freedom
## Define the multivariate distribution
tCop <- ellipCopula("t", param=iTau(ellipCopula("t", df=nu), tau = 0.5),
                    dim=2, df=nu)
set.seed(19)
X <- qexp(rCopula(n = 400, tCop))

## H0 (wrongly): a Normal copula, with correct tau
copH0 <- ellipCopula("normal", param=iTau(ellipCopula("normal"), tau = 0.5))

## create array of pairwise copH0-transformed data columns
cu.u <- pairwiseCcop(pobs(X), copH0)

## compute pairwise matrix of p-values and corresponding colors
pwIT <- pairwiseIndepTest(cu.u, N=200) # (d,d)-matrix of test results

round(pmat <- pviTest(pwIT), 3) # pick out p-values
## .286 and .077
pairsRosenblatt(cu.u, pvalueMat= pmat)



### A shortened version of   demo(gof_graph) -------------------------------

N <- 32 ## too small, for "testing"; realistically, use a larger one:
if(FALSE)
N <- 100

## 5d Gumbel copula ##########

n <- 250 # sample size
d <- 5 # dimension
family <- "Gumbel" # copula family
tau <- 0.5
set.seed(17)
## define and sample the copula (= H0 copula), build pseudo-observations
cop <- getAcop(family)
th <- cop@iTau(tau) # correct parameter value
copH0 <- onacopulaL(family, list(th, 1:d)) # define H0 copula
U. <- pobs(rCopula(n, cop=copH0))

## create array of pairwise copH0-transformed data columns
cu.u <- pairwiseCcop(U., copH0)

## compute pairwise matrix of p-values and corresponding colors
pwIT <- pairwiseIndepTest(cu.u, N=N, verbose=interactive()) # (d,d)-matrix of test results
round(pmat <- pviTest(pwIT), 3) # pick out p-values
## Here (with seed=1):  no significant ones, smallest = 0.0603

## Plots ---------------------

## plain (too large plot symbols here)
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".")

## with title, no subtitle
pwRoto <- "Pairwise Rosenblatt transformed observations"
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", main=pwRoto, sub=NULL)

## two-line title including expressions, and centered
title <- list(paste(pwRoto, "to test"),
              substitute(italic(H[0]:C~~bold("is Gumbel with"~~tau==tau.)),
                         list(tau.=tau)))
main.line <- c(4, 1.4)
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".",
                main=title, main.line=main.line, main.centered=TRUE)

## Q-Q plots -- can, in general, better detect outliers
pairsRosenblatt(cu.u, pvalueMat=pmat, method="QQchisq", cex=0.2)
}
\keyword{hplot}
\keyword{multivariate}
