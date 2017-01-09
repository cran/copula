\name{qqplot2}
\alias{qqplot2}
\title{Q-Q Plot with Rugs and Pointwise Asymptotic Confidence Intervals}
\description{
  A Q-Q plot (possibly) with rugs and pointwise approximate (via the
  Central Limit Theorem) two-sided 1-\eqn{\alpha}{alpha} confidence intervals.
}
\usage{
qqplot2(x, qF, log = "", qqline.args = if(log=="x" || log=="y")
        list(untf=TRUE) else list(),
        rug.args = list(tcl=-0.6*par("tcl")),
        alpha = 0.05, CI.args = list(col="gray40"),
        CI.mtext = list(text=paste0("Pointwise asymptotic ", 100*(1-alpha),
                        "\% confidence intervals"), side=4,
                        cex=0.6*par("cex.main"), adj=0, col="gray40"),
        main = quote(bold(italic(F)~~"Q-Q plot")),
        main.args = list(text=main, side=3, line=1.1, cex=par("cex.main"),
                         font=par("font.main"), adj=par("adj"), xpd=NA),
        xlab = "Theoretical quantiles", ylab = "Sample quantiles",
        file="", width=6, height=6, \dots)
}
\arguments{
  \item{x}{\code{\link{numeric}}.}
  \item{qF}{(theoretical) quantile function against which the Q-Q plot
    is created.}
  \item{log}{\code{\link{character}} string indicating whether log-scale should be
    used; see \code{?plot.default}.}
  \item{qqline.args}{argument \code{\link{list}} passed to
    \code{\link{qqline}()} for creating the Q-Q line. Use
    \code{qqline.args=NULL} to omit the Q-Q line.}
  \item{rug.args}{argument \code{\link{list}} passed to
    \code{\link{rug}()} for creating the rugs. Use
    \code{rug.args=NULL} to omit the rugs.}
  \item{alpha}{significance level.}
  \item{CI.args}{argument \code{\link{list}} passed to
    \code{\link{lines}()} for plotting the confidence intervals. Use
    \code{CI.args=NULL} to omit the confidence intervals.}
  \item{CI.mtext}{argument \code{\link{list}} passed to
    \code{\link{mtext}()} for plotting information about the confidence
    intervals. Use \code{CI.mtext=NULL} to omit the information.}
  \item{main}{title (can be an expression; use "" for no title).}
  \item{main.args}{argument \code{\link{list}} passed to
    \code{\link{mtext}()} for plotting the title. Use
    \code{main.args=NULL} to omit the title.}
  \item{xlab}{x axis label.}
  \item{ylab}{y axis label.}
  \item{file}{file name including the extension \dQuote{.pdf}.}
  \item{width}{width parameter of \code{\link{pdf}()}.}
  \item{height}{height parameter of \code{\link{pdf}()}.}
  \item{\dots}{additional arguments passed to \code{\link{plot}()} based
    for plotting the points.}
}
\value{
  \code{\link{invisible}()}.
}
\details{
  See the source code for how the confidence intervals are constructed precisely.
}
%\author{Marius Hofert}
\seealso{
  \code{\link{plot}()} for the underlying plot function,
  \code{\link{qqline}()} for how the Q-Q line is implemented,
  \code{\link{rug}()} for how the rugs are constructed,
  \code{\link{lines}()} for how the confidence intervals are drawn, and
  \code{\link{mtext}()} for how the title and information about the
  confidence intervals is printed.
  \code{\link{pdf}()} for plotting to pdf.
}
\examples{
n <- 250
df <- 7
set.seed(1)
x <- rchisq(n, df=df)

## Q-Q plot against the true quantiles (of a chi^2_3 distribution)
qqplot2(x, qF = function(p) qchisq(p, df=df),
        main = substitute(bold(italic(chi[NU])~~"Q-Q Plot"), list(NU=df)))

## in log-log scale
qqplot2(x, qF = function(p) qchisq(p, df=df), log="xy",
        main = substitute(bold(italic(chi[NU])~~"Q-Q Plot"), list(NU=df)))

## Q-Q plot against wrong quantiles (of an Exp(1) distribution)
qqplot2(x, qF=qexp, main = quote(bold(Exp(1)~~"Q-Q Plot")))
}
\keyword{hplot}
