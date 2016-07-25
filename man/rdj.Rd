\name{rdj}
\title{Daily Returns of Three Stocks in the Dow Jones}
\alias{rdj}
\docType{data}
\description{
  Five years of daily log-returns (from 1996 to 2000) of Intel
  (INTC), Microsoft (MSFT) and General Electric (GE) stocks. These
  data were analysed in Chapter 5 of McNeil, Frey and
  Embrechts (2005).
}
\usage{data(rdj)}
\format{
  A data frame of 1262 daily log-returns from 1996 to 2000.
  \describe{
    \item{\code{Date}}{the date, of class \code{"\link{Date}"}.}
    \item{\code{INTC}}{daily log-return of the Intel stock}
    \item{\code{MSFT}}{daily log-return of the Microsoft stock}
    \item{\code{GE}}{daily log-return of the General Electric}
  }
}
%\author{Marius Hofert, Ivan Kojadinovic, Martin Maechler and Jun Yan}
\references{
  McNeil, A. J., Frey, R., and Embrechts, P. (2005).
  \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
  Princeton University Press.
}
\examples{
data(rdj)
str(rdj)# 'Date' is of class "Date"

with(rdj, {
   matplot(Date, rdj[,-1], type = "o", xaxt = "n", ylim = .15* c(-1,1),
           main = paste("rdj - data;  n =", nrow(rdj)))
   Axis(Date, side=1)
})
legend("top", paste(1:3, names(rdj[,-1])), col=1:3, lty=1:3, bty="n")

\donttest{% save 0.75 sec for fitting  and many more for gof
x <- rdj[, -1] # '-1' : not the Date
## a t-copula (with a vague inital guess of Rho entries)
tCop <- tCopula(rep(.2, 3), dim=3, dispstr="un", df=10, df.fixed=TRUE)
ft <- fitCopula(tCop, data = pobs(x))
ft
ft@copula # the fitted t-copula as tCopula object
system.time(
g.C <- gofCopula(claytonCopula(dim=3), as.matrix(x), simulation = "mult")
) ## 5.3 sec
system.time(
g.t <- gofCopula(ft@copula, as.matrix(x), simulation = "mult")
) ## 8.1 sec

}% dont...
}
\keyword{datasets}
