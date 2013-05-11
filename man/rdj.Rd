\name{rdj}
\Rdversion{1.1}
\alias{rdj}
\docType{data}
\title{Daily Returns of Three Stocks in the Dow Jones}
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
    \item{\code{DATE}}{a character vector specifying the date}
    \item{\code{INTC}}{daily log-return of the Intel stock}
    \item{\code{MSFT}}{daily log-return of the Microsoft stock}
    \item{\code{GE}}{daily log-return of the General Electric}
  }
}
\references{
  McNeil, A. J., Frey, R., and Embrechts, P. (2005).
  \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
  Princeton University Press.
}
\examples{
data(rdj)
}
\keyword{datasets}
