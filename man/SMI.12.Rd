\name{SMI.12}
\title{SMI Data -- 141 Days in Winter 2011/2012}
\alias{SMI.12}
\alias{lSMI}
\docType{data}
\description{
  \code{SMI.12} contains the close prices of all 20 constituents of
  the Swiss Market Index (SMI) from 2011-09-09 to 2012-03-28.
}
\usage{
data(SMI.12)
}
\format{
  \code{SMI.12} is conceptually a multivariate time series,
  here simply stored as \code{\link{numeric}} \code{\link{matrix}},
  where the \code{\link{rownames}} are dates (of week days).

  The format is:

 num [1:141, 1:20] 16.1 15.7 15.7 16.1 16.6 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:141] "2011-09-09" "2011-09-12" "2011-09-13" "2011-09-14" ...
  ..$ : chr [1:20] "ABBN" "ATLN" "ADEN" "CSGN" ...

  ... from 2011-09-09 to 2012-03-28

  \code{lSMI} is the list of the original data (\emph{before} \code{NA}
  \dQuote{imputation}).
}
%\author{Marius Hofert and Martin Maechler}
\source{
  The data was drawn from Yahoo! Finance.
}
\examples{
data(SMI.12)
## maybe
head(SMI.12)

str(D.12 <- as.Date(rownames(SMI.12)))
summary(D.12)

matplot(D.12, SMI.12, type="l", log = "y",
        main = "The 20 SMI constituents  (2011-09 -- 2012-03)",
        xaxt="n", xlab = "2011  /  2012")
Axis(D, side=1)



if(FALSE) { ##--- This worked up to mid 2012, but no longer ---
 begSMI <- "2011-09-09"
 endSMI <- "2012-03-28"
 ##-- read *public* data ------------------------------
 stopifnot(require(zoo), # -> to access all the zoo methods
           require(tseries))
 symSMI <- c("ABBN.VX","ATLN.VX","ADEN.VX","CSGN.VX","GIVN.VX","HOLN.VX",
	     "BAER.VX","NESN.VX","NOVN.VX","CFR.VX", "ROG.VX", "SGSN.VX",
	     "UHR.VX", "SREN.VX","SCMN.VX","SYNN.VX","SYST.VX","RIGN.VX",
	     "UBSN.VX","ZURN.VX")
 lSMI <- sapply(symSMI, function(sym)
		get.hist.quote(instrument = sym, start= begSMI, end= endSMI,
			       quote = "Close", provider = "yahoo",
			       drop=TRUE))
 ## check if stock data have the same length for each company.
 sapply(lSMI, length)
 ## "concatenate" all:
 SMIo <- do.call(cbind, lSMI)
 ## and fill in the NAs :
 SMI.12 <- na.fill(SMIo, "extend")
 colnames(SMI.12) <- sub("\\\\.VX", "", colnames(SMI.12))
 SMI.12 <- as.matrix(SMI.12)
}##----       --- original download

zoo.there <- "package:zoo" \%in\% search()
if(zoo.there || require("zoo")) {
  stopifnot(identical(SMI.12,
     local({ S <- as.matrix(na.fill(do.call(cbind, lSMI), "extend"))
             colnames(S) <- sub("\\\\.VX", "", colnames(S)); S })))
  if(!zoo.there) detach("package:zoo")
}

}
\keyword{datasets}
