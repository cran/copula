\name{device}
\alias{dev.off.pdf}
\title{Cropping and Font Embedding PDF Device}
\description{
  \code{dev.off.pdf()} is a wrapper of \code{\link{dev.off}()} which is
  meant for closing a pdf device. It also performs
  cropping and font embedding if chosen.
}
\usage{
dev.off.pdf(file="Rplots.pdf", crop=NULL, embedFonts="", ...)
}
\arguments{
  \item{file}{output file name including extension .pdf.}
  \item{crop}{cropping command, can be one of:
    \describe{
      \item{\code{\link{NULL}}}{crop with the command
	\code{"pdfcrop --pdftexcmd pdftex file file 1>/dev/null 2>&1"}.
	This is suitable for Unix; for non-Unix, no cropping is done.}
      \item{\code{\link{character}}}{a string containing the crop
	command.}
      \item{""}{do not crop.}
    }
  }
  \item{embedFonts}{font embedding command, can be one of:
    \describe{
      \item{\code{\link{NULL}}}{embed fonts with the command
	\code{embedFonts(file, options="-dSubsetFonts=true
	  -dEmbedAllFonts=true -dPDFSETTINGS=/printer
	  -dUseCIEColor")}. This is suitable for Unix; for non-Unix,
	  no font embedding is done.}
      \item{\code{\link{character}}}{a string containing a font embedding
	command.}
      \item{""}{do not embed fonts.}
  }}
  \item{\dots}{additional arguments passed to \code{\link{dev.off}()}.}
}
\value{
  \code{\link{invisible}()}.
}
\author{Marius Hofert.}
\seealso{
  \code{\link{dev.off}()} for closing a device,
  \code{\link{embedFonts}()} for font embedding.
  \pkg{sfsmisc}'s \code{\link[sfsmisc]{pdf.end}()} for another approach.
}
\examples{
## typical usage
doPDF <- !dev.interactive(orNone=TRUE)
if(doPDF) pdf(file=(file <- "crop_device.pdf"), width=6, height=6)
plot(1)
if(doPDF) dev.off.pdf(file)
}
\keyword{utilities}
