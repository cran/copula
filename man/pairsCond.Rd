\name{.pairsCond}
\alias{.pairsCond}
\title{Pairs Plot of a cu.u Object (Internal Use)}
\description{
  \code{.pairsCond()} is an internal function for plotting the pairwise
  Rosenblatt transforms, i.e., the pairwise conditional distributions,
  as returned by \code{\link{pairwiseCcop}()}, via the principal function
  \code{\link{pairsRosenblatt}()}.

  The intention is that \code{\link{pairsRosenblatt}()} be called,
  rather than this auxiliary function.
}
\usage{
.pairsCond(gcu.u, panel = points, colList,
     col = par("col"), bg = par("bg"), labels, ...,
     text.panel = textPanel, label.pos = 0.5,
     cex.labels = NULL, font.labels = 1, gap = 0,
     axes = TRUE, panel.border = TRUE, key = TRUE,
     keyOpt = list(space= 2.5, width= 1.5, axis= TRUE,
                   rug.at= numeric(), title= NULL, line= 5),
     main = NULL, main.centered = FALSE,
     line.main = if(is.list(main)) 5/4*par("cex.main")* rev(seq_along(main)) else 2,
     sub = NULL, sub.centered = FALSE, line.sub = 4)
}
\arguments{
  \item{gcu.u}{(n,d,d)-array of pairwise Rosenblatt-transformed u's as returned by
    \code{\link{pairwiseCcop}()}.}
  \item{panel}{panel function, as for \code{\link{pairs}()}.}
  \item{colList}{list of colors and information as returned by
    \code{\link{pairsColList}()}.}
  \item{col}{\emph{instead} of \code{colList}, specifying the points' color.}
  \item{bg}{\emph{instead} of \code{colList}, specifying the constant
    background color.}
  \item{labels}{pairs() argument; can be missing (in which
    case a suitable default is chosen or can be "none" [or
    something else])}

  \item{...}{further arguments, as for \code{\link{pairs}}.  These are
    passed to \code{panel()}, and \code{\link{axis}}, may also contain
    \code{font.main}, \code{cex.main}, and \code{adj}, for title
    adjustments; further, \code{oma} for modifying the default
    \code{\link{par}("oma")}.}

  \item{text.panel, label.pos, cex.labels, font.labels, gap}{see
    \code{\link{pairs}()}.}
  \item{axes}{logical indicating whether axes are drawn.}
  \item{panel.border}{logical indicating whether a border is drawn
    around the pairs (to mimic the behavior of \code{\link{image}()}).}

  \item{key}{logical indicating whether a color key is drawn.}
  \item{keyOpt}{a \code{\link{list}} of options for the color key;
    \describe{
      \item{\code{space}:}{white space in height of characters in
      inch to specify the the distance of the key to the pairs plot.}
      \item{\code{width}:}{key width in height of characters in inch.}
      \item{\code{axis}:}{logical indicating whether an axis for the
	color key is drawn.}
      \item{\code{rug.at}:}{values where rugs are plotted at the key.}
      \item{\code{title}:}{key title.}
      \item{\code{line}:}{key placement (horizontal distance from color
	key in lines).}
    }
  }
  \item{main}{title}
  \item{main.centered}{logical indicating if the title should be
    centered or not; the default \code{FALSE} centers it according to
    the pairs plot, not the whole plotting region.}
  \item{line.main}{title placement (vertical distance from pairs plot in
    lines).}
  \item{sub}{sub-title}
  \item{sub.centered}{logical indicating if the sub-title
    should be centered or not; see \code{main.centered}.}
  \item{line.sub}{sub-title placement, see \code{line.main}.}
}
%% \value{
%%   invisible()
%% }
\note{based on pairs.default() and
  filled.contour() from R-2.14.1 - used in Hofert and Maechler (2013)
}
\author{Marius Hofert and Martin Maechler}
\seealso{
  \code{\link{pairsRosenblatt}()}, the prinicipal function, calling
  \code{.pairsCond()}.
}
\keyword{hplot}
\keyword{multivariate}
