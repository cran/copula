## Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.


## From Matrix package [ ~/R/Pkgs/Matrix/R/Auxiliaries.R ]
chk.s <- function(..., which.call = -1) {
    if(nx <- length(list(...)))
	warning(sprintf(ngettext(nx,
                                 "extra argument %s will be disregarded in\n %s",
                                 "extra arguments %s will be disregarded in\n %s"),
                        sub(")$", '', sub("^list\\(", '', deparse(list(...), control=c()))),
                        deparse(sys.call(which.call), control=c())),
                call. = FALSE, domain=NA)
}


###  'interval'	 class utilities
###  =========================== these are small and simple
###  use require(package= "Intervals") if you want serious interval "work"

interval <- function(ch) {
    ## Purpose: "interval" object constructor from string  "[ a, b)", ...
    ## Author: Martin Maechler, Date: 16 Nov 2009
    stopifnot(is.character(ch), length(ch) == 1L)
    sp <- strsplit(ch, ", *")[[1]]
    if(length(sp) != 2L) stop("'ch' must contain exactly one \",\"")
    L <- gsub(" +", "", sp[1]); bL <- substr(L, 1,1)
    if(!any(iL <- bL == c("(","[","]")))
	stop("interval specification must start with \"[\",  \"(\"  or \"]\"")
    R <- gsub(" +", "", sp[2]); nR <- nchar(R); bR <- substr(R, nR,nR)
    if(!any(iR <- bR == c(")","]","[")))
	stop("interval specification must end with \")\",  \"]\"  or \"[\"")
    new("interval", as.numeric(c(substring(L, 2), substr(R, 1, nR-1))),
	open = c(which(iL) != 2, which(iR) != 2))
}

##' Directly convert numeric (length 2) vector to closed interval
##'
##' @title Closed Interval from Numeric
##' @param x numeric vector of length two
##' @param open logical, of length one or two
##' @return an \code{"interval"} object
##' @author Martin Maechler
num2interval <- function(x, open = FALSE) {
    stopifnot(is.numeric(x), length(x) == 2,
	      is.logical(open), 1 <= (lo <- length(open)), lo <= 2)
    new("interval", as.numeric(x), open = rep(open, length.out=2)[1:2])
}
setAs("numeric", "interval", function(from) num2interval(from))

setMethod("format", "interval",
	  function(x, trim = TRUE, ...) {
    r <- format(x@.Data, trim=trim, ...)
    paste(if(x@open[1]) "(" else "[", r[1],", ", r[2],
	  if(x@open[2]) ")" else "]", sep="")
})

setMethod("show", "interval",
	  function(object) cat("'interval' object  ", format(object), "\n",
	sep=''))

##' Summary group method: range(), min(), max(), [sum(), prod(), any(), all()] :
setMethod("Summary", signature(x = "interval", na.rm = "ANY"),
	  function(x, ..., na.rm) callGeneric(x@.Data, ..., na.rm=na.rm))

setMethod("%in%", signature(x = "numeric", table = "interval"),
	  function(x, table) {
	      op <- table@open
	      (if(op[1]) `<` else `<=`)(table[1], x) &
	      (if(op[2]) `<` else `<=`)(x, table[2])
	  })

##' Not exported .. to be used for fast default checks ==> in ../tests/ and ../man/
##' Now can be numeric via (shell)   export R_copula_check_extra=2.5
doExtras <- function() {
    nz <- FALSE
    if(interactive() || (nz <- nzchar(copX <- Sys.getenv("R_copula_check_extra"))) ||
       identical("true", unname(Sys.getenv("R_PKG_CHECKING_doExtras"))))
        if(nz) { # ==> cX is  logical or numeric
            if(!is.na(cX <- as.logical(copX))) cX else as.numeric(copX)
        } else TRUE
    else FALSE
}

## ===> ../man/corKendall.Rd
##' @title (Fast) Computation of pairwise Kendall's taus
##' @param x data, a n x p matrix (or less efficiently a data.frame).
##' @param checkNA logical indicating if 'x' should be checked for NAs ....
##' @param use a string to determine the treatment of NA's,...
##' @return The p x p matrix K of pairwise Kendall's taus, K[i,j] := tau(x[,i], x[,j])
##' @author Martin Maechler
corKendall <- function(x, checkNA = TRUE,
                       use = if(checkNA && anyNA(x)) "pairwise" else "everything") {
    if(use != "everything")
        cor(x, method="kendall", use = use)
    else cor.fk(x) # from package pcaPP
}

##' format() a 'call' -- used for print(<fitCopula>) and print(<fitMvdc>):
formatCall <- function(cal, className, sep = "\n", collapse = "\n") {
    if(cal[[1L]] == as.symbol(".local"))
	cal[[1L]] <- as.symbol(className)
    if(names(cal[2L]) == "copula")
	names(cal)[2L] <- ""
    paste(deparse(cal), sep=sep, collapse=collapse)
}

## describeCop() -- generic in ./AllGeneric.R
## =============

setMethod("describeCop", c("Copula", "missing"), # -> default kind = "short"
	  function(x, kind, prefix="", ...) describeCop(x, "short", prefix, ...))

## FIXME: This table can be extended to cover more atomic copulas
.copulaNameTab <- matrix(c("galambosCopula",     "Galambos",
                           "gumbelCopula",       "Gumbel",
                           "huslerReissCopula",  "Husler-Reiss",
                           "tawnCopula",         "Tawn",
                           "tevCopula",          "t-ev"), byrow = TRUE, ncol = 2)


setMethod("describeCop", c("copula", "character"),
          function(x, kind = c("short", "very short", "long"), prefix = "", ...) {
    kind <- match.arg(kind)
    cl <- class(x)
    if(!is.na(idx <- match(cl, .copulaNameTab[,1]))) cl <- .copulaNameTab[idx, 2]
    if(kind == "very short") # e.g. for show() which has more parts
        return(paste0(prefix, cl, " copula"))
    ## else
    d <- dim(x)
    ch <- paste0(prefix, cl, " copula, dim. d = ", d)
    switch(kind <- match.arg(kind),
           short = ch,
           long = paste0(ch, "\n", prefix, " param.: ",
                         capture.output(str(x@parameters,
                                            give.head=FALSE))),
           stop("invalid 'kind': ", kind))
})

setMethod("describeCop", "xcopula", # "ANY"
	  function(x, kind, prefix = "", ...) {
	      paste(class(x), "copula: ", describeCop(x@copula, kind=kind, prefix=prefix, ...))
	  })

## *Specific* describeCop() methods in the class  ./<copClass>.R  files

## FIXME: dput(*, control = "namedVector")  could do this
dputNamed <- function(x, add.c=FALSE, ...) {
    x <- format(x, ...)
    stopifnot(is.character(x), is.vector(x))
    if(!is.null(nx <- names(x)))
        x <- paste(nx, x, sep = " = ")
    paste0(if(add.c) "c", "(", paste(x, collapse=", "), ")")
  }

##' @title Get c(.) (expression) by differentiating C(.) wrt to u1, u2, .., u<d>
##' @param cdf expression of cdf C(.)
##' @param d dimension
##' @return Expression of pdf c(.)
cdfExpr2pdfExpr <- function(cdf, d) {
    if (is.null(cdf)) return(NULL)
    for (i in seq_len(d))
      cdf <- D(cdf, paste0("u", i))
    cdf
}

## This function uses the algorithmic expressions stored in the class object
## It is used by khoudrajiExplicitCopula, joeCopula, etc.
.ExplicitCopula.algr <- function(u, copula, log=FALSE, algoNm, ...) {
    stopifnot((!is.null(copula@exprdist$cdf)))
    dim <- dim(copula)
    stopifnot(!is.null(d <- ncol(u)), dim == d)

    colnames(u) <- paste0("u", 1:dim)
    u.df <- data.frame(u)
    params <- getTheta(copula, freeOnly = FALSE, named = TRUE)

    target <- c(eval(attr(copula@exprdist, algoNm), c(u.df, params)))
    if(log) log(target) else target
}
