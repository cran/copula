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


## Catch and save both warnings and errors and in the case of
## a warning, also keep the computed result
tryCatch.W.E <- function(expr){
    W <- NULL
    w.handler <- function(w){ # warning handler
        W <<- w
        invokeRestart("muffleWarning")
    }
    list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
         warning = w.handler), warning = W)
}

##' @title If needed, get file from internet - but do not "error out"
##' @param file
##' @param remoteDIR
##' @param method download method
##' @param mode writing mode,    see ?download.file
##' @param ... potentially further arguments to download.file()
##' @return logical: TRUE if download succeeded
##' @author Martin Maechler (22 Mar 2011)
canGet <- function(file,
                   remoteDIR = "http://nacopula.r-forge.r-project.org/resources",
                   method, mode = "wb", ...)
{
    if(file.exists(file))
        return(TRUE)
    ## else try to down load it
    fullURL <- file.path(remoteDIR, file)
    r <- tryCatch( download.file(fullURL, destfile = file,
                                 method=method, mode=mode, ...),
                  error = function(e) e, warning = function(w) w)
    ok <- !is(r, "condition") && r == 0
    if(!ok && file.exists(file)) ## try to remove a newly created empty file
        tryCatch(file.remove(file), error = function(e){})
    ok
}

if(FALSE)## is now in
    system.file(package="Matrix", "test-tools-1.R")
##' Compute the relative error between target and current vector
##' @title Relative Error (:= 0 when absolute error == 0)
##' @param target
##' @param current
##' @return vector of the same length as target and current
##' @author Martin Maechler
relErrV <- function(target, current) {
    ## relative error, but give 0 when absolute error==0
    ## assert( <length current> is multiple of <length target>) :
    n <- length(target <- as.vector(target))
    if(length(current) %% n)
	stop("length(current) must be a multiple of length(target)")
    RE <- current
    RE[] <- 0
    fr <- current/target
    neq <- is.na(current) | (current != target)
    RE[neq] <- 1 - fr[neq]
    RE
}

##' @title Number of correct digits - recoding "Inf" to 'zeroDigs'
##' @param target  numeric vector of "true" values
##' @param current numeric vector of "approximate" values
##' @param zeroDigs how many correct digits should zero error give
##' @return
##' @author Martin Maechler
nCorrDigits <- function(target, current, zeroDigs = 16) {
    stopifnot(zeroDigs >= -log10(.Machine$double.eps))# = 15.65
    RE <- relErrV(target, current)
    r <- -log10(abs(RE))
    r[RE == 0] <- zeroDigs
    r[is.na(RE) | r < 0] <- 0 # no correct digit, when relErr is NA
    r
}
