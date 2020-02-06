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

source(system.file(package="Matrix", "test-tools-1.R", mustWork=TRUE))
##--> tryCatch.W.E(), showProc.time(), assertError(), relErrV(), ...

##' @title If needed, get file from internet - but do not "error out"
##' @param file
##' @param remoteDIR
##' @param method download method
##' @param mode writing mode,    see ?download.file
##' @param ... potentially further arguments to download.file()
##' @return logical: TRUE if download succeeded
##' @author Martin Maechler (22 Mar 2011)
canGet <- function(file,
                   remoteDIR = "http://copula.r-forge.r-project.org/resources",
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

if(!exists("nCorrDigits", mode="function"))
nCorrDigits <- function(target, current, zeroDigs = 16) {
    stopifnot(zeroDigs >= -log10(.Machine$double.eps))# = 15.65
    RE <- relErrV(target, current)
    r <- -log10(abs(RE))
    r[RE == 0] <- zeroDigs
    r[is.na(RE) | r < 0] <- 0 # no correct digit, when relErr is NA
    r
}

## FIXME?  setTheta() alone should work nowadays
setPar <- function(cop, par) {
    if((is(cop, "tCopula") || is(cop, "tevCopula")) && !cop@df.fixed)
	par <- c(par, df = cop@parameters[length(cop@parameters)])
    setTheta(cop, par, noCheck=TRUE)
}

##' Compare true and numerical derivatives
##' @return max error of d_C_d*() and d_logc_d*() :
comparederiv <- function(cop, u, may.warn=FALSE) {
    c(dCdu = max(abs((copula:::dCdu     (cop, u) -
                      copula:::dCduNumer(cop, u, may.warn=may.warn)))),
      dCdtheta = max(abs(copula:::dCdtheta     (cop, u) -
                         copula:::dCdthetaNumer(cop, u, may.warn=may.warn))),
      dlogcdu = max(abs(copula:::dlogcdu     (cop, u) -
                        copula:::dlogcduNumer(cop, u, may.warn=may.warn))),
      dlogcdtheta = max(abs(copula:::dlogcdtheta     (cop, u) -
                            copula:::dlogcdthetaNumer(cop, u, may.warn=may.warn))))
}
