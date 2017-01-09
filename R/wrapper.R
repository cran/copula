## Copyright (C) 2012 Marius Hofert and Martin Maechler
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


### Wrappers for dealing with elliptical (Gauss, t_nu) and Archimedean copulas

##' @title Copula class for the given copula object
##' @param cop copula object
##' @return "ellipCopula" or "outer_nacopula" depending on the given copula object
##' @author Marius Hofert
copClass <- function(cop)
{
    cls <- class(cop)
    if(is(cop, "copula") && (cls=="normalCopula" || cls=="tCopula")) "ellipCopula" # note: there could be other "copula" objects which are not elliptical
    else if(cls=="outer_nacopula") "outer_nacopula" # can be Archimedean or nested Archimedean
    else stop("not yet supported copula object")
}

##' @title Copula family for the given copula object
##' @param cop copula object (either elliptical or (nested) Archimedean)
##' @return family string
##' @author Marius Hofert and Martin Maechler
copFamily <- function(cop)
{
    cls <- getClass(class(cop)) # so extends( . , "..")  is efficient
    if(extends(cls, "copula")) {
        if(extends(cls, "normalCopula")) "normal"
        else if(extends(cls, "tCopula")) "t"
        else stop("unsupported copula family")
    } else if(extends(cls, "outer_nacopula")) {
        cop@copula@name # could be nested or not
    } else stop("not yet supported copula object")
}

##' @title Copula family for the given copula object
##' @param cop copula object (either elliptical or (nested) Archimedean)
##' @return family string
##' @author Marius Hofert
copFamilyClass <- function(family)
{
    if(family == "normal" || family == "t")
	"ellipCopula"
    else if(family %in% .ac.longNames ||
	    family %in% paste0("opower:", .ac.longNames))
	"outer_nacopula" # note: opower not really supported yet
    else stop("family ", family, " not yet supported")
}

