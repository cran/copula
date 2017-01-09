### --- At the end, when "everything" is defined  ---

fitCopula_methods <- eval(formals(fitCopula_dflt)$method)

.copulaEnv <- new.env(parent = emptyenv(), hash = FALSE)#  e.g., for  once-per-session warnings


if(FALSE) { ## <<--- 2016-07-28 --- The following kills  hasMethod(<genfun>, <signature>)

## the generics for which we may want to have "bail out" methods:
.thisEnv <- environment()# == asNamespace("copula")  but not yet
gg <- local({
    g <- getGenerics(.thisEnv)
    g <- g[g@package == "copula"]
    funs <- lapply(g, function(fnam) get(fnam, mode="function", envir=.thisEnv))
    ## those that have 1st arg 'copula':
    g[vapply(lapply(funs, formals), function(f) names(f[1]) == "copula", logical(1))]
})

## bail out methods for everything --- only called if no method for *sub*class exists:
## i.e., almost always if not yet implemented for 'nacopula':

## a1 <- vapply(gg, function(gn) length(formals(gn)) == 1, logical(1))
## for(gname in gg[a1])
##     setMethod(gname, "Copula", function(copula)
##               stop(gettextf("%s() method for class \"%s\" not yet implemented",
##                             .Generic, class(copula))))
## for(gname in gg[!a1]) { # possibly more than one arg -- use correct argument list
for(gname in gg) { # possibly more than one arg -- use correct argument list
    f <- getDataPart(getGeneric(gname, where=.thisEnv))
    body(f) <- bquote(stop(gettextf("%s() method for class \"%s\" not yet implemented",
                                    .Generic, class(copula))))
    setMethod(gname, "Copula", f)
}
rm(gg, .thisEnv)

}## no longer (2016-07-28)
