expr2R <- function(fname) {
  myexpr <- readLines(fname)
  myexpr <- gsub("Log", "log", myexpr)
  myexpr <- sub("List\\(", "", myexpr)
  myexpr <- sub("\\)$", "", myexpr)
  myexpr <- strsplit(myexpr, ",")[[1]]
  myexpr <- parse(text=myexpr)
#  myexpr <- sapply(myexpr, function(x) deriv(x, "s"))
  myexpr
}




## deriv(x, "s"): in fact pdf.expr doesn't contain "s"; just a trick to
## get the algorithmic expression for pdf, not its derivative
expr2algr2dump <- function(cname) {
  pdf.expr.name <- paste(cname, "Copula.pdf.expr", sep="")
  pdf.algr.name <- sub("expr", "algr", pdf.expr.name)
  genfun.expr.name <- paste(cname, "Copula.genfun.expr", sep="")
  genfun.algr.name <- sub("expr", "algr", genfun.expr.name)
  assign(pdf.expr.name, pdf.expr <- expr2R(pdf.expr.name))
  assign(pdf.algr.name, sapply(pdf.expr, function(x) deriv(x, "s")))
  assign(genfun.expr.name, genfun.expr <- expr2R(genfun.expr.name))
  assign(genfun.algr.name, sapply(genfun.expr, function(x) deriv(x, "u")))
  dname <- paste("../../R/", cname, "Expr.R", sep="")
  dump(c(pdf.expr.name, pdf.algr.name,
         genfun.expr.name, genfun.algr.name),
       file=dname)
}

## this is not working as expected:
##   assign(pdf.algr.name, sapply(as.name(pdf.expr.name), function(x) deriv(x, "s")))

expr2algr2dump("frank")
expr2algr2dump("clayton")
expr2algr2dump("gumbel")
expr2algr2dump("amh")


## frankCopula.pdf.expr <- expr2R("frankCopula.pdf.expr")
## frankCopula.pdf.algr <- sapply(frankCopula.pdf.expr,
##                                function(x) deriv(x, "s"))
## frankCopula.genfun.expr <- expr2R("frankCopula.genfun.expr")
## frankCopula.genfun.algr <- sapply(frankCopula.genfun.expr,
##                                function(x) deriv(x, "s"))

## dump(c("frankCopula.pdf.expr", "frankCopula.pdf.algr",
##        "frankCopula.genfun.expr", "frankCopula.genfun.algr"),
##        file="../../R/frankExpr.R")



## claytonCopula.pdf.expr <- expr2R("claytonCopula.pdf.expr")
## claytonCopula.pdf.algr <- sapply(claytonCopula.pdf.expr,
##                                function(x) deriv(x, "s"))
## dump(c("claytonCopula.pdf.expr", "claytonCopula.pdf.algr"),
##        file="../../R/claytonExpr.R")



## gumbelCopula.pdf.expr <- expr2R("gumbelCopula.pdf.expr")
## gumbelCopula.pdf.algr <- sapply(gumbelCopula.pdf.expr,
##                                function(x) deriv(x, "s"))
## dump(c("gumbelCopula.pdf.expr", "gumbelCopula.pdf.algr"),
##        file="../../R/gumbelExpr.R")


## galambosCopula.expr <- expr2R("galambos.expr")
## galambosCopula.algr <- sapply(galambosCopula.expr,
##                               function(x) deriv(x, "s"))
## names(galambosCopula.expr) <- names(galambosCopula.algr) <- c("cdf", "pdf", "deriv1cdf")
## dump(c("galambosCopula.expr", "galambosCopula.algr"),
##      file="../../R/galambosExpr.R")
