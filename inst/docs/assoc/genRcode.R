## function to generate R source defining functions
genApproxFun <- function(fin, fout, fun1name, fun2name) {
  txy <- read.table(fin, head=TRUE)
  vnames <- names(txy)
  cat("## This file is generated. DO NOT EDIT MANUALLY\n", file=fout)
  cat(vnames[1], " <- c(", "\n",
      paste(txy[,1], collapse=", "),
      ")", "\n", file=fout, append=TRUE)
  cat(vnames[2], " <- c(", "\n",
      paste(txy[,2], collapse=", "),
      ")", "\n", file=fout, append=TRUE)
#### splinefun is not good for sparse grid points like this
  ##  cat(fun1name, " <- splinefun(", 
  cat(fun1name, " <- approxfun(",
      "x = ", vnames[1], ", ",
      "y = ",  vnames[2], ", ",
      ##      " method='natural'",
      ")\n",
      file=fout, append=TRUE)
  ##  cat(fun2name, " <- splinefun(",
  if (!is.null(fun2name))
    cat(fun2name, " <- approxfun(",
        "x = ",  vnames[2], ", ",
        "y = ", vnames[1], ", ",
        ##      " method='natural'",
        ")\n",
        file=fout, append=TRUE)
}
