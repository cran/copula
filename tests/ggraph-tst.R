require(copula)
doX <- FALSE # no "doExtras" -- be fast

if(!dev.interactive(orNone=TRUE)) pdf("ggraph-tst.pdf")

demo(gof_graph)

dev.off()
