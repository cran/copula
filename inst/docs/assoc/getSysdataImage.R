source("../../../R/plackettCopula.asso.R")
source("../../../R/plackettCopula.dTau.PosLogAlp.R")

source("../../../R/claytonCopula.rho.R")
source("../../../R/claytonCopula.rhoDer.R")

source("../../../R/gumbelCopula.rho.R")
source("../../../R/gumbelCopula.rhoDer.R")

save(calibKendallsTauPlackettCopula.tr,
     calibSpearmansRhoClaytonCopula.tr,
     calibSpearmansRhoGumbelCopula.tr,
     kendallsTauDerPosLogAlpPlackettCopula.tr,
     kendallsTauPlackettCopula.tr,
     spearmansRhoClaytonCopula.tr,
     spearmansRhoDerClaytonCopula.tr,
     spearmansRhoDerGumbelCopula.tr,
     spearmansRhoGumbelCopula.tr,
     file = "../../../R/sysdata.rda", compress=TRUE)
