## mouse HPC script for sensitivity analysis

require(MSnbase)
require(pRoloc)
require(pRolocdata)
require(mcclust)
require(pheatmap)

source("machinelearning-framework-tagm-Novelty.R")
source("machinelearning-functions-tagm-Novelty.R")

data("hyperLOPIT2015")
fData(hyperLOPIT2015)$markers[fData(hyperLOPIT2015)$markers %in% c("40S Ribosome","60S Ribosome",
                                                                   "Nucleus - Chromatin", "Nucleus - Non-chromatin")] <- "unknown"



hlmouseTagmNoveltyparams_beta03 <- tagmMcmcTrain_Nov(object = hyperLOPIT2015,
                                                     fcol = "markers",
                                                     numIter = 15000,
                                                     thin = 10,
                                                     burnin = 5000,
                                                     numChains = 6,
                                                     beta0 = 0.0002,
                                                     K_bar = 10)


save(hlmouseTagmNoveltyparams_beta03, file = "/home/omc25/rds/hpc-work/hlmouseTagmNoveltyparams_beta03.rda")