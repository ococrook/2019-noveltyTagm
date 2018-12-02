## HyperLOPIT MOUSE novelty

data("hyperLOPIT2015")
fData(hyperLOPIT2015)$markers[fData(hyperLOPIT2015)$markers %in% c("40S Ribosome","60S Ribosome",
                                                                "Nucleus - Chromatin", "Nucleus - Non-chromatin")] <- "unknown"

hlmouseTagmNoveltyparams <- tagmMcmcTrain_Nov(object = hyperLOPIT2015,
                                             numIter = 20000,
                                             thin = 20,
                                             burnin = 10000,
                                             numChains = 6, 
                                             K_bar = 10)

save(hlmouseTagmNoveltyparams, file = "hlmouseTagmNoveltyparams.rda")