## HyperLOPIT U2OS novelty

data(hyperLOPITU2OS2018)
fData(hyperLOPITU2OS2018)$markers[fData(hyperLOPITU2OS2018)$markers %in% c("RIBOSOME 40S","RIBOSOME 60S",
                                                                   "CHROMATIN", "NUCLEUS")] <- "unknown"

hlU2OSTagmNoveltyparams <- tagmMcmcTrain_Nov(object = hyperLOPITU2OS2018,
                                             numIter = 20000,
                                             thin = 20,
                                             burnin = 10000,
                                             numChains = 6, 
                                             K_bar = 10)

save(hlU2OSTagmNoveltyparams, file = "hlU2OSTagmNoveltyparams.rda")