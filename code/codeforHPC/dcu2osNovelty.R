## LOPIT-DC U2OS novelty

data(lopitdcU2OS2018)
fData(lopitdcU2OS2018)$markers[fData(lopitdcU2OS2018)$markers %in% c("NUCLEUS/CHROMATIN","PROTEASOME",
                                                                           "RIBOSOME")] <- "unknown"

dcU2OSTagmNoveltyparams <- tagmMcmcTrain_Nov(object = lopitdcU2OS2018,
                                            numIter = 20000,
                                            thin = 20,
                                            burnin = 10000,
                                            numChains = 6, 
                                            K_bar = 10)

save(dcU2OSTagmNoveltyparams, file = "dcU2OSTagmNoveltyparams.rda")