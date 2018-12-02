## itzhak 2016 novelty

data("itzhak2016stcSILAC")
fData(itzhak2016stcSILAC)$markers[fData(itzhak2016stcSILAC)$markers %in% c("Large Protein Complex")] <- "unknown"
itzhak2016stcSILAC <- filterNA(itzhak2016stcSILAC)
normalise(itzhak2016stcSILAC, method = "vsn")

itzhak2016TagmNoveltyparams <- tagmMcmcTrain_Nov(object = itzhak2016stcSILAC,
                                            numIter = 20000,
                                            thin = 20,
                                            burnin = 10000,
                                            numChains = 6, 
                                            K_bar = 10)

save(itzhak2016TagmNoveltyparams, file = "itzhak2016TagmNoveltyparams.rda")