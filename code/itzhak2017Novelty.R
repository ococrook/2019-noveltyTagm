## itzhak 2017 novelty

data("itzhak2017")
fData(itzhak2017)$markers[fData(itzhak2017)$markers %in% c("Large Protein Complex")] <- "unknown"
itzhak2017 <- filterNA(itzhak2017[,c(1:18)])
normalise(itzhak2017,method = "vsn")

itzhak2017TagmNoveltyparams <- tagmMcmcTrain_Nov(object = itzhak2017,
                                                 numIter = 20000,
                                                 thin = 20,
                                                 burnin = 10000,
                                                 numChains = 6, 
                                                 K_bar = 10)

save(itzhak2017TagmNoveltyparams, file = "itzhak2017TagmNoveltyparams.rda")