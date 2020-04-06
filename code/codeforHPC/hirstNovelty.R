## hirst novelty

data("hirst2018")
fData(hirst2018)$markers[fData(hirst2018)$markers %in% c("Large Protein Complex")] <- "unknown"

hirstTagmNoveltyparams <- tagmMcmcTrain_Nov(object = hirst2018[,c(11:15, 26:30, 41:45)],
                                                    numIter = 20000,
                                                    thin = 20,
                                                    burnin = 10000,
                                                    numChains = 6, 
                                                    K_bar = 10)

save(hirstTagmNoveltyparams, file = "hirstTagmNoveltyparams.rda")