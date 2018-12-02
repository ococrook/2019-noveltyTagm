## yeast novelty

data("yeast2018")

hlyeastTagmNoveltyparams <- tagmMcmcTrain_Nov(object = yeast2018,
                                               numIter = 20000,
                                               thin = 20,
                                               burnin = 10000,
                                               numChains = 6, 
                                               K_bar = 10)

save(hlyeastTagmNoveltyparams, file = "hlyeastTagmNoveltyparams.rda")