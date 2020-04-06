## Beltran HCMV 24 novelty

data("beltran2016MOCK24")

beltranMOCK24TagmNoveltyparams <- tagmMcmcTrain_Nov(object = beltran2016MOCK24,
                                                    numIter = 20000,
                                                    thin = 20,
                                                    burnin = 10000,
                                                    numChains = 6, 
                                                    K_bar = 10)

save(beltranMOCK24TagmNoveltyparams, file = "beltranMOCK24TagmNoveltyparams.rda")