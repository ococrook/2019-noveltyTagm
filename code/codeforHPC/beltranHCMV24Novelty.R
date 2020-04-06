## Beltran HCMV 24 novelty

data("beltran2016HCMV24")

beltranHCMV24TagmNoveltyparams <- tagmMcmcTrain_Nov(object = beltran2016HCMV24,
                                                    numIter = 20000,
                                                    thin = 20,
                                                    burnin = 10000,
                                                    numChains = 6, 
                                                    K_bar = 10)

save(beltranHCMV24TagmNoveltyparams, file = "beltranHCMV24TagmNoveltyparams.rda")