## HEK phenodisco/Novelty TAGM comparison

data("HEK293T2011")
fData(HEK293T2011)$pd.markers[fData(HEK293T2011)$pd.markers %in% c("Chromatin associated",
                                                                   "Cytosol",
                                                                   "Cytosol/Nucleus",
                                                                   "Lysosome",
                                                                   "Nucleus")] <- "unknown"

hektagmNov <- tagmMcmcTrain_Nov(object = HEK293T2011,
                                  fcol = "pd.markers",
                                  numIter = 20000,
                                  thin = 20,
                                  burnin = 10000,
                                  numChains = 6,
                                  K_bar = 10)

save(hektagmNov, file = "hektagmNov.rda")
