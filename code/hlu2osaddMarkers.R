## HyperLOPIT U2OS novelty

data(hyperLOPITU2OS2018)
load("endosomeMarkers.rda")
levels(fData(hyperLOPITU2OS2018)$markers) <- factor(c(levels(fData(hyperLOPITU2OS2018)$markers), "endosome"),
                                                    levels = c("CHROMATIN","CYTOSOL","ER","GA","LYSOSOME","MITOCHONDRION","NUCLEUS",
                                                    "PEROXISOME","PM","PROTEASOME","RIBOSOME 40S","RIBOSOME 60S","endosome","unknown" ))
fData(hyperLOPITU2OS2018)[endosomeMarkers,]$markers <- "endosome"  #make sure levels are in right order!!

hlU2OSTagmaddMarkers <- tagmMcmcTrain_Nov(object = hyperLOPITU2OS2018,
                                                   numIter = 20000,
                                                   thin = 20,
                                                   burnin = 10000,
                                                   numChains = 6,
                                                   K_bar = 1)

save(hlU2OSTagmaddMarkers, file = "hlU2OSTagmaddMarkers.rda")
