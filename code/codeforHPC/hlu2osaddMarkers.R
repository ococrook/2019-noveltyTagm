## HyperLOPIT U2OS novelty

data(hyperLOPITU2OS2018)
load("endosomeMarkers.rda")
fData(hyperLOPITU2OS2018)$nt.markers <- fData(hyperLOPITU2OS2018)$markers
levels(fData(hyperLOPITU2OS2018)$nt.markers) <- factor(c(levels(fData(hyperLOPITU2OS2018)$markers), "endosome"),
                                                    levels = c("CHROMATIN","CYTOSOL","ER","GA","LYSOSOME","MITOCHONDRION","NUCLEUS",
                                                    "PEROXISOME","PM","PROTEASOME","RIBOSOME 40S","RIBOSOME 60S","endosome","unknown" ))
fData(hyperLOPITU2OS2018)[endosomeMarkers,]$nt.markers <- "endosome"  #make sure levels are in right order!!

hlU2OSTagmaddMarkers <- tagmMcmcTrain(object = hyperLOPITU2OS2018,
                                                numIter = 20000,
                                                thin = 20,
                                                burnin = 10000,
                                                numChains = 6,
                                                fcol = "nt.markers")

save(hlU2OSTagmaddMarkers_2, file = "hlU2OSTagmaddMarkers_2.rda")
