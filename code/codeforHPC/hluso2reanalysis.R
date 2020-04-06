## hyperLOPIT U2OS with endosomes
setStockcol(paste0(getStockcol(), 90))
load("hlU2OSTagmaddMarkers_2.rda")

data(hyperLOPITU2OS2018)
load("endosomeMarkers.rda")
fData(hyperLOPITU2OS2018)$nt.markers <- fData(hyperLOPITU2OS2018)$markers
levels(fData(hyperLOPITU2OS2018)$nt.markers) <- factor(c(levels(fData(hyperLOPITU2OS2018)$markers), "ENDOSOME"),
                                                       levels = c("CHROMATIN","CYTOSOL","ER","GA","LYSOSOME","MITOCHONDRION","NUCLEUS",
                                                                  "PEROXISOME","PM","PROTEASOME","RIBOSOME 40S","RIBOSOME 60S","ENDOSOME","unknown" ))
fData(hyperLOPITU2OS2018)[endosomeMarkers,]$nt.markers <- "ENDOSOME"  #make sure levels are in right order!!


# check convergence
outliers <- mcmc_get_outliers(hlU2OSTagmaddMarkers_2)

u2osdiag <- gelman.diag(outliers[-c(1,4,6)]) # converged yay
u2ostagmNov_conv <- hlU2OSTagmaddMarkers_2[-c(1,4,6)]

u2osTagmNoveltyRes <- tagmNoveltyProcess(object = hyperLOPITU2OS2018, params = u2ostagmNov_conv, fcol = "nt.markers")

u2ostagmNov_conv <- tagmMcmcProcess(u2ostagmNov_conv)
u2ostagm <- tagmPredict(object = hyperLOPITU2OS2018, params = u2ostagmNov_conv, fcol = "nt.markers")
save(u2osTagmNoveltyRes, file = "u2osTagmNoveltyRes_endo.rda")

## attach information to MSnset
u2ostagm <- tagmNoveltyPredict(object = u2ostagm, params = u2osTagmNoveltyRes)

fData(u2ostagm)$tagm.mcmc.allocation[fData(u2ostagm)$tagm.mcmc.allocation == "endosome"] <- "ENDOSOME"

ptsze <- exp(fData(u2ostagm)$tagm.mcmc.probability) - 1.5
ptsze[fData(u2ostagm)$tagm.mcmc.outlier > 10^{-12}] <- 0.01
ptsze[fData(u2ostagm)$tagm.mcmc.allocation == "Phenotype 1"] <- 0.01
plot2D(u2ostagm, fcol = "tagm.mcmc.allocation", cex = ptsze, main = "PCA of U2OS hyperLOPIT data", cex.main = 2, grid = FALSE)
addLegend(u2ostagm, fcol = "tagm.mcmc.allocation", ncol = 5, where = "topleft", cex = 1)

cls <- getStockcol()[as.factor(fData(u2ostagm)$tagm.mcmc.allocation)]
plot(fData(u2ostagm)$tagm.mcmc.probability,
     fData(u2ostagm)$tagm.mcmc.mean.shannon,
     col = cls, pch = 19,
     xlab = "Localisation probability",
     ylab = "Shannon entropy")

plot(u2ostagmNov_conv, "A6NGN9")

unknownLoc <- rownames(u2ostagm)[(fData(u2ostagm)$tagm.mcmc.probability * (1 - fData(u2ostagm)$tagm.mcmc.outlier) < 0.99)]

mf1u2os_endo <- enrichGO(gene = unknownLoc,
                   OrgDb = "org.Hs.eg.db",
                   keyType = "UNIPROT",
                   ont = "MF",
                   universe = rownames(fData(u2ostagm)))

bp1u2os_endo <- enrichGO(gene = unknownLoc,
                         OrgDb = "org.Hs.eg.db",
                         keyType = "UNIPROT",
                         ont = "BP",
                         universe = rownames(fData(u2ostagm)))
tempNames <- rownames(fData(u2ostagm)[fData(u2ostagm)$tagm.mcmc.allocation == "endosome",])[order(fData(u2ostagm)$tagm.mcmc.mean.shannon[fData(u2ostagm)$tagm.mcmc.allocation == "endosome"], decreasing = T)]

library(patchwork)
## following fixes capitalisation
ch <-chains(u2ostagmNov_conv)
rownames(u2ostagmNov_conv@chains@chains[[1]]@ComponentParam@mk)[rownames(u2ostagmNov_conv@chains@chains[[1]]@ComponentParam@mk) == "endosome"] <- "ENDOSOME"

gg1 <- plot(u2ostagmNov_conv, "Q92738") # Interesting
gg2 <- plot(u2ostagmNov_conv, "Q15833" ) # vescile
gg3 <- plot(u2ostagmNov_conv, "P61020") # known RAB5B
gg4 <- plot(u2ostagmNov_conv, "O15498") # maybe
#  "Q92738" "Q15833" "P61020" "O15498"

tempNames2 <- rownames(fData(u2ostagm)[fData(u2ostagm)$tagm.mcmc.allocation == "PM",])[order(fData(u2ostagm)$tagm.mcmc.mean.shannon[fData(u2ostagm)$tagm.mcmc.allocation == "PM"], decreasing = T)]

gg5 <- plot(u2ostagmNov_conv, tempNames2[3]) # probably not
gg6 <- plot(u2ostagmNov_conv, tempNames2[9]) # traffiking from endosome to membrane
gg7 <- plot(u2ostagmNov_conv, tempNames2[18]) # known RAB
gg8 <- plot(u2ostagmNov_conv, tempNames2[34])  # KIF16B
gg9 <- plot(u2ostagmNov_conv, tempNames2[35]) # known PM, LYS, ENDO
# "O00186"   "Q9NZN3"   "P20339-2" "Q96L93-6" "Q8NHG8" 

gg1 + gg3 + gg4 + gg6 + gg7 + gg8 + gg9 + plot_layout(ncol = 3)


## Barplot of different number of proteins
numMarkers <- sum(fData(hyperLOPITU2OS2018)$markers != "unknown")
allocThres <- (fData(u2ostagm)$tagm.mcmc.probability * (1 - fData(u2ostagm)$tagm.mcmc.outlier) > 0.999)
tb <- table((fData(u2ostagm)$tagm.mcmc.allocation[allocThres]))
numAlloc <- sum(tb[!(names(tb) %in% c("ENDOSOME", "Phenotype 1"))])
numAlloc_withend <- sum(tb[!(names(tb) %in% c("Phenotype 1"))])
totalProteins <- nrow(u2ostagm)

# 240 endosome allocations

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
df <- matrix(NA, ncol = 1, nrow = 4)
names(df) <- c("Markers", "TAGM allocation", "Reannotation allocations")
#df[, 1] <- c(numMarkers, totalProteins - numMarkers, 0, 0)
#df[, 2] <- c(numMarkers, numAlloc - numMarkers, totalProteins - (numAlloc), 0)
df[, 1] <- c(numMarkers, numAlloc - numMarkers, numAlloc_withend - (numAlloc), totalProteins - (numAlloc_withend))
df_long <- melt(df)
#df_long$Var1 <- c("a", "d", "b", "c", "a", "b", "d", "c", "a", "b", "c", "d")
df_long$Var1 <- c("a", "b", "c", "d")
df_long$Var1[df_long$Var1 == "a"] <- c("Markers")
df_long$Var1[df_long$Var1 == "b"] <- c("Protein allocations")
df_long$Var1[df_long$Var1 == "c"] <- c("Reannotation allocations")
df_long$Var1[df_long$Var1 == "d"] <- c("Unknown")
df_long$Var2 <- c("TAGM allocations")
gg <- ggplot(df_long, aes(x = Var2, y = value, fill = Var1, width = 0.5)) + geom_bar(stat="identity", position = position_stack(reverse = TRUE)) #+ coord_flip()
gg <- gg + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                 text = element_text(size=20)) + scale_fill_manual(values=cbPalette, name = "Legend")
gg <- gg + ylab("Number of Proteins ") + xlab("Method") + ggtitle("Protein allocations") 
gg



