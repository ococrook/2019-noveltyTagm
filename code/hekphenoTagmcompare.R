## Comparison Phenodisco and Novelty TAGM
setStockcol(paste0(getStockcol(), 90))
load("hektagmNov.rda")

fData(HEK293T2011)$pd.markers[fData(HEK293T2011)$pd.markers %in% c("Chromatin associated",
                                                                   "Cytosol",
                                                                   "Cytosol/Nucleus",
                                                                   "Lysosome",
                                                                   "Nucleus")] <- "unknown"


# check convergence
outliers <- mcmc_get_outliers(hektagmNov)

hekdiag <- gelman.diag(outliers[-c(4)]) # converged yay
hektagmNov_conv <- hektagmNov[-c(4)]
fData(HEK293T2011)$markers <- fData(HEK293T2011)$pd.markers

hekTagmNoveltyRes <- tagmNoveltyProcess(object = HEK293T2011, params = hektagmNov_conv)

hektagmNov_conv <- tagmMcmcProcess(hektagmNov_conv)
hektagm <- tagmPredict(object = HEK293T2011, params = hektagmNov_conv)
save(hekTagmNoveltyRes, file = "hekTagmNoveltyRes.rda")


## attach information to MSnset
fData(hektagm)$tagm.phenotype <- hekTagmNoveltyRes@pooledRes@mp@cl[rownames(fData(hektagm))]
fData(hektagm)$tagm.newcluster.prob <- hekTagmNoveltyRes@tagm.newcluster.prob[rownames(fData(hektagm))]

fData(hektagm)$tagm.newcluster.prob[is.na(fData(hektagm)$tagm.newcluster.prob)] <- 0

ptsze <- exp(fData(hektagm)$tagm.newcluster.prob) - 1
plot2D(hektagm, fcol = "tagm.phenotype", cex = ptsze)
addLegend(hektagm, fcol = "tagm.phenotype", ncol = 3, where = "topright", cex = 0.5)

plot2D(hektagm, fcol = "pd.2013")




toSubsethek <- (fData(hektagm)$tagm.newcluster.prob > 0.95)

# Create annotation 
annot <- data.frame(organelle = fData(hektagm)$tagm.phenotype[toSubsethek])
rownames(annot) <- rownames(fData(hektagm))[toSubsethek]
# Set up for heat map
ann_colors = list(
  organelle = c(getStockcol()[1:nlevels(hekTagmNoveltyRes@pooledRes@mp@cl)])
)
names(ann_colors$organelle) <- levels(annot$organelle)

psmSub <- rownames(fData(hektagm))[toSubsethek]

# Create Heatmap of PSM to big too run subset to 
pheatmap(hekTagmNoveltyRes@pooledRes@psm[psmSub,psmSub], useRaster = TRUE,
         annotation_row = annot, show_rownames = F, show_colnames = F, annotation_colors = ann_colors,
         treeheight_row = 0, treeheight_col = 0, main = "Heatmap of the posterior similarity matrix for the HEK data")

table(fData(hektagm)$tagm.phenotype)
table(fData(hektagm)$pd.2013)

heatmap(table(fData(hektagm)$tagm.phenotype, fData(hektagm)$pd.2013), Colv = NA, Rowv = NA, scale = "col")

cc1Hek <- enrichGO(gene = rownames(hektagm)[fData(hektagm)$tagm.phenotype == "Phenotype 1" & 
                                                 fData(hektagm)$tagm.newcluster.prob > 0.95],
                     OrgDb = "org.Hs.eg.db",
                     keyType = "UNIPROT",
                     ont = "CC",
                     universe = rownames(fData(hektagm))) #nucleaus

cc2Hek <- enrichGO(gene = rownames(hektagm)[fData(hektagm)$tagm.phenotype == "Phenotype 2" &
                                                 fData(hektagm)$tagm.newcluster.prob > 0.95],
                     OrgDb = "org.Hs.eg.db",
                     keyType = "UNIPROT",
                     ont = "CC",
                     universe = rownames(fData(hektagm))) #chromatin

cc3Hek <- enrichGO(gene = rownames(hektagm)[fData(hektagm)$tagm.phenotype == "Phenotype 3" & 
                                                 fData(hektagm)$tagm.newcluster.prob > 0.95],
                     OrgDb = "org.Hs.eg.db",
                     keyType = "UNIPROT",
                     ont = "CC",
                     universe = rownames(fData(hektagm)))

cc4Hek <- enrichGO(gene = rownames(hektagm)[fData(hektagm)$tagm.phenotype == "Phenotype 4" & 
                                              fData(hektagm)$tagm.newcluster.prob > 0.95],
                   OrgDb = "org.Hs.eg.db",
                   keyType = "UNIPROT",
                   ont = "CC",
                   universe = rownames(fData(hektagm)))

cc5Hek <- enrichGO(gene = rownames(hektagm)[fData(hektagm)$tagm.phenotype == "Phenotype 5" & 
                                              fData(hektagm)$tagm.newcluster.prob > 0.95],
                   OrgDb = "org.Hs.eg.db",
                   keyType = "UNIPROT",
                   ont = "CC",
                   universe = rownames(fData(hektagm)))

cc6Hek <- enrichGO(gene = rownames(hektagm)[fData(hektagm)$tagm.phenotype == "Phenotype 6" & 
                                              fData(hektagm)$tagm.newcluster.prob > 0.95],
                   OrgDb = "org.Hs.eg.db",
                   keyType = "UNIPROT",
                   ont = "CC",
                   universe = rownames(fData(hektagm)))

cc8Hek <- enrichGO(gene = rownames(hektagm)[fData(hektagm)$tagm.phenotype == "Phenotype 6" & 
                                              fData(hektagm)$tagm.newcluster.prob > 0.95],
                   OrgDb = "org.Hs.eg.db",
                   keyType = "UNIPROT",
                   ont = "CC",
                   universe = rownames(fData(hektagm)))

cc9Hek <- enrichGO(gene = rownames(hektagm)[fData(hektagm)$tagm.phenotype == "Phenotype 6" & 
                                              fData(hektagm)$tagm.newcluster.prob > 0.95],
                   OrgDb = "org.Hs.eg.db",
                   keyType = "UNIPROT",
                   ont = "CC",
                   universe = rownames(fData(hektagm)))

df1 <- data.frame(table(fData(hektagm)$tagm.phenotype)[order(table(fData(hektagm)$tagm.phenotype))])
colnames(df1) <- c("Phenotype", "NumProteins")
gg <- ggplot(df1, aes(x= Phenotype, y = NumProteins)) 
gg <- gg + geom_bar(stat = "Identity", fill = "steelblue", color = "steelblue") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                panel.background = element_blank(), axis.line = element_line(colour = "black"))
gg <- gg + coord_flip() + ylab("Number of proteins") + ggtitle("The distribution of proteins across phenotypes for Novelty TAGM") 
gg

df2 <- data.frame(table(fData(hektagm)$pd.2013)[order(table(fData(hektagm)$pd.2013))][-13])
colnames(df2) <- c("Phenotype", "NumProteins")
gg <- ggplot(df2, aes(x= Phenotype, y = NumProteins)) 
gg <- gg + geom_bar(stat = "Identity", fill = "steelblue", color = "steelblue") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                        panel.background = element_blank(), axis.line = element_line(colour = "black"))
gg <- gg + coord_flip() + ylab("Number of proteins") + ggtitle("The distribution of proteins across phenotypes for phenodisco") 
gg


## functional enrichment of phenodisco clusters

cc1Hekpd <- enrichGO(gene = rownames(hektagm)[fData(hektagm)$pd.2013 == "Phenotype 1" ],
                   OrgDb = "org.Hs.eg.db",
                   keyType = "UNIPROT",
                   ont = "CC",
                   universe = rownames(fData(hektagm)))

cc2Hekpd <- enrichGO(gene = rownames(hektagm)[fData(hektagm)$pd.2013 == "Phenotype 2" ],
                     OrgDb = "org.Hs.eg.db",
                     keyType = "UNIPROT",
                     ont = "CC",
                     universe = rownames(fData(hektagm)))

cc3Hekpd <- enrichGO(gene = rownames(hektagm)[fData(hektagm)$pd.2013 == "Phenotype 3" ],
                     OrgDb = "org.Hs.eg.db",
                     keyType = "UNIPROT",
                     ont = "CC",
                     universe = rownames(fData(hektagm)))

cc4Hekpd <- enrichGO(gene = rownames(hektagm)[fData(hektagm)$pd.2013 == "Phenotype 4"],
                     OrgDb = "org.Hs.eg.db",
                     keyType = "UNIPROT",
                     ont = "CC",
                     universe = rownames(fData(hektagm)))

cc5Hekpd <- enrichGO(gene = rownames(hektagm)[fData(hektagm)$pd.2013 == "Phenotype 5"],
                     OrgDb = "org.Hs.eg.db",
                     keyType = "UNIPROT",
                     ont = "CC",
                     universe = rownames(fData(hektagm)))

cc6Hekpd <- enrichGO(gene = rownames(hektagm)[fData(hektagm)$pd.2013 == "Phenotype 6"],
                     OrgDb = "org.Hs.eg.db",
                     keyType = "UNIPROT",
                     ont = "CC",
                     universe = rownames(fData(hektagm)))

cc7Hekpd <- enrichGO(gene = rownames(hektagm)[fData(hektagm)$pd.2013 == "Phenotype 7"],
                     OrgDb = "org.Hs.eg.db",
                     keyType = "UNIPROT",
                     ont = "CC",
                     universe = rownames(fData(hektagm)))

cc8Hekpd <- enrichGO(gene = rownames(hektagm)[fData(hektagm)$pd.2013 == "Phenotype 8"],
                     OrgDb = "org.Hs.eg.db",
                     keyType = "UNIPROT",
                     ont = "CC",
                     universe = rownames(fData(hektagm)))

df3 <- data.frame(c(8, 5), c(8,4), row.names = c("Total Phenotypes", "Functional Enriched"))
colnames(df3) <- c("Novelty TAGM", "phenodisco")
df3_long <- melt(as.matrix(df3))

gg2 <- ggplot(df3_long, aes(y = value, x = Var2, fill = Var1)) + geom_bar(stat = "Identity", position=position_dodge())
gg2 <- gg2 + coord_flip() + scale_fill_brewer(palette="Dark2") + scale_y_continuous(breaks = c(0,2,4,6,8,10))
gg2 <- gg2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black"),
             legend.title = element_blank())
gg2 <- gg2 + ylab("Number of novel Phenotypes") + xlab("Method") + ggtitle("Number of functional enriched phenotypes uncovered") 
gg2

head(fData(hektagm)$tagm.mcmc.probability)
fData(hektagm)$tagm.allocation <- "unknown"
levels(fData(hektagm)$tagm.allocation) <- levels(fData(hektagm)$tagm.phenotype)
fData(hektagm)$tagm.allocation[fData(hektagm)$tagm.newcluster.prob > 0.95] <- as.character(fData(hektagm)$tagm.phenotype[fData(hektagm)$tagm.newcluster.prob > 0.95])
fData(hektagm)$tagm.allocation[fData(hektagm)$tagm.mcmc.probability > 0.95] <- as.character(fData(hektagm)$tagm.mcmc.allocation[fData(hektagm)$tagm.mcmc.probability > 0.95])

df4 <- table(fData(hektagm)$tagm.allocation, fData(hektagm)$pd.2013)
df4 <- df4/max(df4)
pheatmap(df4, scale = "none", cluster_rows = F, cluster_cols = F, color = myCols)
