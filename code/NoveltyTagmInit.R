object = tan2009r1
fcol = "markers"
method = "MCMC"
numIter = 2000
burnin = 100L
thin = 5L
mu0 = NULL
lambda0 = 0.01
nu0 = NULL
S0 = NULL
beta0 = NULL
u = 2
v = 10
K_bar = 5

setStockcol(paste0(getStockcol(), 90))

tanTagmparams <- tagmMcmcTrain_Nov(object = tan2009r1, BPPARAM = SerialParam())

## Check orginal pipeline holds
tanTagmparams <- tagmMcmcProcess(tanTagmparams)
tan2009r1 <- tagmPredict(object = tan2009r1, params = tanTagmparams, probJoint = TRUE)

ptsze <- fData(tan2009r1)$tagm.mcmc.probability

plot2D(tan2009r1, fcol = "tagm.mcmc.allocation", cex = ptsze)

## Novelty detection methods
tanChain1 <- pRoloc:::chains(tanTagmparams)[[1]]

tagmNoveltyChainProcess <- function(object,
                               params,
                               fcol = "markers") {
  ## Checks that params is object of class MCMCChain
  stopifnot(inherits(params, "MCMCChain"))
  
    ## Make marker allocation matrix and align rownames
  markersubset <- markerMSnSet(object = object, fcol)
  markers <- getMarkerClasses(markersubset, fcol = fcol)
  markerallocationMat <- matrix(NA, ncol = ncol(params@Component), nrow = nrow(markersubset))
  rownames(markerallocationMat) <- rownames(markersubset)
  
  K <- params@K
  K_markers <- nlevels(fData(markersubset)[, fcol])
  
  
  ## populate marker allocation matrix
  for (j in seq.int(K_markers)) {
    toSubset <- rownames(markersubset)[fData(markersubset)[, fcol] == markers[j]]
    markerallocationMat[toSubset, ] <- rep(j, ncol(params@Component))
  }
  
  # Form component matrix with markers
  combinedComponent <- rbind(params@Component, markerallocationMat)
  
  
  #Compute PSM and summarise using MCclust methods
  modal.k <- which.max(tabulate(apply(params@Component, 2, function(x) length(unique(x)))))
  psm <- mcclust::comp.psm(t(combinedComponent))
  mp <- mcclust::maxpear(psm = psm, max.k = K)
  
  # Match rownames 
  names(mp$cl) <- rownames(combinedComponent)
  rownames(psm) <- colnames(psm) <- c(rownames(unknownMSnSet(object)), rownames(markersubset))
  
  # Finding mapping from cluster to organelle because of label switching
  mapping <- matrix(NA, ncol = K_markers, nrow = 1)
  for (j in seq.int(K_markers)) {
    toSubset <- rownames(combinedComponent
                         [rownames(markersubset), ])[combinedComponent[rownames(markersubset), 1] == j]
    mapping[j] <- mp$cl[toSubset][1]
  }
  mapping <- c(mapping, seq.int(1:K)[!seq.int(1:K) %in% mapping])
  
  # Match up clusters with organelles and phenotypes
  for (j in seq.int(K)) {
    mp$cl[mp$cl == mapping[j]] <- rownames(params@ComponentParam@mk)[j]
  }
  
  # Make factor
  mp$cl <- factor(mp$cl)
 
  ## Make output
  .out <- list(psm = psm,
              mp = mp)
  
  return(.out) 
}

tagmNoveltyProcess <- function(object,
                               params,
                               fcol = "markers") {
  ## Checks that params is object of class MCMCChain
  stopifnot(inherits(params, "MCMCParams"))
  numChains <- length(params)
  
  ## Storage
  psms <- vector("list", length = numChains)
  mps <- vector("list", length = numChains)
  
  ## Marker Classes
  markersubset <- markerMSnSet(object = object, fcol)
  K_markers <- nlevels(fData(markersubset)[, fcol])
  
  for (j in seq_len(numChains)){
    ## For each Chain compute psm and maxpear
    .res <- tagmNoveltyChainProcess(object = object,
                            params = pRoloc:::chains(params)[[j]],
                            fcol = fcol)
    psms[[j]] <- .res$psm
    mps[[j]] <- .res$mp
  }
  
  ## Pool Chains and get psm and maxpear for combined chains
  pooledparams <- mcmc_pool_chains(params)
  .pooledres <- tagmNoveltyChainProcess(object = object,
                                        params = pRoloc:::chains(pooledparams)[[1]],
                                        fcol = fcol)
  psmCombined <- .pooledres$psm
  mpCombined <- .pooledres$mp
  
  pooledparams <- tagmMcmcProcess(pooledparams)
  
  tagm.newcluster.prob <- 1 - rowSums(pooledparams@summary@tagm.joint[, 1:K_markers])
  
  # Make output
  .out <- list(psmCombined = psmCombined,
               mpCombined = mpCombined,
               psms = psms, mps = mps,
               tagm.newcluster.prob = tagm.newcluster.prob)
  
  return(.out)
}  

noveltyRes2 <- tagmNoveltyProcess(object = tan2009r1, params = tanTagmparams)

# Create annotation 
annot <- data.frame(organelle = noveltyRes2$mps[[3]]$cl)
# Set up for heat map
ann_colors = list(
  organelle = c(getStockcol()[1:nlevels(noveltyRes2$mps[[3]]$cl)])
)
names(ann_colors$organelle) <- levels(annot$organelle)

# Create Heatmap of PSM
pheatmap(noveltyRes2$psms[[3]], useRaster = TRUE, annotation_row = annot, show_rownames = F, show_colnames = F, annotation_colors = ann_colors,
         treeheight_row = 0, treeheight_col = 0)

## attach information to MSnset
fData(tan2009r1)$tagm.phenotype <- noveltyRes2$mpCombined$cl[rownames(fData(tan2009r1))]
fData(tan2009r1)$tagm.newcluster.prob <- noveltyRes2$tagm.newcluster.prob[rownames(fData(tan2009r1))]

ptsze <- exp(fData(tan2009r1)$tagm.newcluster.prob) - 1
plot2D(tan2009r1, fcol = "tagm.phenotype", cex = ptsze)
addLegend(tan2009r1, fcol = "tagm.phenotype", ncol = 2, where = "topleft")

cc1 <- enrichGO(gene = rownames(tan2009r1)[fData(tan2009r1)$tagm.phenotype == "Phenotype 1"],
               OrgDb = "org.Dm.eg.db",
               keyType = "UNIPROT",
               ont = "CC",
               universe = rownames(fData(tan2009r1)))

cc2 <- enrichGO(gene = rownames(tan2009r1)[fData(tan2009r1)$tagm.phenotype == "Phenotype 2"],
                OrgDb = "org.Dm.eg.db",
                keyType = "UNIPROT",
                ont = "MF",
                universe = rownames(fData(tan2009r1)))



## Trye dunkley removing ER and Ribosome and vacuole annotations
fData(dunkley2006)$markers[fData(dunkley2006)$markers %in% c("ER lumen", "Ribosome", "TGN")] <- "unknown"


# Dunkley seems a little weird
dunkleyTagmParams <- tagmMcmcTrain_Nov(object = dunkley2006, BPPARAM = SerialParam(), K_bar = 4)
dunkleyNoveltyRes <- tagmNoveltyProcess(object = dunkley2006, params = dunkleyTagmParams)

# Create annotation 
annot <- data.frame(organelle = dunkleyNoveltyRes$mpCombined$cl)
# Set up for heat map
ann_colors = list(
  organelle = c(getStockcol()[1:nlevels(dunkleyNoveltyRes$mpCombined$cl)])
)
names(ann_colors$organelle) <- levels(annot$organelle)

# Create Heatmap of PSM
pheatmap(dunkleyNoveltyRes$psmCombined, useRaster = TRUE, annotation_row = annot, show_rownames = F, show_colnames = F, annotation_colors = ann_colors,
         treeheight_row = 0, treeheight_col = 0)

## attach information to MSnset
fData(dunkley2006)$tagm.phenotype <- dunkleyNoveltyRes$mpCombined$cl[rownames(fData(dunkley2006))]
fData(dunkley2006)$tagm.newcluster.prob <- dunkleyNoveltyRes$tagm.newcluster.prob[rownames(fData(dunkley2006))]

ptsze <- exp(fData(dunkley2006)$tagm.newcluster.prob) - 1
plot2D(dunkley2006, fcol = "tagm.phenotype", cex = ptsze)
addLegend(dunkley2006, fcol = "tagm.phenotype", ncol = 2, where = "topleft", cex = 0.7)


# Working try and remove annotations
fData(HEK293T2011)$markers[fData(HEK293T2011)$markers %in% c("Chromatin associated",
                                                             "Cytosol", "Cytosol/Nucleus", "Nucleus")] <- "unknown"


hekTagmParams <- tagmMcmcTrain_Nov(object = HEK293T2011, BPPARAM = SerialParam(), K_bar = 4)
hekNoveltyRes <- tagmNoveltyProcess(object = HEK293T2011, params = hekTagmParams)


# Create annotation 
annot <- data.frame(organelle = hekNoveltyRes$mpCombined$cl)
# Set up for heat map
ann_colors = list(
  organelle = c(getStockcol()[1:nlevels(hekNoveltyRes$mpCombined$cl)])
)
names(ann_colors$organelle) <- levels(annot$organelle)

# Create Heatmap of PSM
pheatmap(hekNoveltyRes$psmCombined, useRaster = TRUE, annotation_row = annot, show_rownames = F, show_colnames = F, annotation_colors = ann_colors,
         treeheight_row = 0, treeheight_col = 0)

## attach information to MSnset
fData(HEK293T2011)$tagm.phenotype <- hekNoveltyRes$mpCombined$cl[rownames(fData(HEK293T2011))]
fData(HEK293T2011)$tagm.newcluster.prob <- hekNoveltyRes$tagm.newcluster.prob[rownames(fData(HEK293T2011))]

ptsze <- exp(fData(HEK293T2011)$tagm.newcluster.prob) - 1
plot2D(HEK293T2011, fcol = "tagm.phenotype", cex = ptsze)
addLegend(HEK293T2011, fcol = "tagm.phenotype", ncol = 2, where = "topleft", cex = 0.7)
