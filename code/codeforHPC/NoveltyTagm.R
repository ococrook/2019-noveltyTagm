tagmMcmcPredict <- function(object,
                            params,
                            fcol = "markers",
                            probJoint = FALSE,
                            probOutlier = TRUE) {
  stopifnot(inherits(params, "MCMCParams"))
  ## Checks for object and MCMCParams match
  stopifnot(featureNames(unknownMSnSet(object, fcol = fcol))
            == rownames(params@summary@posteriorEstimates))
  
  ## Create marker set and size
  markerSet <- markerMSnSet(object, fcol = fcol)
  markers <- getMarkerClasses(object, fcol = fcol)
  M <- nrow(markerSet)
  K <- chains(params)[[1]]@K
  
  
  ## Get Summary object from MCMCParams maybe better to check
  ## columns exist/pass which objects we need
  .tagm.allocation <- c(as.character(params@summary@posteriorEstimates[,"tagm.allocation"]),
                        as.character(fData(markerSet)[, fcol]))
  .tagm.probability <- c(params@summary@posteriorEstimates[,"tagm.probability"],
                         rep(1, M)) ## set all probabilities of markers to 1.
  .tagm.probability.lowerquantile <- c(params@summary@posteriorEstimates[,"tagm.probability.lowerquantile"],
                                       rep(1, M)) ## set all probabilities of markers to 1.
  .tagm.probability.upperquantile <- c(params@summary@posteriorEstimates[,"tagm.probability.upperquantile"],
                                       rep(1, M)) ## set all probabilities of markers to 1.
  .tagm.mean.shannon <- c(params@summary@posteriorEstimates[,"tagm.mean.shannon"],
                          rep(0, M)) ## set all probabilities of markers to 0
  
  ## Create data frame to store new summaries
  .tagm.summary <- data.frame(tagm.mcmc.allocation = .tagm.allocation ,
                              tagm.mcmc.probability = .tagm.probability,
                              tagm.mcmc.probability.lowerquantile = .tagm.probability.lowerquantile,
                              tagm.mcmc.probability.upperquantile = .tagm.probability.upperquantile,
                              tagm.mcmc.mean.shannon = .tagm.mean.shannon)
  if (probOutlier)
    .tagm.summary$tagm.mcmc.outlier <-
    c(params@summary@posteriorEstimates[, "tagm.probability.Outlier"],
      rep(0, M)) ## set all probabilities of markers to 0
  
  
  ## Check number of rows match and add feature names
  stopifnot(nrow(.tagm.summary) == nrow(object))
  rownames(.tagm.summary) <- c(rownames(params@summary@posteriorEstimates),
                               rownames(markerSet))
  
  ## Append data to fData of MSnSet
  fData(object) <- cbind(fData(object), .tagm.summary[rownames(fData(object)),])
  
  if  (probJoint) {
    ## create allocation matrix for markers
    .probmat <- matrix(0, nrow = nrow(markerSet), ncol = K)
    .class <- fData(markerSet)[, fcol]
    for (j in seq_len(nrow(markerSet))) {
      ## give markers prob 1
      .probmat[j, as.numeric(factor(.class), seq(1, length(unique(.class))))[j]] <- 1
    }
    colnames(.probmat) <- markers
    rownames(.probmat) <- rownames(markerSet)
    .joint <- rbind(params@summary@tagm.joint, .probmat)
    fData(object)$tagm.mcmc.joint <- .joint[rownames(fData(object)), ]
  }
  
  return(object)
}


tagmMcmcTrain_Nov <- function(object,
                          fcol = "markers",
                          method = "MCMC",
                          numIter = 1000L,
                          burnin = 100L,
                          thin = 5L,
                          mu0 = NULL,
                          lambda0 = 0.01,
                          nu0 = NULL,
                          S0 = NULL,
                          beta0 = NULL,
                          u = 2,
                          v = 10,
                          numChains = 4L,
                          K_bar = 5,
                          BPPARAM = BiocParallel::bpparam()) {
  
  ## get expression marker data
  markersubset <- markerMSnSet(object, fcol = fcol)
  markers <- getMarkerClasses(markersubset, fcol = fcol)
  mydata <- exprs(markersubset)
  X <- exprs(unknownMSnSet(object, fcol = fcol))
  
  ## get data dize
  N <- nrow(mydata)
  D <- ncol(mydata)
  K <- as.integer(length(markers) + K_bar)
  
  ## set priors
  if (is.null(nu0))
    nu0 <- D + 2
  
  if (is.null(S0))
    S0 <- diag( colSums(( mydata - mean( mydata)) ^ 2) / N)/( K ^ (1/D))
  
  if (is.null(mu0))
    mu0 <- colMeans( mydata)
  
  if (is.null(beta0))
    beta0 <- rep(1, K)
  
  ## Store Priors
  .priors <- list(mu0 = mu0,
                  lambda0 = lambda0,
                  nu0 = nu0,
                  S0 = S0,
                  beta0 = beta0)
  
  ## chains run in parallel, repeating number of iterations
  .res <- bplapply(rep(numIter, numChains),
                   FUN = tagmMcmcChain_Nov,
                   object = object,
                   fcol = fcol,
                   method = "MCMC",
                   burnin = burnin,
                   thin = thin,
                   mu0 = mu0,
                   lambda0 = lambda0,
                   nu0 = nu0,
                   S0 = S0,
                   beta0 = beta0,
                   u = u,
                   v = v,
                   K_bar = K_bar,
                   BPPARAM = BPPARAM)
  
  ## Construct class MCMCChains
  .ans <- pRoloc:::.MCMCChains(chains = .res)
  
  ## Construct class MCMCParams
  out <- pRoloc:::.MCMCParams(method = "TAGM.MCMC",
                     chains = .ans,
                     priors = .priors,
                     summary = pRoloc:::.MCMCSummary())
  
  return(out)
}





tagmMcmcChain_Nov <- function(object,
                          fcol = "markers",
                          method = "MCMC",
                          numIter = 1000L,
                          burnin = 100L,
                          thin = 5L,
                          mu0 = NULL,
                          lambda0 = 0.01,
                          nu0 = NULL,
                          S0 = NULL,
                          beta0 = NULL,
                          u = 2,
                          v = 10,
                          K_bar = 5) {
  if (burnin >= numIter)
    stop("Burnin is larger than iterations you will not retain any samples")
  
  ## on the fly number of samples to be retained
  retained <- seq.int(burnin + 1L, numIter , thin)
  
  ## get expression marker data
  markersubset <- markerMSnSet(object, fcol = fcol)
  markers <- getMarkerClasses(markersubset, fcol = fcol)
  mydata <- exprs(markersubset)
  X <- exprs(unknownMSnSet(object, fcol = fcol))
  
  ## get data dize
  N <- nrow(mydata)
  D <- ncol(mydata)
  K_markers <- length(markers)
  K <- as.integer(length(markers) + K_bar)
  
  ## set empirical priors
  if (is.null(nu0))
    nu0 <- D + 2
  
  if (is.null(S0))
    S0 <- diag( colSums(( mydata - mean( mydata)) ^ 2) / N)/( K ^ (1/D))
  
  if (is.null(mu0))
    mu0 <- colMeans( mydata)
  
  if (is.null(beta0))
    beta0 <- rep(1, K)
  
  ## save priors
  .priors <- list(mu0 = mu0,
                  lambda0 = lambda0,
                  nu0 = nu0,
                  S0 = S0,
                  beta0 = beta0)
  
  ## create storage for posterior parameters
  mk <- matrix(0, nrow = K, ncol = D)
  lambdak <- matrix(0, nrow = K, ncol = 1)
  nuk <- matrix(0, nrow = K, ncol = 1)
  sk <- array(0, dim = c(K, D, D))
  
  ## create storage for cluster parameters
  xk <- matrix(0, nrow = K, ncol = D)
  
  ## update prior with training data
  nk <- tabulate(fData(markersubset)[, fcol])
  nk <- c(nk, rep(0, K_bar)) # attach pheno-discovery clusters
  
  for (j in seq.int(K_markers))
    xk[j, ] <- colSums(mydata[fData(markersubset)[, fcol] == markers[j], ])/nk[j]
  
  lambdak <- lambda0 + nk
  nuk <- nu0 + nk
  mk <- t((t(nk * xk) + lambda0 * mu0)) / lambdak
  betak <- beta0 + nk
  
  for(j in seq.int(K)) {
    if (j > K_markers) { ## pheno-discovery clusters get prior
      sk[j, , ] <- S0
    } else {
    sk[j, , ] <- S0 + t(mydata[fData(markersubset)[, fcol] == markers[j], ]) %*%
    mydata[fData(markersubset)[, fcol] == markers[j],] +
    lambda0 * mu0 %*% t(mu0) - lambdak[j] * mk[j, ] %*% t(mk[j, ])
    }
  }
  
  ## global parameters
  M <- colMeans(exprs(object))
  V <- cov(exprs(object))/2
  
  ## storage
  Component <- matrix(0, nrow = nrow(X), ncol = numIter)
  ComponentProb <- array(0, c(nrow(X), numIter, K))
  Outlier <- matrix(0, nrow = nrow(X), ncol = numIter)
  OutlierProb <- array(0, c(nrow(X), numIter, 2))
  
  ## initially assigned all unlabelled points to clusters greedily
  for(j in seq.int(K))
    ComponentProb[, 1, j] <- pRoloc:::dmvtCpp(X,
                                     mu_ = mk[j, ],
                                     sigma_ = (1 + lambdak[j]) * sk[j, , ] / (lambdak[j] * (nuk[j] - D + 1)),
                                     df_ = nuk[j] - D + 1,
                                     log_ = TRUE,
                                     ncores_ = 1,
                                     isChol_ = FALSE)
  
  Component[, 1] <- apply(X = ComponentProb[, 1, ], 1, FUN = which.max)
  
  ## randomly reassign half the proteins to empty clusters
  prm <- sample(x = seq.int(1:nrow(X)), size = floor(nrow(X)/2), replace = FALSE)
  Component[prm, 1] <- sample(x = seq.int(1:K_bar), size = floor(nrow(X)/2), replace = TRUE) + K_markers
  
  ## Need statistics for pheno clusters
  nk <- nk + tabulate(Component[prm, 1])
  
  for (j in seq.int(K_markers + 1, K))
    xk[j, ] <- colSums(X[prm,][Component[prm, 1] == j,])/nk[j]
  
  lambdak <- lambda0 + nk
  nuk <- nu0 + nk
  mk <- t((t(nk * xk) + lambda0 * mu0)) / lambdak
  betak <- beta0 + nk
  
  for (j in seq.int(K_markers + 1, K))
    sk[j, , ] <- S0 + t(X[prm,][Component[prm, 1] == j,]) %*% X[prm,][Component[prm, 1] == j,] +
    lambda0 * mu0 %*% t(mu0) - lambdak[j] * mk[j, ] %*% t(mk[j, ])
  
  
  ## initially assign all proteins to global component unless they've been assigned to new phenotype
  Outlier[seq.int(nrow(X))[-prm], 1] <- 0
  Outlier[prm, 1] <- 1
  
  ## initial allocation statistics
  tau1 <- sum(Outlier[, 1] == 1) + N
  tau2 <- sum(Outlier[, 1] == 0)
  
  for (t in seq.int(2L, numIter)) {
    
    if (t%%500 == 0){
      print(t)
    }
    
    ## consider each protein in turn
    for (i in seq.int(nrow(X))) {
      
      ## if assigned to a cluster remove statistics
      if ( Outlier[i, t - 1] == 1) {
        idx <- Component[i, t - 1] ## temporary variable index
        tempS <- mk[idx, ] %*% t(mk[idx, ]) ## temporary scatter, since mk to be update on next line
        mk[idx, ] <- (lambdak[idx] * mk[idx, ] - X[i, ]) / (lambdak[idx] - 1)
        lambdak[idx] <- lambdak[idx] - 1
        nuk[idx] <- nuk[idx] - 1
        nk[idx] <- nk[idx] - 1
        tau1 <- tau1 - 1
        sk[idx, , ] <- sk[idx, , ] - (X[i, ] %*% t(X[i, ])) +
          (lambdak[idx] + 1) * tempS  - lambdak[idx] * mk[idx,] %*% t(mk[idx,])
      } else {
        if (t > 2) {
          tau2 <- tau2 - 1
        }
      }
      
      ## compute probability of belonging to each organelle
      ## precompute terms for speed
      weight <- (nk + betak)/(sum(nk) + sum(betak) - 1) ## Component weights
      sigmak <- ((1 + lambdak) * sk)/(lambdak * (nuk - D + 1)) ## scale matrix
      degf <- nuk - D + 1 ## degrees freedom
      for (j in seq.int(K)) {
        ComponentProb[i, t, j] <- log(weight[j]) + pRoloc:::dmvtCpp(X[i, ,drop = FALSE],
                                                           mu_ = mk[j, ],
                                                           sigma_ = sigmak[j, , ],
                                                           df_ = degf[j],
                                                           log_ = TRUE,
                                                           ncores_ = 1,
                                                           isChol_ = FALSE)
      }
      
      ## normalise with underflow correction
      c <-  max(ComponentProb[i, t, ])
      ComponentProb[i, t , ] <- exp(ComponentProb[i ,t , ] - c) / sum(exp(ComponentProb[i, t, ] - c))
      
      ## sample component
      Component[i, t] <- sample(x = 1:K, size = 1, prob = ComponentProb[i, t , ] )
      
      ## compute outlier allocation
      n <- nrow(object)
      idk <- Component[i, t] ## temporary allocation variable
      OutlierProb[i, t, 1] <- log((tau1 + v)/(n + u + v - 1)) + pRoloc:::dmvtCpp(X[i, ,drop=FALSE],
                                                                        mu_ = mk[idk, ],
                                                                        sigma_ = sigmak[idk,,],
                                                                        df_ = degf[idk],
                                                                        log_ = TRUE,
                                                                        ncores_ = 1,
                                                                        isChol_ = FALSE)
      OutlierProb[i, t, 2] <- log((tau2 + u)/(n + u + v - 1)) + pRoloc:::dmvtCpp(X[i, ,drop = FALSE],
                                                                        mu_ = M,
                                                                        sigma_ = V,
                                                                        df_ = 4,
                                                                        log_ = TRUE,
                                                                        ncores_ = 1,
                                                                        isChol_ = FALSE)
      ## normalise and sample
      OutlierProb[i, t, ] <- exp(OutlierProb[i, t, ])/sum(exp(OutlierProb[i, t, ]))
      Outlier[i, t] <- sample(x = c(1, 0), 1, prob = OutlierProb[i, t, ]) ## reversed sample so 2nd entry is prob of 0.
      
      ## allocation statistics
      if ( Outlier[i, t] == 1) {
        idx <- Component[i, t] ## temporary variable index
        tempS <- mk[idx, ] %*% t(mk[idx, ]) ## temporary scatter, since mk to be update on next line
        mk[idx, ] <- (lambdak[idx] * mk[idx, ] + X[i, ]) / (lambdak[idx] + 1)
        lambdak[idx] <- lambdak[idx] + 1
        nuk[idx] <- nuk[idx] + 1
        nk[idx] <- nk[idx] + 1
        tau1 <- tau1 + 1
        #if ( t == 2) {
        #  tau2 <- tau2 - 1 ## default on first round is phi = 0
        #}
        sk[idx, , ] <- sk[idx,,] + (X[i,] %*% t(X[i,])) +
          (lambdak[idx] - 1) * tempS  - lambdak[idx] * mk[idx,] %*% t(mk[idx,])
      } else {
        #if (t > 2) {
          tau2 <- tau2 + 1
        #}
      }
    } ## end loop over proteins
  } ## end iterations
  
  ## create names for objects
  rownames(mk) <- names(lambdak) <- names(nuk) <- dimnames(sk)[[1]] <- c(markers, paste0("Phenotype ", seq.int(K_bar)))
  colnames(mk) <- dimnames(sk)[[2]] <- dimnames(sk)[[3]] <- sampleNames(object)
  
  ## name storage
  p_names <- rownames(X)
  rownames(Component) <- rownames(ComponentProb) <- rownames(Outlier) <- rownames(OutlierProb) <- p_names
  dimnames(ComponentProb)[[3]] <- c(markers, paste0("Phenotype ", seq.int(K_bar)))
  
  ## save Component parameters
  .ComponentParam <- pRoloc:::.ComponentParam(K = K, D = D,
                                     mk = mk,
                                     lambdak = lambdak,
                                     nuk = nuk,
                                     sk = sk)
  ## apply thinning and burn-in
  .Component <- Component[, retained]
  .ComponentProb <- ComponentProb[, retained, ]
  .Outlier <- Outlier[, retained]
  .OutlierProb <- OutlierProb[, retained, ]
  
  ## make MCMCChain object
  .MCMCChain <- pRoloc:::.MCMCChain(n = length(retained),
                           K = K,
                           N = nrow(X),
                           Component = .Component,
                           ComponentProb = .ComponentProb,
                           Outlier = .Outlier,
                           OutlierProb = .OutlierProb,
                           ComponentParam = .ComponentParam)
  
  return(.MCMCChain)
}
