\dontshow{
  ### Import RAINBOWR
  require(RAINBOWR)


  ### Load example datasets
  data("Rice_Zhao_etal")
  Rice_geno_score <- Rice_Zhao_etal$genoScore
  Rice_geno_map <- Rice_Zhao_etal$genoMap
  Rice_pheno <- Rice_Zhao_etal$pheno


  ### Select one trait for example
  trait.name <- "Flowering.time.at.Arkansas"
  y <- as.matrix(Rice_pheno[1:30, trait.name, drop = FALSE])
  # use first 30 accessions


  ### Remove SNPs whose MAF <= 0.05
  x.0 <- t(Rice_geno_score[, 1:30])
  MAF.cut.res <- MAF.cut(x.0 = x.0, map.0 = Rice_geno_map)
  x <- MAF.cut.res$x
  map <- MAF.cut.res$map


  ### Assume adjacent individuals are regarded as "neighbors"
  xAdj <- array(data = NA, dim = dim(x), dimnames = dimnames(x))
  for (i in 1:nrow(x)) {
    adjs <- (i - 1):(i + 1)
    adjs <- adjs[adjs %in% 1:nrow(x)]
    adjs <- adjs[adjs != i]
    nAdjs <- length(adjs)

    xAdj[i, ] <- x[i, , drop = FALSE] *
      apply(X = x[adjs, , drop = FALSE],
            MARGIN = 2, FUN = mean)
  }


  ### Estimate additive genomic relationship matrix (GRM)
  K.A <- tcrossprod(x) / ncol(x)
  K.Adj <- tcrossprod(xAdj) / ncol(xAdj)  # for neighbor kernel


  ### Modify data
  Z <- design.Z(pheno.labels = rownames(y),
                geno.names = rownames(K.A))  ### design matrix for random effects
  pheno.mat <- y[rownames(Z), , drop = FALSE]
  ZETA <- list(A = list(Z = Z, K = K.A),
               Adj = list(Z = Z, K = K.Adj))


  ### Prepare for covairance structure between two random effects
  K12 <- tcrossprod(x, xAdj) / sqrt(ncol(x) * ncol(xAdj))
  K21 <- tcrossprod(xAdj, x) / sqrt(ncol(x) * ncol(xAdj))

  covList <- rep(list(rep(list(NULL), 2)), 2)
  covList[[1]][[2]] <- K12
  covList[[2]][[1]] <- K21


  ### Solve multi-kernel linear mixed effects model (2 random efects)
  ### conidering covariance structure
  EM3cov.res <- EM3.cov(y = pheno.mat,
                        X0 = NULL,
                        ZETA = ZETA,
                        covList = covList)
  Vu <- EM3cov.res$Vu   ### estimated genetic variance
  Ve <- EM3cov.res$Ve   ### estimated residual variance
  weights <- EM3cov.res$weights   ### estimated proportion of two genetic variances
  herit <- Vu * weights / (Vu + Ve)   ### genomic heritability (additive, neighbor)
  rho <- EM3cov.res$rhosMat[2, 1]   ### correlation parameter

  beta <- EM3cov.res$beta   ### Here, this is an intercept.
  u.each <- EM3cov.res$u.each   ### estimated genotypic values (additive, neighbor)
}


\donttest{
  ### Import RAINBOWR
  require(RAINBOWR)

  ### Load example datasets
  data("Rice_Zhao_etal")
  Rice_geno_score <- Rice_Zhao_etal$genoScore
  Rice_geno_map <- Rice_Zhao_etal$genoMap
  Rice_pheno <- Rice_Zhao_etal$pheno

  ### View each dataset
  See(Rice_geno_score)
  See(Rice_geno_map)
  See(Rice_pheno)

  ### Select one trait for example
  trait.name <- "Flowering.time.at.Arkansas"
  y <- as.matrix(Rice_pheno[, trait.name, drop = FALSE])

  ### Remove SNPs whose MAF <= 0.05
  x.0 <- t(Rice_geno_score)
  MAF.cut.res <- MAF.cut(x.0 = x.0, map.0 = Rice_geno_map)
  x <- MAF.cut.res$x
  map <- MAF.cut.res$map


  ### Assume adjacent individuals are regarded as "neighbors"
  xAdj <- array(data = NA, dim = dim(x), dimnames = dimnames(x))
  for (i in 1:nrow(x)) {
    adjs <- (i - 1):(i + 1)
    adjs <- adjs[adjs %in% 1:nrow(x)]
    adjs <- adjs[adjs != i]
    nAdjs <- length(adjs)

    xAdj[i, ] <- x[i, , drop = FALSE] *
      apply(X = x[adjs, , drop = FALSE],
            MARGIN = 2, FUN = mean)
  }


  ### Estimate additive genomic relationship matrix (GRM)
  K.A <- tcrossprod(x) / ncol(x)
  K.Adj <- tcrossprod(xAdj) / ncol(xAdj)  # for neighbor kernel


  ### Modify data
  Z <- design.Z(pheno.labels = rownames(y),
                geno.names = rownames(K.A))  ### design matrix for random effects
  pheno.mat <- y[rownames(Z), , drop = FALSE]
  ZETA <- list(A = list(Z = Z, K = K.A),
               Adj = list(Z = Z, K = K.Adj))


  ### Prepare for covairance structure between two random effects
  K12 <- tcrossprod(x, xAdj) / sqrt(ncol(x) * ncol(xAdj))
  K21 <- tcrossprod(xAdj, x) / sqrt(ncol(x) * ncol(xAdj))

  covList <- rep(list(rep(list(NULL), 2)), 2)
  covList[[1]][[2]] <- K12
  covList[[2]][[1]] <- K21


  ### Solve multi-kernel linear mixed effects model (2 random efects)
  ### conidering covariance structure
  EM3cov.res <- EM3.cov(y = pheno.mat,
                        X0 = NULL,
                        ZETA = ZETA,
                        covList = covList)
  (Vu <- EM3cov.res$Vu)   ### estimated genetic variance
  (Ve <- EM3cov.res$Ve)   ### estimated residual variance
  (weights <- EM3cov.res$weights)   ### estimated proportion of two genetic variances
  (herit <- Vu * weights / (Vu + Ve))   ### genomic heritability (additive, neighbor)
  (rho <- EM3cov.res$rhosMat[2, 1])   ### correlation parameter

  (beta <- EM3cov.res$beta)   ### Here, this is an intercept.
  u.each <- EM3cov.res$u.each   ### estimated genotypic values (additive, neighbor)
  See(u.each)


  ### Perform genomic prediction with 10-fold cross validation (multi-kernel)
  noNA <- !is.na(c(pheno.mat))   ### NA (missing) in the phenotype data

  phenoNoNA <- pheno.mat[noNA, , drop = FALSE]   ### remove NA
  ZETANoNA <- ZETA
  ZETANoNA <- lapply(X = ZETANoNA, FUN = function (List) {
    List$Z <- List$Z[noNA, ]

    return(List)
  })   ### remove NA


  nFold <- 10    ### # of folds
  nLine <- nrow(phenoNoNA)
  idCV <- sample(1:nLine %% nFold)   ### assign random ids for cross-validation
  idCV[idCV == 0] <- nFold

  yPred <- yPredCov <- rep(NA, nLine)


  for (noCV in 1:nFold) {
    print(paste0("Fold: ", noCV))
    yTrain <- phenoNoNA
    yTrain[idCV == noCV, ] <- NA   ### prepare test data

    EM3.resCV <- EM3.cpp(y = yTrain, X0 = NULL,
                         ZETA = ZETANoNA)   ### prediction
    EM3cov.resCV <- EM3.cov(y = yTrain, X0 = NULL,
                            ZETA = ZETANoNA,
                            covList = covList)   ### prediction
    yTest <- EM3.resCV$y.pred     ### predicted values
    yTestCov <- EM3cov.resCV$y.pred

    yPred[idCV == noCV] <- yTest[idCV == noCV]
    yPredCov[idCV == noCV] <- yTestCov[idCV == noCV]
  }

  ### Plot the results
  plotRange <- range(phenoNoNA, yPred)
  plot(x = phenoNoNA, y = yPred,
       xlim = plotRange, ylim = plotRange,
       xlab = "Observed values", ylab = "Predicted values (EM3.cpp)",
       main = "Results of Genomic Prediction (multi-kernel)",
       cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.3)
  abline(a = 0, b = 1, col = 2, lwd = 2, lty = 2)
  R2 <- cor(x = phenoNoNA[, 1], y = yPred) ^ 2
  text(x = plotRange[2] - 10,
       y = plotRange[1] + 10,
       paste0("R2 = ", round(R2, 3)),
       cex = 1.5)


  plotRange <- range(phenoNoNA, yPred)
  plot(x = phenoNoNA, y = yPredCov,
       xlim = plotRange, ylim = plotRange,
       xlab = "Observed values", ylab = "Predicted values (EM3.cov)",
       main = "Results of Genomic Prediction (multi-kernel with covariance)",
       cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.3)
  abline(a = 0, b = 1, col = 2, lwd = 2, lty = 2)
  R2 <- cor(x = phenoNoNA[, 1], y = yPredCov) ^ 2
  text(x = plotRange[2] - 10,
       y = plotRange[1] + 10,
       paste0("R2 = ", round(R2, 3)),
       cex = 1.5)


  plotRange <- range(yPred, yPredCov)
  plot(x = yPred, y = yPredCov,
       xlim = plotRange, ylim = plotRange,
       xlab = "Predicted values (EM3.cpp)", ylab = "Predicted values (EM3.cov)",
       main = "Comparison of Multi-Kernel Genomic Prediction",
       cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.3)
  abline(a = 0, b = 1, col = 2, lwd = 2, lty = 2)
  R2 <- cor(x = yPred, y = yPredCov) ^ 2
  text(x = plotRange[2] - 10,
       y = plotRange[1] + 10,
       paste0("R2 = ", round(R2, 3)),
       cex = 1.5)
}
