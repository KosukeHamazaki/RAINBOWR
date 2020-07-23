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
  x.0 <- t(Rice_geno_score)
  MAF.cut.res <- MAF.cut(x.0 = x.0, map.0 = Rice_geno_map)
  x <- MAF.cut.res$x
  map <- MAF.cut.res$map
  
  
  ### Estimate additive genomic relationship matrix (GRM) & epistatic relationship matrix
  K.A <- calcGRM(genoMat = x)
  K.AA <- K.A * K.A   ### additive x additive epistatic effects
  
  
  ### Modify data
  Z <- design.Z(pheno.labels = rownames(y),
                geno.names = rownames(K.A))  ### design matrix for random effects
  pheno.mat <- y[rownames(Z), , drop = FALSE]
  ZETA <- list(A = list(Z = Z, K = K.A),
               AA = list(Z = Z, K = K.AA))
  
  
  ### Solve multi-kernel linear mixed effects model (2 random efects)
  EM3.res <- EM3.cpp(y = pheno.mat, X = NULL, ZETA = ZETA)
  Vu <- EM3.res$Vu   ### estimated genetic variance
  Ve <- EM3.res$Ve   ### estimated residual variance
  weights <- EM3.res$weights   ### estimated proportion of two genetic variances
  herit <- Vu * weights / (Vu + Ve)   ### genomic heritability (additive, additive x additive)
  
  beta <- EM3.res$beta   ### Here, this is an intercept.
  u <- EM3.res$u   ### estimated genotypic values (additive, additive x additive)
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
  
  
  ### Estimate additive genomic relationship matrix (GRM) & epistatic relationship matrix
  K.A <- calcGRM(genoMat = x) 
  K.AA <- K.A * K.A   ### additive x additive epistatic effects
  
  
  ### Modify data
  Z <- design.Z(pheno.labels = rownames(y),
                geno.names = rownames(K.A))  ### design matrix for random effects
  pheno.mat <- y[rownames(Z), , drop = FALSE]
  ZETA <- list(A = list(Z = Z, K = K.A),
               AA = list(Z = Z, K = K.AA))
  
  
  ### Solve multi-kernel linear mixed effects model (2 random efects)
  EM3.res <- EM3.cpp(y = pheno.mat, X = NULL, ZETA = ZETA)
  (Vu <- EM3.res$Vu)   ### estimated genetic variance
  (Ve <- EM3.res$Ve)   ### estimated residual variance
  (weights <- EM3.res$weights)   ### estimated proportion of two genetic variances
  (herit <- Vu * weights / (Vu + Ve))   ### genomic heritability (additive, additive x additive)
  
  (beta <- EM3.res$beta)   ### Here, this is an intercept.
  u <- EM3.res$u   ### estimated genotypic values (additive, additive x additive)
  See(u)
  
  
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
  
  yPred <- rep(NA, nLine)
  
  for (noCV in 1:nFold) {
    print(paste0("Fold: ", noCV))
    yTrain <- phenoNoNA
    yTrain[idCV == noCV, ] <- NA   ### prepare test data
    
    EM3.resCV <- EM3.cpp(y = yTrain, X = NULL, ZETA = ZETANoNA)   ### prediction
    yTest <-  EM3.resCV$y.pred     ### predicted values
    
    yPred[idCV == noCV] <- yTest[idCV == noCV]
  }
  
  ### Plot the results
  plotRange <- range(phenoNoNA, yPred)
  plot(x = phenoNoNA, y = yPred,xlim = plotRange, ylim = plotRange,
       xlab = "Observed values", ylab = "Predicted values",
       main = "Results of Genomic Prediction (multi-kernel)",
       cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.3)
  abline(a = 0, b = 1, col = 2, lwd = 2, lty = 2)
  R2 <- cor(x = phenoNoNA[, 1], y = yPred) ^ 2
  text(x = plotRange[2] - 10,
       y = plotRange[1] + 10,
       paste0("R2 = ", round(R2, 3)), 
       cex = 1.5)
}
