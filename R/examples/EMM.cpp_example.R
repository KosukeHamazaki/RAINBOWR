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
  
  ### Estimate genomic relationship matrix (GRM)
  K.A <- calcGRM(genoMat = x)
  
  ### Modify data
  modify.res <- modify.data(pheno.mat = y, geno.mat = x, return.ZETA = TRUE)
  pheno.mat <- modify.res$pheno.modi
  ZETA <- modify.res$ZETA
  
  
  ### Solve linear mixed effects model
  EMM.res <- EMM.cpp(y = pheno.mat, X = NULL, ZETA = ZETA)
  Vu <- EMM.res$Vu   ### estimated genetic variance
  Ve <- EMM.res$Ve   ### estimated residual variance
  herit <- Vu / (Vu + Ve)   ### genomic heritability
  
  beta <- EMM.res$beta   ### Here, this is an intercept.
  u <- EMM.res$u   ### estimated genotypic values
}

### Perform genomic prediction with 10-fold cross validation
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
  
  ### Estimate genomic relationship matrix (GRM)
  K.A <- calcGRM(genoMat = x)
  
  ### Modify data
  modify.res <- modify.data(pheno.mat = y, geno.mat = x, return.ZETA = TRUE)
  pheno.mat <- modify.res$pheno.modi
  ZETA <- modify.res$ZETA
  
  
  ### Solve linear mixed effects model
  EMM.res <- EMM.cpp(y = pheno.mat, X = NULL, ZETA = ZETA)
  (Vu <- EMM.res$Vu)   ### estimated genetic variance
  (Ve <- EMM.res$Ve)   ### estimated residual variance
  (herit <- Vu / (Vu + Ve))   ### genomic heritability
  
  (beta <- EMM.res$beta)   ### Here, this is an intercept.
  u <- EMM.res$u   ### estimated genotypic values
  See(u)
  
  noNA <- !is.na(c(pheno.mat))   ### NA (missing) in the phenotype data
  
  phenoNoNA <- pheno.mat[noNA, , drop = FALSE]   ### remove NA
  ZETANoNA <- ZETA
  ZETANoNA$A$Z <- ZETA$A$Z[noNA, ]   ### remove NA
  
  
  nFold <- 10    ### # of folds
  nLine <- nrow(phenoNoNA)
  idCV <- sample(1:nLine %% nFold)   ### assign random ids for cross-validation
  idCV[idCV == 0] <- nFold
  
  yPred <- rep(NA, nLine)
  
  for (noCV in 1:nFold) {
    yTrain <- phenoNoNA
    yTrain[idCV == noCV, ] <- NA   ### prepare test data
    
    EMM.resCV <- EMM.cpp(y = yTrain, X = NULL, ZETA = ZETANoNA)   ### prediction
    yTest <-  EMM.resCV$beta + EMM.resCV$u   ### predicted values
    
    yPred[idCV == noCV] <- (yTest[noNA])[idCV == noCV]
  }
  
  ### Plot the results
  plotRange <- range(phenoNoNA, yPred)
  plot(x = phenoNoNA, y = yPred,xlim = plotRange, ylim = plotRange,
       xlab = "Observed values", ylab = "Predicted values",
       main = "Results of Genomic Prediction",
       cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.3)
  abline(a = 0, b = 1, col = 2, lwd = 2, lty = 2)
  R2 <- cor(x = phenoNoNA[, 1], y = yPred) ^ 2
  text(x = plotRange[2] - 10,
       y = plotRange[1] + 10,
       paste0("R2 = ", round(R2, 3)), 
       cex = 1.5)
}
