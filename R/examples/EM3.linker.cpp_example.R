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
  
  
  ### Estimate additive genomic relationship matrix (GRM)
  K.A <- calcGRM(genoMat = x)
  
  
  ### Modify data
  Z <- design.Z(pheno.labels = rownames(y),
                geno.names = rownames(K.A))  ### design matrix for random effects
  pheno.mat <- y[rownames(Z), , drop = FALSE]
  ZETA <- list(A = list(Z = Z, K = K.A))
  
  
  ### Including the first 20 SNPs
  W.A <- x[, 1:20]    ### marker genotype data of first 20 SNPs
  
  Zs0 <- list(A.part = Z)
  Ws0 <- list(A.part = W.A)       ### This will be regarded as linear kernel
  ### for the variance-covariance matrix of another random effects.
  
  
  ### Solve multi-kernel linear mixed effects model (2 random efects)
  EM3.linker.res <- EM3.linker.cpp(y0 = pheno.mat, X0 = NULL, ZETA = ZETA,
                                   Zs0 = Zs0, Ws0 = Ws0)
  Vu <- EM3.linker.res$Vu   ### estimated genetic variance
  Ve <- EM3.linker.res$Ve   ### estimated residual variance
  weights <- EM3.linker.res$weights   ### estimated proportion of two genetic variances
  herit <- Vu * weights / (Vu + Ve)   ### genomic heritability (all chromosomes, chromosome 12)
  
  beta <- EM3.linker.res$beta   ### Here, this is an intercept.
  u <- EM3.linker.res$u   ### estimated genotypic values (all chromosomes, chromosome 12)
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
  
  
  ### Estimate additive genomic relationship matrix (GRM)
  K.A <- calcGRM(genoMat = x)
  
  
  ### Modify data
  Z <- design.Z(pheno.labels = rownames(y),
                geno.names = rownames(K.A))  ### design matrix for random effects
  pheno.mat <- y[rownames(Z), , drop = FALSE]
  ZETA <- list(A = list(Z = Z, K = K.A))
  
  
  ### Including the additional linear kernel for chromosome 12
  chrNo <- 12
  W.A <- x[, map$chr == chrNo]    ### marker genotype data of chromosome 12
  
  Zs0 <- list(A.part = Z)
  Ws0 <- list(A.part = W.A)       ### This will be regarded as linear kernel
  ### for the variance-covariance matrix of another random effects.
  
  
  ### Solve multi-kernel linear mixed effects model (2 random efects)
  EM3.linker.res <- EM3.linker.cpp(y0 = pheno.mat, X0 = NULL, ZETA = ZETA,
                                   Zs0 = Zs0, Ws0 = Ws0)
  (Vu <- EM3.linker.res$Vu)   ### estimated genetic variance
  (Ve <- EM3.linker.res$Ve)   ### estimated residual variance
  (weights <- EM3.linker.res$weights)   ### estimated proportion of two genetic variances
  (herit <- Vu * weights / (Vu + Ve))   ### genomic heritability (all chromosomes, chromosome 12)
  
  (beta <- EM3.linker.res$beta)   ### Here, this is an intercept.
  u <- EM3.linker.res$u   ### estimated genotypic values (all chromosomes, chromosome 12)
  See(u)
}
