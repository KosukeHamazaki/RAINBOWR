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
  # use first 30 acessions
  
  
  ### Remove SNPs whose MAF <= 0.05
  x.0 <- t(Rice_geno_score)
  MAF.cut.res <- MAF.cut(x.0 = x.0, map.0 = Rice_geno_map)
  x <- MAF.cut.res$x
  map <- MAF.cut.res$map
  
  
  ### Estimate genomic relationship matrix (GRM)
  K.A <- calcGRM(genoMat = x) 
  
  
  ### Modify data
  modify.data.res <- modify.data(pheno.mat = y, geno.mat = x, map = map,
                                 return.ZETA = TRUE, return.GWAS.format = TRUE)
  pheno.GWAS <- modify.data.res$pheno.GWAS
  geno.GWAS <- modify.data.res$geno.GWAS
  ZETA <- modify.data.res$ZETA
  
  
  ### Check epistatic effects (by regarding 161 SNPs as one SNP-set)
  epistasis.res <- RGWAS.epistasis(pheno = pheno.GWAS, geno = geno.GWAS, ZETA = ZETA,
                                   n.PC = 4, test.method = "score", gene.set = NULL,
                                   window.size.half = 40, window.slide = 81,
                                   plot.epi.3d = FALSE, plot.epi.2d = FALSE,
                                   verbose = FALSE, count = FALSE, time = FALSE)
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


  ### Estimate genomic relationship matrix (GRM)
  K.A <- calcGRM(genoMat = x) 


  ### Modify data
  modify.data.res <- modify.data(pheno.mat = y, geno.mat = x, map = map,
                                 return.ZETA = TRUE, return.GWAS.format = TRUE)
  pheno.GWAS <- modify.data.res$pheno.GWAS
  geno.GWAS <- modify.data.res$geno.GWAS
  ZETA <- modify.data.res$ZETA


  ### View each data for RAINBOWR
  See(pheno.GWAS)
  See(geno.GWAS)
  str(ZETA)


  ### Check epistatic effects (by regarding 11 SNPs as one SNP-set)
  epistasis.res <- RGWAS.epistasis(pheno = pheno.GWAS, geno = geno.GWAS, ZETA = ZETA,
                                   n.PC = 4, test.method = "score", gene.set = NULL,
                                   window.size.half = 5, window.slide = 11)

  See(epistasis.res$scores$scores)
}
