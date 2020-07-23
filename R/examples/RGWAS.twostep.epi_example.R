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
  y <- Rice_pheno[1:50, trait.name, drop = FALSE]
  # use first 50 accessions
  
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
  
  
  
  ### Perform two-step epistasis GWAS (single-snp GWAS -> Check epistasis for significant markers)
  twostep.epi.res <- RGWAS.twostep.epi(pheno = pheno.GWAS, geno = geno.GWAS, ZETA = ZETA,
                                       n.PC = 4, test.method = "LR", gene.set = NULL,
                                       window.size.half = 10, window.slide = 21,
                                       plot.Manhattan.1 = FALSE, plot.qq.1 = FALSE,
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
  y <- Rice_pheno[, trait.name, drop = FALSE]
  
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
  
  
  
  
  ### Perform two-step epistasis GWAS (single-snp GWAS -> Check epistasis for significant markers)
  twostep.epi.res <- RGWAS.twostep.epi(pheno = pheno.GWAS, geno = geno.GWAS, ZETA = ZETA,
                                       n.PC = 4, test.method = "LR", gene.set = NULL,
                                       window.size.half = 10, window.slide = 21)
  
  See(twostep.epi.res$epistasis$scores)
}
