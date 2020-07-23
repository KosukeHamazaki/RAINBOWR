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
  modify.data.res <- modify.data(pheno.mat = y, geno.mat = x, map = map,
                                 return.ZETA = TRUE, return.GWAS.format = TRUE)
  pheno.GWAS <- modify.data.res$pheno.GWAS
  geno.GWAS <- modify.data.res$geno.GWAS
  ZETA <- modify.data.res$ZETA
  
  
  ### Perform two step SNP-set GWAS (single-snp GWAS -> SNP-set GWAS for significant markers)
  twostep.SNP_set.res <- RGWAS.twostep(pheno = pheno.GWAS, geno = geno.GWAS, ZETA = ZETA,
                                       kernel.percent = 0.2, n.PC = 4, test.method.2 = "LR",
                                       kernel.method = "linear", gene.set = NULL,
                                       test.effect.2 = "additive", window.size.half = 3,
                                       window.slide = 2, plot.qq.1 = FALSE,
                                       plot.Manhattan.1 = FALSE, plot.qq.2 = FALSE,
                                       plot.Manhattan.2 = FALSE, verbose = FALSE,
                                       count = FALSE, time = FALSE)
  
  ### Column 4 contains -log10(p) values for markers with the first method (single-SNP GWAS)
  ### Column 5 contains -log10(p) values for markers with the second method (SNP-set GWAS)
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
  
  
  ### Perform two step SNP-set GWAS (single-snp GWAS -> SNP-set GWAS for significant markers)
  twostep.SNP_set.res <- RGWAS.twostep(pheno = pheno.GWAS, geno = geno.GWAS, ZETA = ZETA,
                                       kernel.percent = 0.2, n.PC = 4, test.method.2 = "LR",
                                       kernel.method = "linear", gene.set = NULL,
                                       test.effect.2 = "additive", window.size.half = 3,
                                       window.slide = 2)
  
  See(twostep.SNP_set.res$D)
  ### Column 4 contains -log10(p) values for markers with the first method (single-SNP GWAS)
  ### Column 5 contains -log10(p) values for markers with the second method (SNP-set GWAS)
}
