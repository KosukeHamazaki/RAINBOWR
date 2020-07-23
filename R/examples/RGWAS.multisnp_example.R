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
  geno.GWAS <- modify.data.res$geno.GWAS[1:300, ]  ### first 300 SNPs
  ZETA <- modify.data.res$ZETA
  
  
  
  ### Perform SNP-set GWAS (by regarding 41 SNPs as one SNP-set)
  SNP_set.res <- RGWAS.multisnp(pheno = pheno.GWAS, geno = geno.GWAS,
                                ZETA = ZETA, n.PC = 4, test.method = "LR",
                                kernel.method = "linear", gene.set = NULL,
                                test.effect = "additive", window.size.half = 10,
                                window.slide = 21, plot.Manhattan = FALSE,
                                plot.qq = FALSE, verbose = FALSE,
                                count = FALSE, time = FALSE)
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
  
  
  ### Perform SNP-set GWAS (by regarding 21 SNPs as one SNP-set)
  SNP_set.res <- RGWAS.multisnp(pheno = pheno.GWAS, geno = geno.GWAS,
                                ZETA = ZETA, n.PC = 4, test.method = "LR",
                                kernel.method = "linear", gene.set = NULL,
                                test.effect = "additive", window.size.half = 10,
                                window.slide = 21)
  See(SNP_set.res$D)  ### Column 4 contains -log10(p) values for markers
  
  ### Perform SNP-set GWAS 2 (by regarding 11 SNPs as one SNP-set with sliding window)
  ### It will take almost 25 minutes...
  SNP_set.res2 <- RGWAS.multisnp(pheno = pheno.GWAS, geno = geno.GWAS,
                                 ZETA = ZETA, n.PC = 4, test.method = "LR",
                                 kernel.method = "linear", gene.set = NULL,
                                 test.effect = "additive", window.size.half = 5,
                                 window.slide = 1)
  See(SNP_set.res2$D)  ### Column 4 contains -log10(p) values for markers
}
