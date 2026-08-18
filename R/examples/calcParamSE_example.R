\dontrun{
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
  ### and additive x additive epistatic relationship matrix
  K.A <- calcGRM(genoMat = x)
  K.AA <- K.A * K.A

  ### Modify data
  Z <- design.Z(
    pheno.labels = rownames(y),
    geno.names = rownames(K.A)
  )
  pheno.mat <- y[rownames(Z), , drop = FALSE]


  ### Single-kernel model
  ZETA.A <- list(
    A = list(Z = Z, K = K.A)
  )

  EMM.res <- EMM.cpp(
    y = pheno.mat,
    X = NULL,
    ZETA = ZETA.A,
    REML = TRUE
  )

  EMM.SE.res <- calcParamSE(
    y = pheno.mat,
    MM.res = EMM.res,
    ZETA = ZETA.A,
    X0 = NULL,
    REML = TRUE
  )

  EMM.SE.res$parameterTable


  ### Multi-kernel model
  ZETA <- list(
    A = list(Z = Z, K = K.A),
    AA = list(Z = Z, K = K.AA)
  )

  EM3.res <- EM3.cpp(
    y = pheno.mat,
    X0 = NULL,
    ZETA = ZETA,
    REML = TRUE
  )

  EM3.SE.res <- calcParamSE(
    y = pheno.mat,
    MM.res = EM3.res,
    ZETA = ZETA,
    X0 = NULL,
    REML = TRUE
  )

  EM3.SE.res$parameterTable
}
