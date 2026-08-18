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

  ### Assume adjacent individuals are regarded as "neighbors"
  xAdj <- array(
    data = NA,
    dim = dim(x),
    dimnames = dimnames(x)
  )

  for (i in seq_len(nrow(x))) {
    adjs <- (i - 1):(i + 1)
    adjs <- adjs[adjs %in% seq_len(nrow(x))]
    adjs <- adjs[adjs != i]

    xAdj[i, ] <- x[i, , drop = FALSE] *
      apply(
        X = x[adjs, , drop = FALSE],
        MARGIN = 2,
        FUN = mean
      )
  }

  ### Estimate additive genomic relationship matrix (GRM)
  ### and neighbor relationship matrix
  K.A <- tcrossprod(x) / ncol(x)
  K.Adj <- tcrossprod(xAdj) / ncol(xAdj)

  ### Modify data
  Z <- design.Z(
    pheno.labels = rownames(y),
    geno.names = rownames(K.A)
  )
  pheno.mat <- y[rownames(Z), , drop = FALSE]

  ZETA <- list(
    A = list(Z = Z, K = K.A),
    Adj = list(Z = Z, K = K.Adj)
  )

  ### Prepare covariance structures between random effects
  K12 <- tcrossprod(x, xAdj) /
    sqrt(ncol(x) * ncol(xAdj))

  K21 <- tcrossprod(xAdj, x) /
    sqrt(ncol(x) * ncol(xAdj))

  covList <- rep(
    list(rep(list(NULL), 2)),
    2
  )

  covList[[1]][[2]] <- K12
  covList[[2]][[1]] <- K21

  ### Solve multi-kernel linear mixed-effects model
  ### considering covariance structures
  EM3cov.res <- EM3.cov(
    y = pheno.mat,
    X0 = NULL,
    ZETA = ZETA,
    covList = covList,
    REML = TRUE
  )

  ### Calculate approximate standard errors
  EM3cov.SE.res <- calcParamSE.cov(
    y = pheno.mat,
    EM3.cov.res = EM3cov.res,
    ZETA = ZETA,
    covList = covList,
    X0 = NULL,
    REML = TRUE
  )

  EM3cov.SE.res$parameterTable
}
