#' Check epistatic effects by kernel-based GWAS (genome-wide association studies)
#'
#' @param pheno Data frame where the first column is the line name (gid). The remaining columns should be a phenotype to test.
#' @param geno Data frame with the marker names in the first column. The second and third columns contain the chromosome and map position.
#'        Columns 4 and higher contain the marker scores for each line, coded as {-1, 0, 1} = {aa, Aa, AA}.
#' @param ZETA A list of covariance (relationship) matrix (K: \eqn{m \times m}) and its design matrix (Z: \eqn{n \times m}) of random effects.
#' Please set names of list "Z" and "K"! You can use more than one kernel matrix.
#' For example,
#'
#' ZETA = list(A = list(Z = Z.A, K = K.A), D = list(Z = Z.D, K = K.D))
#' \describe{
#' \item{Z.A, Z.D}{Design matrix (\eqn{n \times m}) for the random effects. So, in many cases, you can use the identity matrix.}
#' \item{K.A, K.D}{Different kernels which express some relationships between lines.}
#' }
#' For example, K.A is additive relationship matrix for the covariance between lines, and K.D is dominance relationship matrix.
#' @param covariate A \eqn{n \times 1} vector or a \eqn{n \times p _ 1} matrix. You can insert continuous values, such as other traits or genotype score for special markers.
#' This argument is regarded as one of the fixed effects.
#' @param covariate.factor A \eqn{n \times p _ 2} dataframe. You should assign a factor vector for each column.
#' Then RGWAS changes this argument into model matrix, and this model matrix will be included in the model as fixed effects.
#' @param structure.matrix You can use structure matrix calculated by structure analysis when there are population structure.
#' You should not use this argument with n.PC > 0.
#' @param n.PC Number of principal components to include as fixed effects. Default is 0 (equals K model).
#' @param min.MAF Specifies the minimum minor allele frequency (MAF).
#' If a marker has a MAF less than min.MAF, it is assigned a zero score.
#' @param n.core Setting n.core > 1 will enable parallel execution on a machine with multiple cores (use only at UNIX command line).
#' @param test.method RGWAS supports two methods to test effects of each SNP-set.
#' \describe{
#' \item{"LR"}{Likelihood-ratio test, relatively slow, but accurate (default).}
#' \item{"score"}{Score test, much faster than LR, but sometimes overestimate -log10(p).}
#' }
#' @param dominance.eff If this argument is TRUE, dominance effect is included in the model,
#' and additive x dominance and dominance x dominance are also tested as epistatic effects.
#' When you use inbred lines, please set this argument FALSE.
#' @param haplotype If the number of lines of your data is large (maybe > 100), you should set haplotype = TRUE.
#'             When haplotype = TRUE, haplotype-based kernel will be used for calculating -log10(p).
#'             (So the dimension of this gram matrix will be smaller.)
#'             The result won't be changed, but the time for the calculation will be shorter.
#' @param num.hap When haplotype = TRUE, you can set the number of haplotypes which you expect.
#'           Then similar arrays are considered as the same haplotype, and then make kernel(K.SNP) whose dimension is num.hap x num.hap.
#'           When num.hap = NULL (default), num.hap will be set as the maximum number which reflects the difference between lines.
#' @param window.size.half This argument decides how many SNPs (around the SNP you want to test) are used to calculated K.SNP.
#' More precisely, the number of SNPs will be 2 * window.size.half + 1.
#' @param window.slide This argument determines how often you test markers. If window.slide = 1, every marker will be tested.
#' If you want to perform SNP set by bins, please set window.slide = 2 * window.size.half + 1.
#' @param chi0.mixture RAINBOWR assumes the deviance is considered to follow a x chisq(df = 0) + (1 - a) x chisq(df = r).
#' where r is the degree of freedom.
#' The argument chi0.mixture is a (0 <= a < 1), and default is 0.5.
#' @param optimizer The function used in the optimization process. We offer "optim", "optimx", and "nlminb" functions.
#' @param gene.set If you have information of gene (or haplotype block), you can use it to perform kernel-based GWAS.
#'            You should assign your gene information to gene.set in the form of a "data.frame" (whose dimension is (the number of gene) x 2).
#'            In the first column, you should assign the gene name. And in the second column, you should assign the names of each marker,
#'            which correspond to the marker names of "geno" argument.
#' @param plot.epi.3d If TRUE, draw 3d plot
#' @param plot.epi.2d If TRUE, draw 2d plot
#' @param main.epi.3d The title of 3d plot. If this argument is NULL, trait name is set as the title.
#' @param main.epi.2d The title of 2d plot. If this argument is NULL, trait name is set as the title.
#' @param saveName When drawing any plot, you can save plots in png format. In saveName, you should substitute the name you want to save.
#' When saveName = NULL, the plot is not saved.
#' @param verbose If this argument is TRUE, messages for the current steps will be shown.
#' @param verbose2 If this argument is TRUE, welcome message will be shown.
#' @param count When count is TRUE, you can know how far RGWAS has ended with percent display.
#' @param time When time is TRUE, you can know how much time it took to perform RGWAS.
#'
#'
#' @return
#' \describe{
#' \item{$map}{Map information for SNPs which are tested epistatic effects.}
#' \item{$scores}{\describe{
#' \item{$scores}{This is the matrix which contains -log10(p) calculated by the test about epistasis effects.}
#' \item{$x, $y}{The information of the positions of SNPs detected by regular GWAS.
#'  These vectors are used when drawing plots. Each output correspond to the repliction of row and column of scores.}
#' \item{$z}{This is a vector of $scores.  This vector is also used when drawing plots.}
#' }
#' }
#' }
#'
#'
#' @references
#' Storey, J.D. and Tibshirani, R. (2003) Statistical significance for genomewide studies. Proc Natl Acad Sci. 100(16): 9440-9445.
#'
#' Yu, J. et al. (2006) A unified mixed-model method for association mapping that accounts for multiple levels of relatedness. Nat Genet. 38(2): 203-208.
#'
#' Kang, H.M. et al. (2008) Efficient Control of Population Structure in Model Organism Association Mapping. Genetics. 178(3): 1709-1723.
#'
#' Endelman, J.B. (2011) Ridge Regression and Other Kernels for Genomic Selection with R Package rrBLUP. Plant Genome J. 4(3): 250.
#'
#' Endelman, J.B. and Jannink, J.L. (2012) Shrinkage Estimation of the Realized Relationship Matrix. G3 Genes, Genomes, Genet. 2(11): 1405-1413.
#'
#' Su, G. et al. (2012) Estimating Additive and Non-Additive Genetic Variances and Predicting Genetic Merits Using Genome-Wide Dense Single Nucleotide Polymorphism Markers. PLoS One. 7(9): 1-7.
#'
#' Zhou, X. and Stephens, M. (2012) Genome-wide efficient mixed-model analysis for association studies. Nat Genet. 44(7): 821-824.
#'
#' Listgarten, J. et al. (2013) A powerful and efficient set test for genetic markers that handles confounders. Bioinformatics. 29(12): 1526-1533.
#'
#' Lippert, C. et al. (2014) Greater power and computational efficiency for kernel-based association testing of sets of genetic variants. Bioinformatics. 30(22): 3206-3214.
#'
#' Jiang, Y. and Reif, J.C. (2015) Modeling epistasis in genomic selection. Genetics. 201(2): 759-768.
#'
#'
#' @example R/examples/RGWAS.epistasis_example.R
#'
#'
#'
#'
RGWAS.epistasis <- function(pheno, geno, ZETA = NULL, covariate = NULL, covariate.factor = NULL,
                            structure.matrix = NULL, n.PC = 0, min.MAF = 0.02, n.core = 1,
                            test.method = "LR", dominance.eff = TRUE, haplotype = TRUE, num.hap = NULL,
                            window.size.half = 5, window.slide = 1, chi0.mixture = 0.5, optimizer = "nlminb",
                            gene.set = NULL, plot.epi.3d = TRUE, plot.epi.2d = TRUE,
                            main.epi.3d = NULL, main.epi.2d = NULL, saveName = NULL, verbose = TRUE,
                            verbose2 = FALSE, count = TRUE, time = TRUE){

  #### The start of the RGWAS function ####
  start <- Sys.time()


  #### Some settings to perform RGWAS ####
  if(verbose2){
    welcome_to_RGWAS()
  }

  ### For phenotype ###
  n.sample.pheno <- nrow(pheno)
  n.pheno <- ncol(pheno) - 1
  pheno.ix <- 2:ncol(pheno)
  pheno.names <- colnames(pheno)[2:ncol(pheno)]
  lines.name.pheno <- as.character(pheno[, 1])

  ### For covariate ###
  X0 <- matrix(1, n.sample.pheno, 1)
  colnames(X0) <- "Intercept"
  rownames(X0) <- lines.name.pheno

  if(!is.null(covariate)){
    p1 <- ncol(covariate)
    X0 <- cbind(X0, scale(covariate))
  }

  if(!is.null(covariate.factor)){
    p2 <- ncol(covariate.factor)
    for (i in 1:p2) {
      cov.fac.now <- covariate.factor[, i]
      if (length(unique(cov.fac.now)) > 1) {
        model.mat.now <- model.matrix(~ x - 1, data.frame(x = cov.fac.now))
        colnames(model.mat.now) <- paste0("cov.fac.", i, "_", 1:length(unique(cov.fac.now)))
        X0 <- cbind(X0, model.mat.now[, -length(unique(cov.fac.now))])
      }
    }
  }

  if(!is.null(structure.matrix)){
    colnames(structure.matrix) <- paste0("subpop", 1:ncol(structure.matrix))
    X0 <- cbind(X0, structure.matrix)

    n.PC <- 0
  }


  ### For genotype ###
  geno <- geno[order(geno[, 2], geno[, 3]), ]
  lines.name.geno <- colnames(geno)[-c(1:3)]
  M0 <- t(geno[, -c(1:3)])
  map <- geno[, 1:3]
  marker <- as.character(map[, 1])
  chr <- map[, 2]
  chr.tab <- table(chr)
  chr.max <- max(chr)
  chr.cum <- cumsum(chr.tab)
  pos <- map[, 3]
  cum.pos <- pos
  if(length(chr.tab) != 1){
    for(i in 1:(chr.max - 1)){
      cum.pos[(chr.cum[i] + 1): (chr.cum[i + 1])] <- pos[(chr.cum[i] + 1):(chr.cum[i + 1])] + cum.pos[chr.cum[i]]
    }
  }
  n.mark <- ncol(M0)
  rownames(M0) <- lines.name.geno


  ### Match phenotype and genotype ###
  pheno.mat <- as.matrix(pheno[, -1, drop = FALSE])
  rownames(pheno.mat) <- lines.name.pheno

  modification.res <- modify.data(pheno.mat = pheno.mat, geno.mat = M0,
                                  pheno.labels = NULL, geno.names = NULL,
                                  map = NULL, return.ZETA = TRUE, return.GWAS.format = FALSE)

  pheno.mat.modi <- modification.res$pheno.modi
  match.modi <- match(rownames(pheno.mat.modi), pheno[, 1])
  pheno.match <- pheno[match.modi, ]

  M <- modification.res$geno.modi

  n.line <- nrow(M)
  X <- as.matrix(X0[match.modi, ])


  if(is.null(ZETA)){
    ZETA <- modification.res$ZETA
  }else{
    ZETA.check <- any(unlist(lapply(ZETA, function(x) {
      (is.null(rownames(x$Z))) | (is.null(colnames(x$Z)))
    })))

    if(ZETA.check){
      stop("No row names or column names for design matrix Z!!
           Please fill them with row : line (variety) names for phenotypes.
           and column : line (variety) names for genotypes.")
    }

    ZETA <- lapply(ZETA, function(x){
      Z.match.pheno.no <- match(rownames(pheno.mat.modi), rownames(x$Z))
      Z.match.geno.no <- match(rownames(M), rownames(x$Z))

      Z.modi <- x$Z[Z.match.pheno.no, Z.match.geno.no]
      K.modi <- x$K[Z.match.geno.no, Z.match.geno.no]

      return(list(Z = Z.modi, K = K.modi))
    })
    }
  K.A <- ZETA[[1]]$K
  Z.A <- ZETA[[1]]$Z


  ### For covariates (again) ###
  if(n.PC > 0){
    eigen.K.A <- eigen(K.A)
    eig.K.vec <- eigen.K.A$vectors

    PC.part <- Z.A %*% eig.K.vec[, 1:n.PC]
    colnames(PC.part) <- paste0("n.PC_", 1:n.PC)

    X <- cbind(X, PC.part)
  }
  X <- make.full(X)


  ### Some settings ###
  trait.names <- colnames(pheno)[pheno.ix]
  if(is.null(gene.set)){
    n.scores.each <- (chr.tab + (window.slide - 1)) %/% window.slide
    n.scores <- sum(n.scores.each)
  }else{
    n.scores <- length(unique(gene.set[, 1]))
  }

  if(n.pheno != 1){
    all.scores.epi <- rep(list(NA), n.pheno)
  }else{
    all.scores.epi <- NULL
  }

  if (n.pheno == 0) {
    stop("No phenotypes.")
  }



  ##### START RGWAS for each phenotype #####
  for(pheno.no in 1:n.pheno){
    trait.name <- trait.names[pheno.no]
    if (verbose) {
      print(paste("GWAS for trait:", trait.name))
    }
    y0 <- pheno.match[, pheno.ix[pheno.no]]
    not.NA <- which(!is.na(y0))
    y <- y0[not.NA]

    n <- length(y)

    X.now <- as.matrix(X[not.NA, ])
    ZETA.now <- lapply(ZETA, function(x) list(Z = x$Z[not.NA, ], K = x$K))

    if (sum(is.na(M)) == 0) {
      M.now <- Z.A[not.NA, ] %*% M
    } else {
      M.now <- M[apply(Z.A[not.NA, ], 1, function(x) which(x == 1)), ]
    }


    #### Calculate LL for the null hypothesis at first ####
    spI <- diag(n)
    S <- spI - tcrossprod(X.now %*% solve(crossprod(X.now)), X.now)

    if(length(ZETA) > 1){
      EMM.res0 <- EM3.cpp(y = y, X0 = X.now, ZETA = ZETA.now,
                          n.thres = 450, REML = TRUE, pred = FALSE)

      weights <- EMM.res0$weights
    }else{
      EMM.res0 <- EMM.cpp(y = y, X = X.now, ZETA = ZETA.now,
                           n.thres = 450, REML = TRUE)
      weights <- 1
    }

    ZKZt.list <- NULL
    ZKZt <- matrix(0, nrow = n, ncol = n)
    for(ZKZt.no in 1:length(ZETA)){
      Z.now <- ZETA.now[[ZKZt.no]]$Z
      K.now <- ZETA.now[[ZKZt.no]]$K
      ZKZt.now <- tcrossprod(Z.now %*% K.now, Z.now)
      ZKZt.weighted <- ZKZt.now * weights[ZKZt.no]

      ZKZt.list <- c(ZKZt.list, list(ZKZt.weighted))
      ZKZt <- ZKZt + ZKZt.weighted
    }

    if(test.method == "LR"){
      LL0 <- EMM.res0$LL

      spectral.res <- spectralG.cpp(ZETA = ZETA.now, X = X.now, weights = weights,
                                    return.G = TRUE, return.SGS = TRUE, spectral.method = "eigen")
      eigen.G <- spectral.res[[1]]
      eigen.SGS <- spectral.res[[2]]
    }else{
      if(test.method == "score"){
        LL0 <- EMM.res0$LL
        Vu <- EMM.res0$Vu
        Ve <- EMM.res0$Ve

        Gu <- tcrossprod(ZKZt)
        Ge <- diag(n)
        V0 <- Vu * Gu + Ve * Ge


        P0 <- MASS::ginv(S %*% V0 %*% S)
      }else{
        stop("We only support 'LR' (likelihood-ratio test) and 'score' (score test)!")
      }
    }





    #### Calculating the value of -log10(p) for each SNPs ####
    if ((n.core > 1) & requireNamespace("parallel", quietly = TRUE)) {
      warning("Sorry. n.core > 1 options have not been implemented yet. We will use n.core = 1 instead.")
      if(test.method == "LR"){
        scores.epi <- score.calc.epistasis.LR(M.now = M.now, y = y, X.now = X.now, ZETA.now = ZETA.now,
                                              eigen.SGS = eigen.SGS, eigen.G = eigen.G, map = map,
                                              haplotype = haplotype, num.hap = num.hap, window.size.half = window.size.half,
                                              window.slide = window.slide, chi0.mixture = chi0.mixture,
                                              gene.set = gene.set,  dominance.eff = dominance.eff, optimizer = optimizer,
                                              min.MAF = min.MAF, count = count)
      }else{
        scores.epi <- score.calc.epistasis.score(M.now = M.now, ZETA.now = ZETA.now, y = y,
                                                 Gu = Gu, Ge = Ge, P0 = P0, map = map, haplotype = haplotype,
                                                 num.hap = num.hap, window.size.half = window.size.half,
                                                 window.slide = window.slide, chi0.mixture = chi0.mixture,
                                                 gene.set = gene.set, dominance.eff = dominance.eff,
                                                 min.MAF = min.MAF, count = count)
      }
    }else {
      if(test.method == "LR"){
        scores.epi <- score.calc.epistasis.LR(M.now = M.now, y = y, X.now = X.now, ZETA.now = ZETA.now,
                                              eigen.SGS = eigen.SGS, eigen.G = eigen.G, map = map,
                                              haplotype = haplotype, num.hap = num.hap, window.size.half = window.size.half,
                                              window.slide = window.slide, chi0.mixture = chi0.mixture,
                                              gene.set = gene.set,  dominance.eff = dominance.eff, optimizer = optimizer,
                                              min.MAF = min.MAF, count = count)
      }else{
        scores.epi <- score.calc.epistasis.score(M.now = M.now, ZETA.now = ZETA.now, y = y,
                                                 Gu = Gu, Ge = Ge, P0 = P0, map = map, haplotype = haplotype,
                                                 num.hap = num.hap, window.size.half = window.size.half,
                                                 window.slide = window.slide, chi0.mixture = chi0.mixture,
                                                 gene.set = gene.set, dominance.eff = dominance.eff,
                                                 min.MAF = min.MAF, count = count)
      }
    }

    if(is.null(gene.set)){
      window.centers <- as.numeric(rownames(scores.epi))
      map2 <- map[window.centers, ]
      cum.pos2 <- cum.pos[window.centers]
    }else{
      map20 <- genesetmap(map = map, gene.set = gene.set, cumulative = TRUE)
      map2 <- map20[, 1:3]
      cum.pos.set.mean <- c(map20[, 4])
      cum.pos2 <- cum.pos.set.mean
    }


    x.3d <- rep(cum.pos2, n.scores)
    y.3d <- rep(cum.pos2, each = n.scores)
    z.3d <- c(scores.epi)

    epi.res <- list(scores = scores.epi, x = x.3d, y = y.3d, z = z.3d)

    if(n.pheno >= 2){
      all.scores.epi[[pheno.no]] <- epi.res
    }else{
      all.scores.epi <- epi.res
    }

    if(is.null(main.epi.3d)){
      main.epi.3d <- trait.name
    }
    if(is.null(main.epi.2d)){
      main.epi.2d <- trait.name
    }

    if (verbose) {
      print("Now Plotting (3d plot for epistasis). Please Wait.")
    }
    manhattan3(input = epi.res, cum.pos = cum.pos, plot.epi.3d = plot.epi.3d,
               plot.epi.2d = plot.epi.2d,  main.epi.3d = main.epi.3d,
               main.epi.2d = main.epi.2d, saveName = saveName)

  }



  end <- Sys.time()

  if(time){
    print(end - start)
  }

  return(list(map = map2, scores = all.scores.epi))
}
