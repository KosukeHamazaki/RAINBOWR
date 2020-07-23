#' Testing multiple SNPs simultaneously for GWAS
#'
#' @description This function performs SNP-set GWAS (genome-wide association studies),
#' which tests multiple SNPs (single nucleotide polymorphisms) simultaneously. The model of SNP-set GWAS is
#'
#' \deqn{y = X \beta + Q v +  Z _ {c} u _ {c} +  Z _ {r} u _ {r} + \epsilon,}
#'
#' where \eqn{y} is the vector of phenotypic values,
#' \eqn{X \beta} and \eqn{Q v} are the terms of fixed effects,
#' \eqn{Z _ {c} u _ {c}} and \eqn{Z _ {c} u _ {c}} are the term of random effects and \eqn{e} is the vector of residuals.
#' \eqn{X \beta} indicates all of the fixed effects other than population structure, and often this term also plays
#' a role as an intercept. \eqn{Q v} is the term to correct the effect of population structure.
#' \eqn{Z _ {c} u _ {c}} is the term of polygenetic effects, and suppose that \eqn{u _ {c}}
#' follows the multivariate normal distribution whose variance-covariance
#' matrix is the genetic covariance matrix. \eqn{u _ {c} \sim MVN (0, K _ {c} \sigma_{c}^{2})}.
#' \eqn{Z _ {r} u _ {r}} is the term of effects for SNP-set of interest, and suppose that \eqn{u _ {r}}
#' follows the multivariate normal distribution whose variance-covariance
#' matrix is the Gram matrix (linear, exponential, or gaussian kernel)
#' calculated from marker genotype which belong to that SNP-set.
#' Therefore, \eqn{u _ {r} \sim MVN (0, K _ {r} \sigma_{r}^{2})}.
#' Finally, the residual term is assumed to identically and independently follow
#' a normal distribution as shown in the following equation.
#' \eqn{e \sim MVN (0, I \sigma_{e}^{2})}.
#'
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
#' @param test.method RGWAS supports two methods to test effects of each SNP-set.
#' \describe{
#' \item{"LR"}{Likelihood-ratio test, relatively slow, but accurate (default).}
#' \item{"score"}{Score test, much faster than LR, but sometimes overestimate -log10(p).}
#' }
#' @param n.core Setting n.core > 1 will enable parallel execution on a machine with multiple cores (use only at UNIX command line).
#' @param kernel.method It determines how to calculate kernel. There are three methods.
#' \describe{
#' \item{"gaussian"}{It is the default method. Gaussian kernel is calculated by distance matrix.}
#' \item{"exponential"}{When this method is selected, exponential kernel is calculated by distance matrix.}
#' \item{"linear"}{When this method is selected, linear kernel is calculated by NOIA methods for additive GRM.}
#'}
#' So local genomic relation matrix is regarded as kernel.
#' @param kernel.h The hyper parameter for gaussian or exponential kernel.
#' If kernel.h = "tuned", this hyper parameter is calculated as the median of off-diagonals of distance matrix of genotype data.
#' @param haplotype If the number of lines of your data is large (maybe > 100), you should set haplotype = TRUE.
#'             When haplotype = TRUE, haplotype-based kernel will be used for calculating -log10(p).
#'             (So the dimension of this gram matrix will be smaller.)
#'             The result won't be changed, but the time for the calculation will be shorter.
#' @param num.hap When haplotype = TRUE, you can set the number of haplotypes which you expect.
#'           Then similar arrays are considered as the same haplotype, and then make kernel(K.SNP) whose dimension is num.hap x num.hap.
#'           When num.hap = NULL (default), num.hap will be set as the maximum number which reflects the difference between lines.
#' @param test.effect Effect of each marker to test. You can choose "test.effect" from "additive", "dominance" and "additive+dominance".
#' You also can choose more than one effect, for example, test.effect = c("additive", "aditive+dominance")
#' @param window.size.half This argument decides how many SNPs (around the SNP you want to test) are used to calculated K.SNP.
#' More precisely, the number of SNPs will be 2 * window.size.half + 1.
#' @param window.slide This argument determines how often you test markers. If window.slide = 1, every marker will be tested.
#' If you want to perform SNP set by bins, please set window.slide = 2 * window.size.half + 1.
#' @param chi0.mixture RAINBOWR assumes the deviance is considered to follow a x chisq(df = 0) + (1 - a) x chisq(df = r).
#' where r is the degree of freedom.
#' The argument chi0.mixture is a (0 <= a < 1), and default is 0.5.
#' @param gene.set If you have information of gene (or haplotype block), you can use it to perform kernel-based GWAS.
#'            You should assign your gene information to gene.set in the form of a "data.frame" (whose dimension is (the number of gene) x 2).
#'            In the first column, you should assign the gene name. And in the second column, you should assign the names of each marker,
#'            which correspond to the marker names of "geno" argument.
#' @param weighting.center In kernel-based GWAS, weights according to the Gaussian distribution (centered on the tested SNP) are taken into account when calculating the kernel if Rainbow = TRUE.
#'           If weighting.center = FALSE, weights are not taken into account.
#' @param weighting.other You can set other weights in addition to weighting.center. The length of this argument should be equal to the number of SNPs.
#'           For example, you can assign SNP effects from the information of gene annotation.
#' @param sig.level Significance level for the threshold. The default is 0.05.
#' @param method.thres Method for detemining threshold of significance. "BH" and "Bonferroni are offered.
#' @param plot.qq If TRUE, draw qq plot.
#' @param plot.Manhattan If TRUE, draw manhattan plot.
#' @param plot.method If this argument = 1, the default manhattan plot will be drawn.
#' If this argument = 2, the manhattan plot with axis based on Position (bp) will be drawn.
#'  Also, this plot's color is changed by all chromosomes.
#' @param plot.col1 This argument determines the color of the manhattan plot.
#'  You should substitute this argument as color vector whose length is 2.
#'  plot.col1[1] for odd chromosomes and plot.col1[2] for even chromosomes
#' @param plot.col2 Color of the manhattan plot. color changes with chromosome and it starts from plot.col2 + 1
#' (so plot.col2 = 1 means color starts from red.)
#' @param plot.type  This argument determines the type of the manhattan plot. See the help page of "plot".
#' @param plot.pch This argument determines the shape of the dot of the manhattan plot. See the help page of "plot".
#' @param saveName When drawing any plot, you can save plots in png format. In saveName, you should substitute the name you want to save.
#' When saveName = NULL, the plot is not saved.
#' @param main.qq The title of qq plot. If this argument is NULL, trait name is set as the title.
#' @param main.man The title of manhattan plot. If this argument is NULL, trait name is set as the title.
#' @param plot.add.last If saveName is not NULL and this argument is TRUE, then you can add lines or dots to manhattan plots.
#' However, you should also write "dev.off()" after adding something.
#' @param optimizer The function used in the optimization process. We offer "optim", "optimx", and "nlminb" functions.
#' @param return.EMM.res When return.EMM.res = TRUE, the results of equation of mixed models are included in the result of RGWAS.
#' @param thres If thres = TRUE, the threshold of the manhattan plot is included in the result of RGWAS.
#' When return.EMM.res or thres is TRUE, the results will be "list" class.
#' @param verbose If this argument is TRUE, messages for the current steps will be shown.
#' @param verbose2 If this argument is TRUE, welcome message will be shown.
#' @param count When count is TRUE, you can know how far RGWAS has ended with percent display.
#' @param time When time is TRUE, you can know how much time it took to perform RGWAS.
#'
#'
#' @return
#' \describe{
#' \item{$D}{Dataframe which contains the information of the map you input and the results of RGWAS (-log10(p)) which correspond to the map.
#' If there are more than one test.effects, then multiple lists for each test.effect are returned respectively.}
#' \item{$thres}{A vector which contains the information of threshold determined by FDR = 0.05.}
#' \item{$EMM.res}{This output is a list which contains the information about the results of "EMM" perfomed at first in regular GWAS.
#'  If you want to know details, see the description for the function "EMM1" or "EMM2".}
#' }
#'
#'
#' @details P-value for each SNP-set is calculated by performing the LR test
#' or the score test (Lippert et al., 2014).
#'
#' In the LR test, first, the function solves the multi-kernel mixed model and
#' calaculates the maximum restricted log likelihood.
#' Then it performs the LR test by using the fact that the deviance
#'
#' \deqn{D = 2 \times (LL _ {alt} - LL _ {null})}
#'
#' follows the chi-square distribution.
#'
#' In the score test, the maximization of the likelihood is only performed for the null model.
#' In other words, the function calculates the score statistic
#' without solving the multi-kernel mixed model for each SNP-set.
#' Then it performs the score test by using the fact that the score statistic follows the chi-square distribution.
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
#' Zhou, X. and Stephens, M. (2012) Genome-wide efficient mixed-model analysis for association studies. Nat Genet. 44(7): 821-824.
#'
#' Listgarten, J. et al. (2013) A powerful and efficient set test for genetic markers that handles confounders. Bioinformatics. 29(12): 1526-1533.
#'
#' Lippert, C. et al. (2014) Greater power and computational efficiency for kernel-based association testing of sets of genetic variants. Bioinformatics. 30(22): 3206-3214.
#'
#' @example R/examples/RGWAS.multisnp_example.R
#'
#'
#'
RGWAS.multisnp <- function(pheno, geno, ZETA = NULL, covariate = NULL, covariate.factor = NULL,
                           structure.matrix = NULL, n.PC = 0, min.MAF = 0.02, test.method = "LR", n.core = 1,
                           kernel.method = "linear", kernel.h = "tuned", haplotype = TRUE, num.hap = NULL,
                           test.effect = "additive", window.size.half = 5, window.slide = 1, chi0.mixture = 0.5,
                           gene.set = NULL, weighting.center = TRUE, weighting.other = NULL,
                           sig.level = 0.05, method.thres = "BH", plot.qq = TRUE, plot.Manhattan = TRUE, plot.method = 1,
                           plot.col1 = c("dark blue", "cornflowerblue"), plot.col2 = 1,
                           plot.type = "p", plot.pch = 16, saveName = NULL, main.qq = NULL,
                           main.man = NULL, plot.add.last = FALSE, return.EMM.res = FALSE, optimizer = "nlminb",
                           thres = TRUE, verbose = TRUE, verbose2 = FALSE, count = TRUE, time = TRUE){

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
  if((kernel.method == "linear") & (length(test.effect) >= 2)){
    all.scores <- rep(list(matrix(0, nrow = n.scores, ncol = n.pheno)), length(test.effect))
    names(all.scores) <- test.effect

    for(test.effect.no in 1:length(test.effect)){
      colnames(all.scores[[test.effect.no]]) <- trait.names
    }
    thresholds <- matrix(NA, nrow = length(test.effect), ncol = n.pheno)
    rownames(thresholds) <- test.effect
    colnames(thresholds) <- trait.names
  }else{
    all.scores <- matrix(0, nrow = n.scores, ncol = n.pheno)
    colnames(all.scores) <- trait.names
    thresholds <- matrix(NA, nrow = 1, ncol = n.pheno)
    rownames(thresholds) <- kernel.method
    colnames(thresholds) <- trait.names
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

    X.now <- X[not.NA, , drop = FALSE]
    ZETA.now <- lapply(ZETA, function(x) list(Z = x$Z[not.NA, ], K = x$K))

    if (sum(is.na(M)) == 0) {
      M.now <- Z.A[not.NA, ] %*% M
    } else {
      M.now <- M[apply(Z.A[not.NA, ], 1, function(x) which(x == 1)), ]
    }

    p <- ncol(X.now)
    m <- ncol(Z.A)

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
        if(test.method == "LR"){
          scores <- score.calc.LR.MC(M.now = M.now, y = y, X.now = X.now, ZETA.now = ZETA.now, LL0 = LL0,
                                     eigen.SGS = eigen.SGS, eigen.G = eigen.G, n.core = n.core, map = map,
                                     kernel.method = kernel.method, kernel.h = kernel.h, haplotype = haplotype,
                                     num.hap = num.hap, test.effect = test.effect, window.size.half = window.size.half,
                                     window.slide = window.slide, chi0.mixture = chi0.mixture, optimizer = optimizer,
                                     weighting.center = weighting.center, weighting.other = weighting.other,
                                     gene.set = gene.set, min.MAF = min.MAF, count = count)
        }else{
          scores <- score.calc.score.MC(M.now = M.now, ZETA.now = ZETA.now, y = y,
                                        LL0 = LL0, Gu = Gu, Ge = Ge, P0 = P0, n.core = n.core, map = map,
                                        kernel.method = kernel.method, kernel.h = kernel.h, haplotype = haplotype,
                                        num.hap = num.hap, test.effect = test.effect, window.size.half = window.size.half,
                                        window.slide = window.slide, chi0.mixture = chi0.mixture,
                                        weighting.center = weighting.center, weighting.other = weighting.other,
                                        gene.set = gene.set, min.MAF = min.MAF, count = count)
        }
    }else {
      if(test.method == "LR"){
        scores <- score.calc.LR(M.now = M.now, y = y, X.now = X.now, ZETA.now = ZETA.now, LL0 = LL0,
                                eigen.SGS = eigen.SGS, eigen.G = eigen.G, map = map, optimizer = optimizer,
                                kernel.method = kernel.method, kernel.h = kernel.h, haplotype = haplotype,
                                num.hap = num.hap, test.effect = test.effect, window.size.half = window.size.half,
                                window.slide = window.slide, chi0.mixture = chi0.mixture,
                                weighting.center = weighting.center, weighting.other = weighting.other,
                                gene.set = gene.set, min.MAF = min.MAF, count = count)
      }else{
        scores <- score.calc.score(M.now = M.now, ZETA.now = ZETA.now, y = y,
                                   LL0 = LL0, Gu = Gu, Ge = Ge, P0 = P0, map = map,
                                   kernel.method = kernel.method, kernel.h = kernel.h, haplotype = haplotype,
                                   num.hap = num.hap, test.effect = test.effect, window.size.half = window.size.half,
                                   window.slide = window.slide, chi0.mixture = chi0.mixture,
                                   weighting.center = weighting.center, weighting.other = weighting.other,
                                   gene.set = gene.set, min.MAF = min.MAF, count = count)
      }
    }

    if(is.null(gene.set)){
      window.centers <- as.numeric(rownames(scores))
      map2 <- map[window.centers, ]
    }else{
      if (verbose) {
        print("Now generating map for gene set. Please wait.")
      }
      map20 <- genesetmap(map = map, gene.set = gene.set, cumulative = TRUE)
      map2 <- map20[, 1:3]
      cum.pos.set.mean <- c(map20[, 4])
    }

    if((kernel.method == "linear") & (length(test.effect) >= 2)){
      for(test.effect.no in 1:length(test.effect)){
        if (plot.qq) {
          if (verbose) {
            print("Now Plotting (Q-Q plot). Please Wait.")
          }
          if(is.null(saveName)){
            if (length(grep("RStudio", names(dev.cur()))) == 0) {
              if (dev.cur() == dev.next()) {
                dev.new()
              }
              else {
                dev.set(dev.next())
              }
            }
            qq(scores[, test.effect.no])
            if(is.null(main.qq)){
              title(main = paste(trait.name, test.effect[test.effect.no]))
            }else{
              title(main = main.qq)
            }
          }else{
            png(paste0(saveName, trait.name, "_qq_kernel_", test.effect[test.effect.no], ".png"))
            qq(scores[, test.effect.no])
            if(is.null(main.qq)){
              title(main = paste(trait.name, test.effect[test.effect.no]))
            }else{
              title(main = main.qq)
            }
            dev.off()
          }
        }


        if (plot.Manhattan) {
          if (verbose) {
            print("Now Plotting (Manhattan plot). Please Wait.")
          }
          if(is.null(saveName)){
            if (length(grep("RStudio", names(dev.cur()))) == 0) {
              if (dev.cur() == dev.next()) {
                dev.new()
              }
              else {
                dev.set(dev.next())
              }
            }
            if(plot.method == 1){
              manhattan(input = cbind(map2, scores[, test.effect.no]), sig.level = sig.level, method.thres = method.thres, plot.col1 = plot.col1,
                        plot.type = plot.type, plot.pch = plot.pch)
            }else{
              manhattan2(input = cbind(map2, scores[, test.effect.no]), sig.level = sig.level, method.thres = method.thres, plot.col2 = plot.col2,
                         plot.type = plot.type, plot.pch = plot.pch, cum.pos = cum.pos.set.mean)
            }
            if(is.null(main.man)){
              title(main = paste(trait.name, test.effect[test.effect.no]))
            }else{
              title(main = main.man)
            }
          }else{
            png(paste0(saveName, trait.name, "_manhattan_kernel", test.effect[test.effect.no], ".png"), width = 800)
            if(plot.method == 1){
              manhattan(input = cbind(map2, scores[, test.effect.no]), sig.level = sig.level, method.thres = method.thres, plot.col1 = plot.col1,
                        plot.type = plot.type, plot.pch = plot.pch)
            }else{
              manhattan2(input = cbind(map2, scores[, test.effect.no]), sig.level = sig.level, method.thres = method.thres, plot.col2 = plot.col2,
                         plot.type = plot.type, plot.pch = plot.pch, cum.pos = cum.pos.set.mean)
            }
            if(is.null(main.man)){
              title(main = paste(trait.name, test.effect[test.effect.no]))
            }else{
              title(main = main.man)
            }
            if(!(plot.add.last & (pheno.no == n.pheno))){
              dev.off()
            }
          }
        }
        all.scores[[test.effect.no]][, pheno.no] <- scores[, test.effect.no]
        threshold <- try(CalcThreshold(cbind(map2, scores[, test.effect.no]), sig.level = sig.level, method = method.thres), silent = TRUE)
        if("try-error" %in% class(threshold)){
          threshold <- NA
        }
        thresholds[test.effect.no, pheno.no] <- threshold

      }
    }else{
      if (plot.qq) {
        if (verbose) {
          print("Now Plotting (Q-Q plot). Please Wait.")
        }
        if(is.null(saveName)){
          if (length(grep("RStudio", names(dev.cur()))) == 0) {
            if (dev.cur() == dev.next()) {
              dev.new()
            }
            else {
              dev.set(dev.next())
            }
          }
          qq(scores)
          if(is.null(main.qq)){
            title(main = trait.name)
          }else{
            title(main = main.qq)
          }
        }else{
          png(paste0(saveName, trait.name, "_qq_kernel.png"))
          qq(scores)
          if(is.null(main.qq)){
            title(main = trait.name)
          }else{
            title(main = main.qq)
          }
          dev.off()
        }
      }


      if (plot.Manhattan) {
        if (verbose) {
          print("Now Plotting (Manhattan plot). Please Wait.")
        }
        if(is.null(saveName)){
          if (length(grep("RStudio", names(dev.cur()))) == 0) {
            if (dev.cur() == dev.next()) {
              dev.new()
            }
            else {
              dev.set(dev.next())
            }
          }
          if(plot.method == 1){
            manhattan(input = cbind(map2, scores), sig.level = sig.level, method.thres = method.thres, plot.col1 = plot.col1,
                      plot.type = plot.type, plot.pch = plot.pch)
          }else{
            manhattan2(input = cbind(map2, scores), sig.level = sig.level, method.thres = method.thres, plot.col2 = plot.col2,
                       plot.type = plot.type, plot.pch = plot.pch, cum.pos = cum.pos.set.mean)
          }
          if(is.null(main.man)){
            title(main = trait.name)
          }else{
            title(main = main.man)
          }
        }else{
          png(paste0(saveName, trait.name, "_manhattan_kernel.png"), width = 800)
          if(plot.method == 1){
            manhattan(input = cbind(map2, scores), sig.level = sig.level, method.thres = method.thres, plot.col1 = plot.col1,
                      plot.type = plot.type, plot.pch = plot.pch)
          }else{
            manhattan2(input = cbind(map2, scores), sig.level = sig.level, method.thres = method.thres, plot.col2 = plot.col2,
                       plot.type = plot.type, plot.pch = plot.pch, cum.pos = cum.pos.set.mean)
          }
          if(is.null(main.man)){
            title(main = trait.name)
          }else{
            title(main = main.man)
          }
          if(!(plot.add.last & (pheno.no == n.pheno))){
            dev.off()
          }
        }
      }
      all.scores[, pheno.no] <- scores
      threshold <- try(CalcThreshold(cbind(map2, scores), sig.level = sig.level, method = method.thres), silent = TRUE)
      if("try-error" %in% class(threshold)){
        threshold <- NA
      }
      thresholds[, pheno.no] <- threshold
    }
  }

  if((kernel.method == "linear") & (length(test.effect) >= 2)){
    res.Data <- lapply(all.scores, function(x){
      cbind(map2, x)
    })
  }else{
    res.Data <- cbind(map2, all.scores)
  }




  if(thres){
    end <- Sys.time()

    if(time){
      print(end - start)
    }

    if(return.EMM.res){
      return(list(D = res.Data, thres = thresholds,
                  EMM.res = EMM.res0))
    }else{
      return(list(D = res.Data, thres = thresholds))
    }
  }else{
    end <- Sys.time()

    if(time){
      print(end - start)
    }

    if(return.EMM.res){
      return(list(D = res.Data, EMM.res = EMM.res0))
    }else{
      return(res.Data)
    }
  }
}
