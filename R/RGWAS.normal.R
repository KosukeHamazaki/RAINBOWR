#' Perform normal GWAS (test each single SNP)
#'
#' @description This function performs single-SNP GWAS (genome-wide association studies). The model of GWAS is
#'
#' \deqn{y = X \beta + S _ {i} \alpha _ {i} + Q v +  Z u + \epsilon,}
#'
#' where \eqn{y} is the vector of phenotypic values,
#' \eqn{X \beta}, \eqn{S _ {i} \alpha _ {i}}, \eqn{Q v} are the terms of fixed effects,
#' \eqn{Z u} is the term of random effects and \eqn{e} is the vector of residuals.
#' \eqn{X \beta} indicates all of the fixed effects other than the effect of SNPs
#' to be tested and of population structure, and often this term also plays
#' a role as an intercept. For \eqn{S _ {i} \alpha _ {i}}, \eqn{S _ {i}}
#' is the ith marker of genotype data and \eqn{\alpha _ {i}} is the effect of that marker.
#' \eqn{Q v} is the term to correct the effect of population structure.
#' \eqn{Z u} is the term of polygenetic effects, and suppose that \eqn{u}
#' follows the multivariate normal distribution whose variance-covariance
#'  matrix is the genetic covariance matrix. \eqn{u \sim MVN (0, K \sigma_{u}^{2})}.
#' Finally, the residual term is assumed to identically and independently follow
#'  a normal distribution as shown in the following equation.
#'  \eqn{e \sim MVN (0, I \sigma_{e}^{2})}.
#'
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
#' @param P3D When P3D = TRUE, variance components are estimated by REML only once, without any markers in the model.
#' When P3D = FALSE, variance components are estimated by REML for each marker separately.
#' @param n.core Setting n.core > 1 will enable parallel execution on a machine with multiple cores.
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
#' \item{$D}{Dataframe which contains the information of the map you input and the results of RGWAS (-log10(p)) which correspond to the map.}
#' \item{$thres}{A vector which contains the information of threshold determined by FDR = 0.05.}
#' \item{$EMM.res}{This output is a list which contains the information about the results of "EMM" perfomed at first in regular GWAS.
#'  If you want to know details, see the description for the function "EMM1" or "EMM2".}
#' }
#'
#'
#' @details P-value for each marker is calculated by performing F-test
#' against the F-value as follows (Kennedy et al., 1992).
#'
#' \deqn{F = \frac { ( L' \hat { b } )' [ L' ( X' H ^ { - 1 } X ) ^ { - 1 }
#' L ] ^ { - 1 } ( L' \hat { b } ) } { f \hat { \sigma }_ { u } ^ { 2 } },}
#'
#' where \eqn{b} is the vector of coefficients of the fixed effects, which combines
#' \eqn{\beta}, \eqn{\alpha _ {i}}, \eqn{v} in the horizontal direction and \eqn{L}
#' is a matrix to indicate which effects in \eqn{b} are tested.
#' \eqn{H} is calculated by dividing the estimated variance-covariance
#' matrix for the phenotypic values by \eqn{\sigma _ { u } ^ { 2 }},
#' and is calculated by \eqn{H = Z K Z' + \hat{\lambda} I}.
#' \eqn{\hat{\lambda}} is the maximum likelihood estimator
#' of the ratio between the residual variance and the additive genetic variance.
#' \eqn{\hat{b}} is the maximum likelihood estimator of \eqn{b}
#' and is calculated by \eqn{\hat { b } = ( X' H ^ { - 1 } X ) ^ { - 1 } X' H ^ { - 1 } y }.
#' \eqn{f} is the number of the fixed effects to be tested, and
#' \eqn{\hat { \sigma }_ { u } ^ { 2 }} is estimated by the following formula.
#' \deqn{\hat { \sigma }_ { u } ^ { 2 } = \frac { ( y - X  \hat { b } )' H ^ { - 1 } ( y - X  \hat { b } ) } { n - p },}
#' where \eqn{n} is the sample size and \eqn{p} is the number of the all fixed effects.
#' We calculated each p-value using the fact that the above F-value follows
#' the F distribution with the degree of freedom (\eqn{f},\eqn{n - p}).
#'
#' @references
#' Kennedy, B.W., Quinton, M. and van Arendonk, J.A. (1992) Estimation of effects of single genes on quantitative traits. J Anim Sci. 70(7): 2000-2012.
#'
#' Storey, J.D. and Tibshirani, R. (2003) Statistical significance for genomewide studies. Proc Natl Acad Sci. 100(16): 9440-9445.
#'
#' Yu, J. et al. (2006) A unified mixed-model method for association mapping that accounts for multiple levels of relatedness. Nat Genet. 38(2): 203-208.
#'
#' Kang, H.M. et al. (2008) Efficient Control of Population Structure in Model Organism Association Mapping. Genetics. 178(3): 1709-1723.
#'
#' Kang, H.M. et al. (2010) Variance component model to account for sample structure in genome-wide association studies. Nat Genet. 42(4): 348-354.
#'
#' Zhang, Z. et al. (2010) Mixed linear model approach adapted for genome-wide association studies. Nat Genet. 42(4): 355-360.
#'
#' Endelman, J.B. (2011) Ridge Regression and Other Kernels for Genomic Selection with R Package rrBLUP. Plant Genome J. 4(3): 250.
#'
#' Endelman, J.B. and Jannink, J.L. (2012) Shrinkage Estimation of the Realized Relationship Matrix. G3 Genes, Genomes, Genet. 2(11): 1405-1413.
#'
#' Zhou, X. and Stephens, M. (2012) Genome-wide efficient mixed-model analysis for association studies. Nat Genet. 44(7): 821-824.
#'
#' @example R/examples/RGWAS.normal_example.R
#'
#'
#'
RGWAS.normal <- function(pheno, geno, ZETA = NULL, covariate = NULL, covariate.factor = NULL,
                         structure.matrix = NULL, n.PC = 0, min.MAF = 0.02, P3D = TRUE, n.core = 1,
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
  all.scores <- matrix(0, nrow = n.mark, ncol = n.pheno)
  trait.names <- colnames(pheno)[pheno.ix]
  colnames(all.scores) <- trait.names

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
    
    #### Calculate Hinv at first ####
    if(P3D){
      if(length(ZETA) > 1){
        EMM.res0 <- EM3.cpp(y = y, X0 = X.now, ZETA = ZETA.now,
                            n.thres = 450, REML = TRUE, pred = FALSE)
      }else{
        EMM.res0 <- EMM.cpp(y = y, X = X.now, ZETA = ZETA.now,
                            n.thres = 450, REML = TRUE)
      }

      Hinv <- EMM.res0$Hinv
      eigen.G <- NULL
      if (verbose) {
        print("Variance components estimated. Testing markers.")
      }
    }else{
      spI <- diag(n)
      if(length(ZETA) > 1){
        EMM.res0 <- EM3.cpp(y = y, X0 = X.now, ZETA = ZETA.now,
                            n.thres = 450, REML = TRUE, pred = FALSE)
        weights <- EMM.res0$weights
        eigen.G <- spectralG.cpp(ZETA = ZETA.now, X = X.now, weights = weights,
                                 return.G = TRUE, return.SGS = FALSE)[[1]]
      }else{
        eigen.G <- spectralG.cpp(ZETA = ZETA.now, X = X.now, return.G = TRUE,
                                 return.SGS = FALSE)[[1]]
      }
      Hinv <- NULL
    }




    #### Calculating the value of -log10(p) for each SNPs ####
    if ((n.core > 1) & requireNamespace("parallel", quietly = TRUE)) {
      scores <- score.calc.MC(M.now = M.now, ZETA.now = ZETA.now, y = y,
                   X.now = X.now, Hinv = Hinv, P3D = P3D, eigen.G = eigen.G,
                   optimizer = optimizer, min.MAF = min.MAF, count = count)
    } else {
      scores <- score.calc(M.now, ZETA.now = ZETA.now, y = y, X.now = X.now, Hinv = Hinv,
                           optimizer = optimizer, P3D = P3D, eigen.G = eigen.G, min.MAF = min.MAF, count = count)
    }

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
        png(paste0(saveName, trait.name, "_qq.png"))
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
          manhattan(input = cbind(map, scores), sig.level = sig.level, method.thres = method.thres, plot.col1 = plot.col1,
                    plot.type = plot.type, plot.pch = plot.pch)
        }else{
          manhattan2(input = cbind(map, scores), sig.level = sig.level, method.thres = method.thres, plot.col2 = plot.col2,
                     plot.type = plot.type, plot.pch = plot.pch)
        }
        if(is.null(main.man)){
          title(main = trait.name)
        }else{
          title(main = main.man)
        }
      }else{
        png(paste0(saveName, trait.name, "_manhattan.png"), width = 800)
        if(plot.method == 1){
          manhattan(input = cbind(map, scores), sig.level = sig.level, method.thres = method.thres, plot.col1 = plot.col1,
                    plot.type = plot.type, plot.pch = plot.pch)
        }else{
          manhattan2(input = cbind(map, scores), sig.level = sig.level, method.thres = method.thres, plot.col2 = plot.col2,
                     plot.type = plot.type, plot.pch = plot.pch)
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
  }

  if(thres){
    thresholds <- rep(NA, n.pheno)
    for(pheno.no in 1:n.pheno){
      thresholds[pheno.no] <- try(CalcThreshold(cbind(map, all.scores[, pheno.no]), sig.level = sig.level, method = method.thres), silent = TRUE)
    }

    end <- Sys.time()

    if(time){
      print(end - start)
    }

    if(return.EMM.res){
      return(list(D = data.frame(map, all.scores), thres = thresholds,
                  EMM.res = EMM.res0))
    }else{
      return(list(D = data.frame(map, all.scores), thres = thresholds))
    }
  }else{
    end <- Sys.time()

    if(time){
      print(end - start)
    }

    if(return.EMM.res){
      return(list(D = data.frame(map, all.scores), EMM.res = EMM.res0))
    }else{
      return(data.frame(map, all.scores))
    }
  }
}
