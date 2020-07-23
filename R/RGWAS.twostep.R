#' Perform normal GWAS (genome-wide association studies) first, then perform SNP-set GWAS for relatively significant markers
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
#' @param check.size This argument determines how many SNPs (around the SNP detected by normal GWAS) you will recalculate -log10(p).
#' @param check.gene.size This argument determines how many genes (around the genes detected by normal GWAS) you will recalculate -log10(p).
#' This argument is valid only when you assign "gene.set" argument.
#' @param kernel.percent This argument determines how many SNPs are detected by normal GWAS.
#' For example, when kernel.percent = 0.1, SNPs whose value of -log10(p) is in the top 0.1 percent are chosen as candidate for recalculation by SNP-set GWAS.
#' @param GWAS.res.first If you have already performed normal GWAS and have the result, you can skip performing normal GWAS.
#' @param P3D When P3D = TRUE, variance components are estimated by REML only once, without any markers in the model.
#' When P3D = FALSE, variance components are estimated by REML for each marker separately.
#' @param test.method.1 RGWAS supports two methods to test effects of each SNP-set for 1st GWAS.
#' \describe{
#' \item{"normal"}{Normal GWAS (default).}
#' \item{"score"}{Score test, much faster than LR, but sometimes overestimate -log10(p).}
#' }
#' @param test.method.2 RGWAS supports two methods to test effects of each SNP-set for 2nd GWAS.
#' \describe{
#' \item{"LR"}{Likelihood-ratio test, relatively slow, but accurate (default).}
#' \item{"score"}{Score test, much faster than LR, but sometimes overestimate -log10(p).}
#' }
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
#' @param test.effect.1 Effect of each marker to test for 1st GWAS. You can choose "test.effect" from "additive", "dominance" and "additive+dominance".
#' you can assign only one test effect for the 1st GWAS!
#' @param test.effect.2 Effect of each marker to test for 2nd GWAS. You can choose "test.effect" from "additive", "dominance" and "additive+dominance".
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
#' @param plot.qq.1 If TRUE, draw qq plot for normal GWAS.
#' @param plot.Manhattan.1 If TRUE, draw manhattan plot for normal GWAS.
#' @param plot.qq.2 If TRUE, draw qq plot for SNP-set GWAS.
#' @param plot.Manhattan.2 If TRUE, draw manhattan plot for SNP-set GWAS.
#' @param plot.method If this argument = 1, the default manhattan plot will be drawn.
#' If this argument = 2, the manhattan plot with axis based on Position (bp) will be drawn.
#' Also, this plot's color is changed by all chromosomes.
#' @param plot.col1 This argument determines the color of the manhattan plot.
#'  You should substitute this argument as color vector whose length is 2.
#'  plot.col1[1] for odd chromosomes and plot.col1[2] for even chromosomes
#' @param plot.col2 Color of the manhattan plot. color changes with chromosome and it starts from plot.col2 + 1
#' (so plot.col2 = 1 means color starts from red.)
#' @param plot.col3 Color of the points of manhattan plot which are added after the reestimation by SNP-set method.
#' You should substitute this argument as color vector whose length is 2.
#' plot.col3[1] for odd chromosomes and plot.col3[2] for even chromosomes.
#' @param plot.type  This argument determines the type of the manhattan plot. See the help page of "plot".
#' @param plot.pch This argument determines the shape of the dot of the manhattan plot. See the help page of "plot".
#' @param saveName When drawing any plot, you can save plots in png format. In saveName, you should substitute the name you want to save.
#' When saveName = NULL, the plot is not saved.
#' @param main.qq.1 The title of qq plot for normal GWAS. If this argument is NULL, trait name is set as the title.
#' @param main.man.1 The title of manhattan plot for normal GWAS. If this argument is NULL, trait name is set as the title.
#' @param main.qq.2 The title of qq plot for SNP-set GWAS. If this argument is NULL, trait name is set as the title.
#' @param main.man.2 The title of manhattan plot for SNP-set GWAS. If this argument is NULL, trait name is set as the title.
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
#' -log10(p) by normal GWAS and recalculated -log10(p) by SNP-set GWAS will be obtained.
#' If there are more than one test.effects, then multiple lists for each test.effect are returned respectively.
#' }
#' \item{$thres}{A vector which contains the information of threshold determined by FDR = 0.05.}
#' \item{$EMM.res}{This output is a list which contains the information about the results of "EMM" perfomed at first in normal GWAS.
#'  If you want to know details, see the description for the function "EMM1" or "EMM2".}
#' }
#'
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
#' Listgarten, J. et al. (2013) A powerful and efficient set test for genetic markers that handles confounders. Bioinformatics. 29(12): 1526-1533.
#'
#' Lippert, C. et al. (2014) Greater power and computational efficiency for kernel-based association testing of sets of genetic variants. Bioinformatics. 30(22): 3206-3214.
#'
#'
#' @example R/examples/RGWAS.twostep_example.R
#'
#'
#'
RGWAS.twostep <- function(pheno, geno, ZETA = NULL, covariate = NULL, covariate.factor = NULL,
                          structure.matrix = NULL, n.PC = 0, min.MAF = 0.02, n.core = 1,
                          check.size = 40, check.gene.size = 4, kernel.percent = 0.1, GWAS.res.first = NULL,
                          P3D = TRUE, test.method.1 = "normal", test.method.2 = "LR",
                          kernel.method = "linear", kernel.h = "tuned", haplotype = TRUE,
                          num.hap = NULL, test.effect.1 = "additive", test.effect.2 = "additive",
                          window.size.half = 5, window.slide = 1, chi0.mixture = 0.5, optimizer = "nlminb",
                          gene.set = NULL, weighting.center = TRUE, weighting.other = NULL,
                          sig.level = 0.05, method.thres = "BH", plot.qq.1 = TRUE, plot.Manhattan.1 = TRUE,
                          plot.qq.2 = TRUE, plot.Manhattan.2 = TRUE, plot.method = 1,
                          plot.col1 = c("dark blue", "cornflowerblue"), plot.col2 = 1,
                          plot.col3 = c("red3", "orange3"), plot.type = "p",
                          plot.pch = 16, saveName = NULL, main.qq.1 = NULL,
                          main.man.1 = NULL, main.qq.2 = NULL, main.man.2 = NULL,
                          plot.add.last = FALSE, return.EMM.res = FALSE, thres = TRUE,
                          verbose = TRUE, verbose2 = FALSE, count = TRUE, time = TRUE){

  start <- Sys.time()
  if(is.null(GWAS.res.first)){
    if (verbose) {
      print("The 1st step: Performing 1st GWAS (for screening)!")
    }
    if(length(test.effect.1) >= 2){
      stop("Sorry, you can assign only one test effect for the 1st GWAS!!")
    }
    if(test.method.1 == "normal"){
      GWAS.res.first <- RGWAS.normal(pheno = pheno, geno = geno, ZETA = ZETA, covariate = covariate,
                                     covariate.factor = covariate.factor, structure.matrix = structure.matrix,
                                     n.PC = n.PC, min.MAF = min.MAF, P3D = P3D, n.core = n.core,
                                     sig.level = sig.level, method.thres = method.thres, plot.qq = plot.qq.1, plot.Manhattan = plot.Manhattan.1,
                                     plot.method = plot.method, plot.col1 = plot.col1, plot.col2 = plot.col2,
                                     plot.type = plot.type, plot.pch = plot.pch, saveName = saveName, optimizer = optimizer,
                                     main.qq = main.qq.1, main.man = main.man.1, plot.add.last = FALSE, return.EMM.res = FALSE,
                                     thres = FALSE, verbose = verbose, verbose2 = verbose2, count = count, time = time)
    }else{
      GWAS.res.first <- RGWAS.multisnp(pheno = pheno, geno = geno, ZETA = ZETA, covariate = covariate,
                                       covariate.factor = covariate.factor, structure.matrix = structure.matrix,
                                       n.PC = n.PC, min.MAF = min.MAF, test.method = test.method.1, n.core = n.core,
                                       kernel.method = kernel.method, kernel.h = kernel.h, haplotype = haplotype,
                                       num.hap = num.hap, test.effect = test.effect.1, window.size.half = window.size.half,
                                       window.slide = window.slide, chi0.mixture = chi0.mixture, gene.set = gene.set,
                                       weighting.center = weighting.center, weighting.other = weighting.other,
                                       sig.level = sig.level, method.thres = method.thres, plot.qq = FALSE, plot.Manhattan = FALSE,
                                       plot.method = plot.method, plot.col1 = plot.col1, plot.col2 = plot.col2,
                                       plot.type = plot.type, plot.pch = plot.pch, saveName = saveName,
                                       main.qq = main.qq.2, main.man = main.man.2, plot.add.last = FALSE,
                                       return.EMM.res = FALSE, optimizer = optimizer,
                                       thres = FALSE, verbose = verbose, verbose2 = verbose2, count = count, time = time)
    }
  }else{
    if (verbose) {
      print("The 1st step has already finished because you input 'GWAS.res.first'.")
    }
  }

  n.pheno <- ncol(GWAS.res.first) - 3
  trait.names <- colnames(GWAS.res.first)[4:(4 + n.pheno - 1)]
  map <- geno[, 1:3]

  if((kernel.method == "linear") & (length(test.effect.2) >= 2)){
    thresholds <- matrix(NA, nrow = length(test.effect.2), ncol = n.pheno)
    thresholds.correction <- matrix(NA, nrow = length(test.effect.2), ncol = n.pheno)
    rownames(thresholds) <- rep("normal", length(test.effect.2))
    rownames(thresholds.correction) <- test.effect.2
    colnames(thresholds) <- colnames(thresholds.correction) <- trait.names
  }else{
    thresholds <- thresholds.correction <- matrix(NA, nrow = 1, ncol = n.pheno)
    rownames(thresholds) <- "normal"
    rownames(thresholds.correction) <- kernel.method
    colnames(thresholds) <- colnames(thresholds.correction) <- trait.names
  }

  if((kernel.method == "linear") & (length(test.effect.2) >= 2)){
    res.all <- rep(list(GWAS.res.first), length(test.effect.2))
  }else{
    res.all <- GWAS.res.first
  }



  for(pheno.no in 1:n.pheno){
    trait.name <- trait.names[pheno.no]
    pheno.now <- pheno[, c(1, pheno.no + 1)]
    GWAS.res.first.now <- GWAS.res.first[, c(1:3, pheno.no + 3)]
    pval.first <- GWAS.res.first[, pheno.no + 3]
    ord.pval.first <- order(pval.first, decreasing = TRUE)
    ord.pval.ker.percent.0 <- ord.pval.first[1:round(length(pval.first) * (kernel.percent / 100), 0)]
    ord.pval.ker.percent <- as.numeric(rownames(GWAS.res.first)[ord.pval.ker.percent.0])

    if(is.null(gene.set)){
      check.obj <- "SNPs"
      check.size.half <- check.size / 2
      checks.mat <- matrix(NA, nrow = length(ord.pval.ker.percent) * (check.size + 1), ncol = length(ord.pval.ker.percent))
      for(check.no in 1:length(ord.pval.ker.percent)){
        check <- sort(ord.pval.ker.percent)[check.no]
        checks.now <- (check - check.size.half):(check + check.size.half)
        checks.mat[, check.no] <- checks.now
      }
      checks <- unique(c(checks.mat))
      checks <- checks[(checks >= 1) & (checks <= max(as.numeric(rownames(GWAS.res.first))))]
      n.checks <- length(checks)
      pseudo.chr <- rep(NA, n.checks)
      pseudo.chr[1] <- 1
      for(k in 2:n.checks){
        pseudo.chr.now <- pseudo.chr[k - 1]

        check.diff <- checks[k] - checks[k - 1]
        if(check.diff == 1){
          pseudo.chr[k] <- pseudo.chr.now
        }else{
          pseudo.chr[k] <- pseudo.chr.now + 1
        }
      }
      pseudo.marker <- as.character(map[checks, 1])
      pseudo.pos <- map[checks, 3]
      pseudo.map <- data.frame(marker = pseudo.marker, chr = pseudo.chr,
                               pos = pseudo.pos)
      rownames(pseudo.map) <- checks
      M.check <- geno[checks, -c(1:3)]
      geno.check <- cbind(pseudo.map, M.check)

      gene.set.now <- NULL
    }else{
      check.obj <- "genes"
      if(test.method.1 != "normal"){
        check.size.half <- check.gene.size / 2
        checks.mat <- matrix(NA, nrow = length(ord.pval.ker.percent) * (check.gene.size + 1), ncol = length(ord.pval.ker.percent))
        for(check.no in 1:length(ord.pval.ker.percent)){
          check <- sort(ord.pval.ker.percent)[check.no]
          checks.now <- (check - check.size.half):(check + check.size.half)
          checks.mat[, check.no] <- checks.now
        }

        checks <- unique(c(checks.mat))
        checks <- checks[(checks >= 1) & (checks <= max(as.numeric(rownames(GWAS.res.first))))]
        n.checks <- length(checks)

        gene.names <- as.character(unique(gene.set[, 1]))
        gene.names.now <- gene.names[checks]
      }else{
        check.size.half <- check.size / 2
        checks.mat <- matrix(NA, nrow = length(ord.pval.ker.percent) * (check.size + 1), ncol = length(ord.pval.ker.percent))
        for(check.no in 1:length(ord.pval.ker.percent)){
          check <- sort(ord.pval.ker.percent)[check.no]
          checks.now <- (check - check.size.half):(check + check.size.half)
          checks.mat[, check.no] <- checks.now
        }
        checks <- unique(c(checks.mat))
        checks <- checks[(checks >= 1) & (checks <= max(as.numeric(rownames(GWAS.res.first))))]

        match.gene.list <- match(as.character(gene.set[, 2]), as.character(map[checks, 1]))
        gene.names.now <- unique(as.character(gene.set[!is.na(match.gene.list), 1]))
        n.checks <- length(gene.names.now)
      }
      gene.set.now <- gene.set[as.character(gene.set[, 1]) %in% gene.names.now, ]

      geno.check <- geno
    }


    if (verbose) {
      print(paste("The 2nd step: Recalculating -log10(p) of", trait.name, "for", n.checks, check.obj, "by kernel-based (mutisnp) GWAS."))
    }
    RGWAS.multisnp.res.0 <- RGWAS.multisnp(pheno = pheno.now, geno = geno.check, ZETA = ZETA, covariate = covariate,
                                         covariate.factor = covariate.factor, structure.matrix = structure.matrix,
                                         n.PC = n.PC, min.MAF = min.MAF, test.method = test.method.2, n.core = n.core,
                                         kernel.method = kernel.method, kernel.h = kernel.h, haplotype = haplotype,
                                         num.hap = num.hap, test.effect = test.effect.2, window.size.half = window.size.half,
                                         window.slide = window.slide, chi0.mixture = chi0.mixture, gene.set = gene.set.now,
                                         weighting.center = weighting.center, weighting.other = weighting.other,
                                         sig.level = sig.level, method.thres = method.thres, plot.qq = FALSE, plot.Manhattan = FALSE,
                                         plot.method = plot.method, plot.col1 = plot.col1, plot.col2 = plot.col2,
                                         plot.type = plot.type, plot.pch = plot.pch, saveName = saveName,
                                         main.qq = main.qq.2, main.man = main.man.2, plot.add.last = FALSE,
                                         return.EMM.res = TRUE, optimizer = optimizer,
                                         thres = FALSE, verbose = verbose, count = count, time = time)

    RGWAS.multisnp.res <- RGWAS.multisnp.res.0$D
    EMM.res0 <- RGWAS.multisnp.res.0$EMM.res

    if((kernel.method == "linear") & (length(test.effect.2) >= 2)){
      GWAS.res.merge.list <- lapply(RGWAS.multisnp.res, function(x){
        colnames(x) <- colnames(GWAS.res.first.now)
        if(is.null(gene.set)){
          x[, 1:3] <- map[match(rownames(x), rownames(map)), ]
        }
        GWAS.res.merge.0 <- rbind(x, GWAS.res.first.now)
        GWAS.res.merge <- GWAS.res.merge.0[!duplicated(as.character(GWAS.res.merge.0[, 1])), ]
        ord.GWAS.res.merge <- order(GWAS.res.merge[, 2], GWAS.res.merge[, 3])
        res.correction <- GWAS.res.merge[ord.GWAS.res.merge, ]
        check.here <- match(1:nrow(x), ord.GWAS.res.merge)
        return(list(res = res.correction, check = check.here))
      })


      res.corrections <- rep(list(NA), length(test.effect.2))
      for(test.effect.no in 1:length(test.effect.2)){
        res.correction <- (GWAS.res.merge.list[[test.effect.no]])[[1]]
        res.corrections[[test.effect.no]] <- res.correction
        check.here <- (GWAS.res.merge.list[[test.effect.no]])[[2]]
        pval.correction <- res.correction[, 4]

        if (plot.qq.2) {
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
            qq(pval.correction)
            if(is.null(main.qq.2)){
              title(main = trait.name)
            }else{
              title(main = main.qq.2)
            }
          }else{
            png(paste0(saveName, trait.name, "_qq_kernel.png"))
            qq(pval.correction)
            if(is.null(main.qq.2)){
              title(main = trait.name)
            }else{
              title(main = main.qq.2)
            }
            dev.off()
          }
        }


        if (plot.Manhattan.2) {
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
              manhattan(input = res.correction, sig.level = sig.level, method.thres = method.thres, plot.col1 = plot.col1,
                        plot.type = plot.type, plot.pch = plot.pch)
              if(!is.null(plot.col3)){
                manhattan.plus(input = res.correction, checks = check.here,
                               plot.col1 = plot.col1, plot.col3 = plot.col3,
                               plot.type = plot.type, plot.pch = plot.pch)
              }
            }else{
              manhattan2(input = res.correction, sig.level = sig.level, method.thres = method.thres, plot.col2 = plot.col2,
                         plot.type = plot.type, plot.pch = plot.pch)
            }
            if(is.null(main.man.2)){
              title(main = trait.name)
            }else{
              title(main = main.man.2)
            }
          }else{
            png(paste0(saveName, trait.name, "_manhattan_kernel.png"), width = 800)
            if(plot.method == 1){
              manhattan(input = res.correction, sig.level = sig.level, method.thres = method.thres, plot.col1 = plot.col1,
                        plot.type = plot.type, plot.pch = plot.pch)
              if(!is.null(plot.col3)){
                manhattan.plus(input = res.correction, checks = check.here,
                               plot.col1 = plot.col1, plot.col3 = plot.col3,
                               plot.type = plot.type, plot.pch = plot.pch)
              }
            }else{
              manhattan2(input = res.correction, sig.level = sig.level, method.thres = method.thres, plot.col2 = plot.col2,
                         plot.type = plot.type, plot.pch = plot.pch)
            }
            if(is.null(main.man.2)){
              title(main = trait.name)
            }else{
              title(main = main.man.2)
            }
            if(!(plot.add.last & (pheno.no == n.pheno))){
              dev.off()
            }
          }
        }

        threshold <- try(CalcThreshold(GWAS.res.first.now, sig.level = sig.level, method = method.thres), silent = TRUE)
        threshold.correction <- try(CalcThreshold(res.correction, sig.level = sig.level, method = method.thres), silent = TRUE)
        if("try-error" %in% class(threshold)){
          threshold <- NA
        }
        if("try-error" %in% class(threshold.correction)){
          threshold.correction <- NA
        }
        thresholds[test.effect.no, pheno.no] <- threshold
        thresholds.correction[test.effect.no, pheno.no] <- threshold.correction
      }

      for(test.effect.no in 1:length(test.effect.2)){
        colnames(res.corrections[[test.effect.no]])[1:3] <-
          colnames(res.all[[test.effect.no]])[1:3] <- c("marker", "chrom", "pos")
        res.all[[test.effect.no]] <- merge(res.all[[test.effect.no]],
                                           res.corrections[[test.effect.no]],
                                           by.x = c("marker", "chrom", "pos"),
                                           by.y = c("marker", "chrom", "pos"),
                                           all.x = T, all.y = T)
        colnames(res.all[[test.effect.no]])[ncol(res.all[[test.effect.no]])] <-
          paste0(trait.name, "_correction")

        res.all[[test.effect.no]] <- (res.all[[test.effect.no]])[order(res.all[[test.effect.no]][, 2],
                                                                       res.all[[test.effect.no]][, 3]), ]
      }
    }else{
      colnames(RGWAS.multisnp.res) <- colnames(GWAS.res.first.now)
      if(is.null(gene.set)){
        RGWAS.multisnp.res[, 1:3] <- map[match(rownames(RGWAS.multisnp.res), rownames(map)), ]
      }
      GWAS.res.merge.0 <- rbind(RGWAS.multisnp.res, GWAS.res.first.now)
      GWAS.res.merge <- GWAS.res.merge.0[!duplicated(as.character(GWAS.res.merge.0[, 1])), ]
      ord.GWAS.res.merge <- order(GWAS.res.merge[, 2], GWAS.res.merge[, 3])
      res.correction <- GWAS.res.merge[ord.GWAS.res.merge, ]
      check.here <- match(1:nrow(RGWAS.multisnp.res), ord.GWAS.res.merge)
      pval.correction <- res.correction[, 4]



      if (plot.qq.2) {
        if (verbose) {
          print("Now Plotting (Q-Q plot). Please Wait.")
        }
        if(is.null(saveName)){
          if (length(grep("RStudio", names(dev.cur()))) == 0) {
            if (dev.cur() == dev.next()) {
              dev.new()
            } else {
              dev.set(dev.next())
            }
          }
          qq(pval.correction)
          if(is.null(main.qq.2)){
            title(main = trait.name)
          }else{
            title(main = main.qq.2)
          }
        }else{
          png(paste0(saveName, trait.name, "_qq_kernel.png"))
          qq(pval.correction)
          if(is.null(main.qq.2)){
            title(main = trait.name)
          }else{
            title(main = main.qq.2)
          }
          dev.off()
        }
      }


      if (plot.Manhattan.2) {
        if (verbose) {
          print("Now Plotting (Manhattan plot). Please Wait.")
        }
        if(is.null(saveName)){
          if (length(grep("RStudio", names(dev.cur()))) == 0) {
            if (dev.cur() == dev.next()) {
              dev.new()
            } else {
              dev.set(dev.next())
            }
          }
          if(plot.method == 1){
            manhattan(input = res.correction, sig.level = sig.level, method.thres = method.thres, plot.col1 = plot.col1,
                      plot.type = plot.type, plot.pch = plot.pch)
            if(!is.null(plot.col3)){
              manhattan.plus(input = res.correction, checks = check.here,
                             plot.col1 = plot.col1, plot.col3 = plot.col3,
                             plot.type = plot.type, plot.pch = plot.pch)
            }
          }else{
            manhattan2(input = res.correction, sig.level = sig.level, method.thres = method.thres, plot.col2 = plot.col2,
                       plot.type = plot.type, plot.pch = plot.pch)
          }
          if(is.null(main.man.2)){
            title(main = trait.name)
          }else{
            title(main = main.man.2)
          }
        }else{
          png(paste0(saveName, trait.name, "_manhattan_kernel.png"), width = 800)
          if(plot.method == 1){
            manhattan(input = res.correction, sig.level = sig.level, method.thres = method.thres, plot.col1 = plot.col1,
                      plot.type = plot.type, plot.pch = plot.pch)
            if(!is.null(plot.col3)){
              manhattan.plus(input = res.correction, checks = check.here,
                             plot.col1 = plot.col1, plot.col3 = plot.col3,
                             plot.type = plot.type, plot.pch = plot.pch)
            }
          }else{
            manhattan2(input = res.correction, sig.level = sig.level, method.thres = method.thres, plot.col2 = plot.col2,
                       plot.type = plot.type, plot.pch = plot.pch)
          }
          if(is.null(main.man.2)){
            title(main = trait.name)
          }else{
            title(main = main.man.2)
          }
          if(!(plot.add.last & (pheno.no == n.pheno))){
            dev.off()
          }
        }
      }

      threshold <- try(CalcThreshold(GWAS.res.first[, c(1:3, pheno.no + 3)], sig.level = sig.level, method = method.thres), silent = TRUE)
      threshold.correction <- try(CalcThreshold(res.correction, sig.level = sig.level, method= method.thres), silent = TRUE)
      if("try-error" %in% class(threshold)){
        threshold <- NA
      }
      if("try-error" %in% class(threshold.correction)){
        threshold.correction <- NA
      }
      thresholds[, pheno.no] <- threshold
      thresholds.correction[, pheno.no] <- threshold.correction
      colnames(res.correction)[1:3] <- colnames(res.all)[1:3] <- c("marker", "chrom", "pos")
      res.all <- merge(res.all, res.correction, by.x = c("marker", "chrom", "pos"),
                       by.y = c("marker", "chrom", "pos"), all.x = T, all.y = T)
      colnames(res.all)[ncol(res.all)] <- paste0(trait.name, "_correction")
      res.all <- res.all[order(res.all[, 2], res.all[, 3]), ]
    }
  }


  thresholds.list <- list(first = thresholds, second = thresholds.correction)

  if(thres){
    end <- Sys.time()

    if(time){
      print(end - start)
    }

    if(return.EMM.res){
      return(list(D = res.all, thres = thresholds.list,
                  EMM.res = EMM.res0))
    }else{
      return(list(D = res.all, thres = thresholds.list))
    }
  }else{
    end <- Sys.time()

    if(time){
      print(end - start)
    }

    if(return.EMM.res){
      return(list(D = res.all, EMM.res = EMM.res0))
    }else{
      return(res.all)
    }
  }
}
