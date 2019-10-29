#'  Perform normal GWAS (genome-wide association studies) first, then check epistatic effects for relatively significant markers
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
#' @param check.size.epi This argument determines how many SNPs (around the SNP detected by normal GWAS) you will check epistasis.
#' @param epistasis.percent This argument determines how many SNPs are detected by normal GWAS.
#' For example, when epistasis.percent = 0.1, SNPs whose value of -log10(p) is in the top 0.1 percent are chosen as candidate for checking epistasis.
#' @param check.epi.max It takes a lot of time to check epistasis, so you can decide the maximum number of SNPs to check epistasis.
#' @param your.check  Because there are less SNPs that can be tested in epistasis than in kernel-based GWAS, you can select which SNPs you want to test.
#'  If you use this argument, please set the number where SNPs to be tested are located in your data (so not position).
#'  In the default setting, your_check = NULL and epistasis between SNPs detected by GWAS will be tested.
#' @param GWAS.res.first If you have already performed regular GWAS and have the result, you can skip performing normal GWAS.
#' @param P3D When P3D = TRUE, variance components are estimated by REML only once, without any markers in the model.
#' When P3D = FALSE, variance components are estimated by REML for each marker separately.
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
#' @param optimizer The function used in the optimization process. We offer "optim", "optimx", and "nlminb" functions.
#' @param chi0.mixture RAINBOWR assumes the deviance is considered to follow a x chisq(df = 0) + (1 - a) x chisq(df = r).
#' where r is the degree of freedom.
#' The argument chi0.mixture is a (0 <= a < 1), and default is 0.5.
#' @param gene.set If you have information of gene (or haplotype block), you can use it to perform kernel-based GWAS.
#'            You should assign your gene information to gene.set in the form of a "data.frame" (whose dimension is (the number of gene) x 2).
#'            In the first column, you should assign the gene name. And in the second column, you should assign the names of each marker,
#'            which correspond to the marker names of "geno" argument.
#' @param sig.level Significance level for the threshold. The default is 0.05.
#' @param method.thres Method for detemining threshold of significance. "BH" and "Bonferroni are offered.
#' @param plot.qq.1 If TRUE, draw qq plot for normal GWAS.
#' @param plot.Manhattan.1 If TRUE, draw manhattan plot for normal GWAS.
#' @param plot.epi.3d If TRUE, draw 3d plot
#' @param plot.epi.2d If TRUE, draw 2d plot
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
#' @param main.qq.1 The title of qq plot for normal GWAS. If this argument is NULL, trait name is set as the title.
#' @param main.man.1 The title of manhattan plot for normal GWAS. If this argument is NULL, trait name is set as the title.
#' @param main.epi.3d The title of 3d plot. If this argument is NULL, trait name is set as the title.
#' @param main.epi.2d The title of 2d plot. If this argument is NULL, trait name is set as the title.
#' @param verbose If this argument is TRUE, messages for the current steps will be shown.
#' @param verbose2 If this argument is TRUE, welcome message will be shown.
#' @param count When count is TRUE, you can know how far RGWAS has ended with percent display.
#' @param time When time is TRUE, you can know how much time it took to perform RGWAS.
#'
#'
#' @return
#' \describe{
#' \item{$first}{The results of first normal GWAS will be returned.}
#' \item{$epistasis}{
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
#' }
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
#'
#' @example R/examples/RGWAS.twostep.epi_example.R
#'
#'
RGWAS.twostep.epi <- function(pheno, geno, ZETA = NULL, covariate = NULL, covariate.factor = NULL,
                              structure.matrix = NULL, n.PC = 0, min.MAF = 0.02, n.core = 1,
                              check.size.epi = 4, epistasis.percent = 0.05, check.epi.max = 200,
                              your.check = NULL, GWAS.res.first = NULL, P3D = TRUE, test.method = "LR",
                              dominance.eff = TRUE, haplotype = TRUE, num.hap = NULL, optimizer = "nlminb",
                              window.size.half = 5, window.slide = 1, chi0.mixture = 0.5,
                              gene.set = NULL, sig.level = 0.05, method.thres = "BH", plot.qq.1 = TRUE, plot.Manhattan.1 = TRUE,
                              plot.epi.3d = TRUE, plot.epi.2d = TRUE, plot.method = 1,
                              plot.col1 = c("dark blue", "cornflowerblue"), plot.col2 = 1,
                              plot.type = "p", plot.pch = 16, saveName = NULL, main.qq.1 = NULL,
                              main.man.1 = NULL, main.epi.3d = NULL, main.epi.2d = NULL,
                              verbose = TRUE, verbose2 = FALSE, count = TRUE, time = TRUE){

  start <- Sys.time()
  if(is.null(GWAS.res.first)){
    if (verbose) {
      print("The 1st step: Performing normal GWAS!")
    }
    GWAS.res.first <- RGWAS.normal(pheno = pheno, geno = geno, ZETA = ZETA, covariate = covariate,
                                   covariate.factor = covariate.factor, structure.matrix = structure.matrix,
                                   n.PC = n.PC, min.MAF = min.MAF, P3D = P3D, n.core = n.core,
                                   sig.level = sig.level, method.thres = method.thres, plot.qq = plot.qq.1, plot.Manhattan = plot.Manhattan.1,
                                   plot.method = plot.method, plot.col1 = plot.col1, plot.col2 = plot.col2,
                                   plot.type = plot.type, plot.pch = plot.pch, saveName = saveName, optimizer = optimizer,
                                   main.qq = main.qq.1, main.man = main.man.1, plot.add.last = FALSE, return.EMM.res = FALSE,
                                   thres = FALSE, verbose = verbose, verbose2 = verbose2, count = count, time = time)
  }else{
    if (verbose){
        print("The 1st step has already finished because you input 'GWAS.res.first'.")
      }
  }


  n.pheno <- ncol(GWAS.res.first) - 3
  trait.names <- colnames(GWAS.res.first)[4:(4 + n.pheno - 1)]
  map <- geno[, 1:3]
  chr <- map[, 2]
  chr.tab <- table(chr)
  chr.max <- max(chr)
  chr.cum <- cumsum(chr.tab)
  pos <- map[,3]
  cum.pos <- pos
  if(length(chr.tab) != 1){
    for(i in 1:(chr.max - 1)){
      cum.pos[(chr.cum[i] + 1): (chr.cum[i + 1])] <- pos[(chr.cum[i] + 1):(chr.cum[i + 1])] + cum.pos[chr.cum[i]]
    }
  }


  if(n.pheno != 1){
    all.epi.res <- rep(list(NA), n.pheno)
  }else{
    all.epi.res <- NULL
  }



  for(pheno.no in 1:n.pheno){
    trait.name <- trait.names[pheno.no]
    pheno.now <- pheno[, c(1, pheno.no + 1)]

    if(is.null(your.check)){
      pval.first <- GWAS.res.first[, pheno.no + 3]
      ord.pval.first <- order(pval.first, decreasing = TRUE)
      check.epi.no.0 <- round(length(pval.first) * (epistasis.percent / 100), 0) * (check.size.epi + 1)
      check.epi.no <- ifelse(check.epi.no.0 >= check.epi.max, check.epi.max, check.epi.no.0)
      ord.pval.epi.percent <- ord.pval.first[1:check.epi.no]

      check.size.epi.half <- check.size.epi / 2
      checks.mat <- matrix(NA, nrow = length(ord.pval.epi.percent) * (check.size.epi + 1),
                           ncol = length(ord.pval.epi.percent))
      for(check.no in 1:length(ord.pval.epi.percent)){
        check <- sort(ord.pval.epi.percent)[check.no]
        checks.mat[, check.no] <- (check - check.size.epi.half):(check + check.size.epi.half)
      }
      checks <- unique(c(checks.mat))
      checks <- checks[(checks >= 1) & (checks <= length(pos))]
    }else{
      checks <- your.check
      checks <- checks[(checks >= 1) & (checks <= length(pos))]
    }

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
    M.check <- geno[checks, -c(1:3)]
    geno.check <- cbind(pseudo.map, M.check)

    if (verbose) {
      print(paste("The 2nd step: Calculating -log10(p) of epistatic effects of", trait.name, "for", n.checks, "x", n.checks,"SNPs."))
    }
    RGWAS.epistasis.res <- RGWAS.epistasis(pheno = pheno.now, geno = geno.check, ZETA = ZETA, covariate = covariate,
                                           covariate.factor = covariate.factor, structure.matrix = structure.matrix,
                                           n.PC = n.PC, min.MAF = min.MAF, n.core = n.core,
                                           test.method = test.method, dominance.eff = dominance.eff, haplotype = haplotype,
                                           num.hap = num.hap, window.size.half = window.size.half, window.slide = window.slide,
                                           chi0.mixture = chi0.mixture, gene.set = gene.set, optimizer = optimizer,
                                           plot.epi.3d = FALSE, plot.epi.2d = FALSE, main.epi.3d = main.epi.3d,
                                           main.epi.2d = main.epi.2d, saveName = saveName, verbose = verbose,
                                           verbose2 = verbose2, count = count, time = time)


    check.tests <- as.numeric(rownames(RGWAS.epistasis.res$map))
    n.check.tests <- length(check.tests)

    scores.epi <- RGWAS.epistasis.res$scores$scores
    rownames(scores.epi) <- colnames(scores.epi) <- check.tests

    cum.pos2 <- cum.pos[check.tests]
    x.3d <- rep(cum.pos2, n.check.tests)
    y.3d <- rep(cum.pos2, each = n.check.tests)
    z.3d <- c(scores.epi)

    epi.res <- list(scores = scores.epi, x = x.3d, y = y.3d, z = z.3d)

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
               plot.epi.2d = plot.epi.2d, main.epi.3d = main.epi.3d,
               main.epi.2d = main.epi.2d, saveName = saveName)

    if(n.pheno >= 2){
      all.epi.res[[pheno.no]] <- epi.res
      names(all.epi.res) <- trait.names
    }else{
      all.epi.res <- epi.res
    }
  }



  end <- Sys.time()

  if(time){
    print(end - start)
  }

  return(list(first = GWAS.res.first, epistasis = all.epi.res))
}
