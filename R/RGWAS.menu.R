#' Print the R code which you should perform for RAINBOWR GWAS
#' @description Print the R code which you should perform for RAINBOWR (Reliable Association INference By Optimizing Weights with R).
#' @return The R code which you should perform for RAINBOWR GWAS
#'
#'
RGWAS.menu <- function(){
  test.methods <- c("'LR'", "'score'")
  gene.sets <- c("gene.set", "NULL")
  kernel.methods <- c("'linear'", "'gaussian'", "'exponential'")
  test.effects <- c("'additive'", "'dominance'", "'additive+dominance'")
  weighting.centers <- c("TRUE", "FALSE")
  dominance.effs <- c("FALSE", "TRUE")

  Q1 <- menu(c("Normal GWAS; testing each single SNP",
               "Kernel-based GWAS; testing multiple SNPs simultaneously",
               "Check epistatic effects"),
             title = "Which GWAS do you want to perform?")

  if(Q1 != 1){
    Q2 <- menu(c("One step method; perform kerenel method for all SNPs",
                 "Two step method; perform kerenel method for SNPs which were detected by normal GWAS"),
               title = "Which method do you want to perform?")

    Q3 <- menu(c("LR; likelihood-ratio test, relatively slow, but accurate (default)",
                 "score; score test, much faster than LR, but sometimes overestimate -log10(p)"),
               title = "Kernel-based GWAS have two methods to test SNPs. Which do you want to use?")

    Q4 <- menu(c("Yes; when you input your gene information, see Arguments of 'gene.set' for details of its format",
                 "No; perform tests for each SNP-set containing a fixed number of SNPs"),
               title = "Do you have some gene information? RGWAS can perform gene-set based GWAS.")

    Q5 <- menu(c("Yes; RGWAS will only test for additive effects",
                 "No; you can choose which effects to test"),
               title = "Do you use inbred lines?")

    if(Q1 == 2){
      Q6 <- menu(c("linear; fastest especially when the number of SNPs in SNP-set is small (default)",
                   "gaussian; slower than linear (especially with LR test), but it may catch non-linear effects in SNP-set",
                   "exponential; slower than linear (especially with LR test), but it may have high detection sensitivity"),
                 title = "Which kernel do you want to use?")

      if(Q5 == 2){
      Q7 <- menu(c("additive; test additive effects (default)",
                   "dominance; test dominance effects, however score test sometimes overestimates this effect",
                   "additive+dominance; test both additive and dominance effects simultaneouly"),
                 title = "Which genetic effects do you want to test?")
      }else{
        Q7 <- 1
      }

      Q8 <- menu(c("Yes; weighting is recommended when you perform kernel-based GWAS by sliding windows (default)",
                   "No; weighting is not recommended when you perform kernel based GWAS by separated bins"),
                 title = "RGWAS can weight the SNP at the center of the SNP-set. Do you want to weight?")
    }
  }



  if(Q1 == 1){
    func.name <- "RGWAS.normal"
    output <- paste0(func.name, "(pheno, geno, ZETA = ZETA, covariate =  NULL, covariate.factor = NULL, ",
                     "structure.matrix = NULL, n.PC = 0, min.MAF = 0.02, P3D = TRUE, n.core = 1, ",
                     "sig.level = 0.05, plot.qq = TRUE, plot.Manhattan = TRUE, plot.method = 1, ",
                     "plot.col1 = c('dark blue', 'cornflowerblue'), plot.col2 = 1, ",
                     "plot.type = 'p', plot.pch = 16, saveName = NULL, main.qq = NULL, ",
                     "main.man = NULL, plot.add.last = FALSE, return.EMM.res = FALSE, ",
                     "thres = TRUE, verbose = FALSE, count = TRUE, time = TRUE)")
  }else{
    if(Q1 == 2){
      if(Q2 == 1){
        func.name <- "RGWAS.multisnp"
        output <- paste0(func.name, "(pheno, geno, ZETA = ZETA, covariate =  NULL, covariate.factor = NULL, ",
                         "structure.matrix = NULL, n.PC = 0, min.MAF = 0.02, test.method = ", test.methods[Q3],
                         ", n.core = 1, kernel.method = ", kernel.methods[Q6], ", kernel.h = 'tuned', ",
                         "haplotype = TRUE, num.hap = NULL, test.effect = ", test.effects[Q7],
                         ", window.size.half = 5, window.slide = 1, chi0.mixture = 0.5, ", "gene.set = ",
                         gene.sets[Q4], ", weighting.center = ", weighting.centers[Q8], ", weighting.other = NULL, ",
                         "sig.level = 0.05, plot.qq = TRUE, plot.Manhattan = TRUE, plot.method = 1, ",
                         "plot.col1 = c('dark blue', 'cornflowerblue'), plot.col2 = 1, ",
                         "plot.type = 'p', plot.pch = 16, saveName = NULL, main.qq = NULL, ",
                         "main.man = NULL, plot.add.last = FALSE, return.EMM.res = FALSE, ",
                         "thres = TRUE, verbose = FALSE, count = TRUE, time = TRUE)")
      }else{
        func.name <- "RGWAS.twostep"
        output <- paste0(func.name, "(pheno, geno, ZETA = ZETA, covariate =  NULL, covariate.factor = NULL, ",
                         "structure.matrix = NULL, n.PC = 0, min.MAF = 0.02, n.core = 1, ",
                         "check.size = 40, check.gene.size = 4, kernel.percent = 0.1, GWAS.res.first = NULL, ",
                         "P3D = TRUE, test.method.1 = 'normal', test.method.2 = ", test.methods[Q3],
                         ", kernel.method = ", kernel.methods[Q6], ", kernel.h = 'tuned', ",
                         "haplotype = TRUE, num.hap = NULL, test.effect.1 = ", test.effects[Q7],
                         ", test.effect.2 = ", test.effects[Q7],
                         ", window.size.half = 5, window.slide = 1, chi0.mixture = 0.5, ", "gene.set = ",
                         gene.sets[Q4], ", weighting.center = ", weighting.centers[Q8], ", weighting.other = NULL, ",
                         "sig.level = 0.05, plot.qq.1 = TRUE, plot.Manhattan.1 = TRUE, ",
                         "plot.qq.2 = TRUE, plot.Manhattan.2 = TRUE, plot.method = 1, ",
                         "plot.col1 = c('dark blue', 'cornflowerblue'), plot.col2 = 1, ",
                         "plot.col3 = c('red3', 'orange3'), plot.type = 'p', plot.pch = 16, saveName = NULL, ",
                         "main.qq.1 = NULL, main.man.1 = NULL, main.qq.2 = NULL, main.man.2 = NULL, ",
                         "plot.add.last = FALSE, return.EMM.res = FALSE, ",
                         "thres = TRUE, verbose = FALSE, count = TRUE, time = TRUE)")
      }
    }else{
      if(Q2 == 1){
        func.name <- "RGWAS.epistasis"
        output <- paste0(func.name, "(pheno, geno, ZETA = ZETA, covariate =  NULL, covariate.factor = NULL, ",
                         "structure.matrix = NULL, n.PC = 0, min.MAF = 0.02, n.core = 1, test.method = ",
                         test.methods[Q3], ", haplotype = TRUE, dominance.eff = ",
                         dominance.effs[Q5], ", num.hap = NULL, window.size.half = 5, window.slide = 1, ",
                         "chi0.mixture = 0.5, ", "gene.set = ", gene.sets[Q4],
                         ", plot.epi.3d = TRUE, plot.epi.2d = TRUE, main.epi.3d = NULL, ",
                         "main.epi.2d = NULL, saveName = NULL, ",
                         "verbose = FALSE, count = TRUE, time = TRUE)")
      }else{
        func.name <- "RGWAS.twostep.epi"
        output <- paste0(func.name, "(pheno, geno, ZETA = ZETA, covariate =  NULL, covariate.factor = NULL, ",
                         "structure.matrix = NULL, n.PC = 0, min.MAF = 0.02, n.core = 1, ",
                         "check.size.epi = 4, epistasis.percent = 0.1, check.epi.max = 200, ",
                         "your.check = NULL, GWAS.res.first = NULL, ",
                         "P3D = TRUE, test.method = ", test.methods[Q3], ", haplotype = TRUE, ",
                         "dominance.eff = ", dominance.effs[Q5], ", num.hap = NULL, window.size.half = 5, ",
                         "window.slide = 1, chi0.mixture = 0.5, gene.set = ", gene.sets[Q4],
                         ", sig.level = 0.05, plot.qq.1 = TRUE, plot.Manhattan.1 = TRUE, ",
                         "plot.epi.3d = TRUE, plot.epi.2d = TRUE, plot.method = 1, ",
                         "plot.col1 = c('dark blue', 'cornflowerblue'), plot.col2 = 1, ",
                         "plot.type = 'p', plot.pch = 16, saveName = NULL, ",
                         "main.qq.1 = NULL, main.man.1 = NULL, main.epi.3d = NULL, main.epi.2d = NULL, ",
                         "verbose = FALSE, count = TRUE, time = TRUE)")
      }
    }
  }


  return(output)
}
