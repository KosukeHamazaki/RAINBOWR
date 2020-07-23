#' Function to greet to users
#'
#'
#' @return Show welcome messages
#'
#'
#'
# welcome_to_RGWAS <- function(){
#   cat("#-------------------------------------------------------------------------------------------------# \n")
#   cat("#          ^   ^   - - - -\\       Welcome to RAINBOWR GWAS!!!                                      #\n")
#   cat("#        ( -   - )      (@_/      Reliable Association INference By Optimizing the Weight         #\n")
#   cat("#       //   -   \\\\  - - -        Function for performing normal and kernel-based GWAS.           #\n")
#   cat("#-------------------------------------------------------------------------------------------------# \n")
# }
welcome_to_RGWAS <- function(){
  if (getPackageName() == ".GlobalEnv") {
    version <- "devel"
  }
  else {
    version <- as.character(packageVersion(getPackageName()))
  }
  
  cat("#------------------------ Reliable Association INference By Optimizing Weights -------------------------#\n")
  cat("#  _____         --      _____  __      _  ____     __    _     _     _                                 #\n")
  cat("#  |  __ \\      /  \\    |_   _||  \\    | ||  _ \\ /  __ \\ | |   | |   | |  Welcome to RAINBOW GWAS!!!    #\n")
  cat("#  | |__) |    / __ \\     | |  | |\\ \\  | || |_) || |  | || |   | |   | |  Version:",
      version, "               #\n")
  cat("#  |  ___/    / /__\\ \\    | |  | | \\ \\ | ||  _ <|| |  | | \\ \\  | |  / /                                 #\n")
  cat("#  | |  \\ \\  / ______ \\  _| |_ | |  \\ \\| || |_) || |__| |  | |_| |_| |    email: hamazaki@ut-biomet.org #\n")
  cat("#  |_|   \\_\\/ /      \\ \\|_____||_|   \\ __||____/  \\ __ /    \\ _ _ _ /                                   #\n")
  cat("#--------------------- Function for performing normal and kernel-based SNP-set GWAS --------------------#\n")
}




#' Function to calculate threshold for GWAS
#'
#' @description Calculate thresholds for the given GWAS (genome-wide association studies) result by the Benjamini-Hochberg method or Bonferroni method.
#'
#' @param input Data frame of GWAS results where the first column is the marker names,
#' the second and third column is the chromosome amd map position, and the forth column is -log10(p) for each marker.
#' @param sig.level  Significance level for the threshold. The default is 0.05. You can also assign vector of sinificance levels.
#' @param method Two methods are offered:
#' 
#' "BH" : Benjamini-Hochberg method. To control FDR, use this method.
#' "Bonf" : Bonferroni method. To perform simple correction of multiple testing, use this method.
#' 
#' You can also assign both of them by 'method = c("BH", "Bonf")'
#'
#' @return The value of the threshold. If there is no threshold, it returns NA.
#'
#' @references Benjamini, Y. and Hochberg, Y. (1995) Controlling the false discovery rate:
#'  a practical and powerful approach to multiple testing. J R Stat Soc. 57(1): 289-300.
#'
#' Storey, J.D. and Tibshirani, R. (2003) Statistical significance for genomewide studies. Proc Natl Acad Sci. 100(16): 9440-9445.
#'
#'
#'
CalcThreshold <- function(input, sig.level = 0.05, method = "BH") {
  # define a function
  qvalue_tmp <- function(p) {
    smooth.df <- 3
    
    if(min(p) < 0 || max(p) > 1) {
      stop("P-values not in valid range.")
      return(0)
    }
    
    lambda <- seq(0, 0.90, 0.05)
    m <- length(p)
    
    pi0 <- rep(0, length(lambda))
    for(i in 1:length(lambda)) {
      pi0[i] <- mean(p >= lambda[i]) / (1 - lambda[i])
    }
    
    spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
    pi0 <- predict(spi0, x = max(lambda))$y
    pi0 <- min(pi0, 1)
    
    if(pi0 <= 0) {
      stop("The estimated pi0 <= 0. Check that you have valid p-values.")
      return(0)
    }
    
    #The estimated q-values calculated here
    u <- order(p)
    
    # ranking function which returns number of observations less than or equal
    qvalue.rank <- function(x) {
      idx <- sort.list(x)
      
      fc <- factor(x)
      nl <- length(levels(fc))
      bin <- as.integer(fc)
      tbl <- tabulate(bin)
      cs <- cumsum(tbl)
      
      tbl <- rep(cs, tbl)
      tbl[idx] <- tbl
      
      return(tbl)
    }
    
    v <- qvalue.rank(p)
    
    qvalue <- pi0 * m * p / v
    qvalue[u[m]] <- min(qvalue[u[m]], 1)
    for(i in (m-1):1) {
      qvalue[u[i]] <- min(qvalue[u[i]], qvalue[u[i + 1]], 1)
    }
    
    return(qvalue)
  }
  
  input <- input[!is.na(input[, 4]), , drop = FALSE]
  input <- input[order(input[, 2], input[, 3]), ]
  
  method[!(method %in% c("BH", "Bonf"))] <- "BH"
  methods <- rep(method, each = length(sig.level))
  sig.levels <- rep(sig.level, length(method))
  
  n.thres <- length(methods)
  
  thresholds <- rep(NA, n.thres)
  for(thres.no in 1:n.thres){
    method.now <- methods[thres.no]
    sig.level.now <- sig.levels[thres.no]
    
    if(method.now == "BH"){
      # input should be a result object of GWAS in {rrBLUP} package
      q.ans <- qvalue_tmp(10 ^ (- input[, 4]))
      temp <- cbind(q.ans, input[, 4])
      temp <- temp[order(temp[, 1]), ]
      if (temp[1, 1] < sig.level.now) {
        temp2 <- tapply(temp[, 2], temp[, 1], mean)
        qvals <- as.numeric(rownames(temp2))
        x <- which.min(abs(qvals - sig.level.now))
        first <- max(1, x - 2)
        last <- min(x + 2, length(qvals))
        if ((last - first) < 4) {
          last <- first + 3
        }
        
        if(sum(is.na(qvals[first:last])) == 1){
          qvals[last] <- mean(qvals[first + 1] + qvals[first + 2])
          temp2[last] <- mean(temp2[first + 1] + temp2[first + 2])
        }
        
        if(sum(is.na(qvals[first:last])) == 2){
          qvals[(last - 1):last] <- quantile(qvals[first:(first + 1)], probs = c(1 / 3, 2 / 3))
          temp2[(last - 1):last] <- quantile(temp2[first:(first + 1)], probs = c(1 / 3, 2 / 3))
        }
        
        qvals <- sort(qvals)
        temp2 <- temp2[order(qvals)]
        
        splin <- smooth.spline(x = qvals[first:last], y=temp2[first:last], df = 3)
        threshold <- predict(splin, x = sig.level.now)$y
      }else{
        threshold <- NA
      }
    }
    
    if(method.now == "Bonf"){
      n.mark <- nrow(input)
      threshold <- -log10(sig.level.now / n.mark)
    }
    
    thresholds[thres.no] <- threshold
  }
  names(thresholds) <- paste0(methods, "_", sig.levels)
  return(thresholds)
}


#' Function to calculate genomic relationship matrix (GRM)
#'
#' @param genoMat A \eqn{N \times M} matrix of marker genotype
#' @param methodGRM Method to calculate genomic relationship matrix (GRM). We offer the following methods;
#' "addNOIA", "domNOIA", "A.mat", "linear", "gaussian", "exponential", "correlation".
#' For NOIA methods, please refer to Vitezica et al. 2017.
#' @param kernel.h The hyper parameter for gaussian or exponential kernel.
#' If kernel.h = "tuned", this hyper parameter is calculated as the median of off-diagonals of distance matrix of genotype data.
#' @param returnWMat  If this argument is TRUE, we will return W matrix instead of GRM.
#' Here, W satisfies \eqn{GRM = W W ^ {T}}. W corresponds to H matix in Vitezica et al. 2017.
#' @param probaa Probability of being homozygous for the reference allele for each marker.
#' If NULL (default), it will be calculated from genoMat.
#' @param probAa Probability of being heterozygous for the reference and alternative alleles for each marker
#' If NULL (default), it will be calculated from genoMat.
#' @return genomic relationship matrix (GRM)
#'
#' @references 
#' 
#' Vitezica, Z.G., Legarra, A., Toro, M.A. and Varona, L. (2017) Orthogonal Estimates of Variances for Additive, Dominance, and Epistatic Effects in Populations. Genetics. 206(3): 1297-1307.
#' 
#' Endelman, J.B. and Jannink, J.L. (2012) Shrinkage Estimation of the Realized Relationship Matrix. G3 Genes, Genomes, Genet. 2(11): 1405-1413.
#'
calcGRM = function(genoMat,
                   methodGRM = "addNOIA",
                   kernel.h = "tuned",
                   returnWMat = FALSE,
                   probaa = NULL,
                   probAa = NULL) {
  supportedMethods <- c("addNOIA", "domNOIA", "A.mat", "linear",
                        "gaussian", "exponential", "correlation")
  stopifnot(methodGRM %in% supportedMethods)
  
  nInd <- nrow(genoMat)
  nMarkers <- ncol(genoMat)
  mrkNames <- colnames(genoMat)
  
  methodNOIA <- stringr::str_detect(string = methodGRM,
                                    pattern = "NOIA")
  if (methodNOIA) {
    if (is.null(probaa)) {
      probaa <- apply(genoMat == -1, 2, mean)
    }
    if (is.null(probAa)) {
      probAa <- apply(genoMat == 0, 2, mean)
    }
    if (methodGRM == "addNOIA") {
      replaceaa <- - (2 - probAa - 2 * probaa)
      replaceAa <- - (1 - probAa - 2 * probaa)
      replaceAA <- - (- probAa - 2 * probaa)
    } else if (methodGRM == "domNOIA") {
      probAA <- 1 - probaa - probAa
      denominator <- probAA + probaa - (probAA - probaa) ^ 2
      replaceaa <- - 2 * probAA * probAa / denominator
      replaceAa <- 4 * probAA * probaa /denominator
      replaceAA <- - 2 * probaa * probAa /denominator
    }
    
    HMat <- sapply(1:nMarkers, function(mrkNo) {
      HMatEachMrk <- genoMat[, mrkNo]
      HMatEachMrk[HMatEachMrk == -1] <- replaceaa[mrkNo]
      HMatEachMrk[HMatEachMrk == 0] <- replaceAa[mrkNo]
      HMatEachMrk[HMatEachMrk == 1] <- replaceAA[mrkNo]
      
      return(HMatEachMrk)
    })
    colnames(HMat) <- mrkNames
    
    HHt <- tcrossprod(HMat)
    GRM <- HHt * nInd / sum(diag(HHt))
  } else if (methodGRM == "A.mat") {
    GRM <- rrBLUP::A.mat(X = genoMat)
  } else if (methodGRM == "linear") {
    HHt <- tcrossprod(genoMat)
    GRM <- HHt * nInd / sum(diag(HHt))
  } else if (methodGRM == "gaussian") {
    distMat <- as.matrix(dist(genoMat)) / sqrt(ncol(genoMat))
    if ("character" %in% class(kernel.h)) {
      hinv <- median((distMat ^ 2)[upper.tri(distMat ^ 2)])
      h <- 1 / hinv
    } else if ("numeric" %in% class(kernel.h)) {
      h <- kernel.h
    } 
    
    GRM <- exp(- h * distMat ^ 2)
  } else if (methodGRM == "exponential") {
    distMat <- as.matrix(dist(genoMat)) / sqrt(ncol(genoMat))
    if ("character" %in% class(kernel.h)) {
      hinv <- median((distMat ^ 2)[upper.tri(distMat ^ 2)])
      h <- 1 / hinv
    } else if ("numeric" %in% class(kernel.h)) {
      h <- kernel.h
    } 
    
    GRM <- exp(- h * distMat)
  } else if (methodGRM == "correlation") {
    GRM <- cor(t(genoMat))
  }
  
  
  if (methodNOIA & returnWMat) {
    WMat <- HMat * sqrt(nInd / sum(HMat * HMat))
    return(WMat)
  } else {
    return(GRM)
  }
}

#' Function to generate design matrix (Z)
#'
#'
#' @param pheno.labels A vector of genotype (line; accesion; variety) names which correpond to phenotypic values.
#' @param geno.names  A vector of genotype (line; accesion; variety) names for marker genotype data (duplication is not recommended).
#'
#' @return Z of \eqn{y = X\beta + Zu + e}. Design matrix, which is useful for GS or GWAS.
#'
#'
#'
design.Z <- function(pheno.labels, geno.names) {
  pheno.labels <- as.character(pheno.labels)
  geno.names <- as.character(geno.names)
  n.geno <- length(geno.names)
  
  match.pheno_geno <- match(pheno.labels, geno.names)
  
  if(any(is.na(match.pheno_geno))){
    warning(paste("The following lines have phenotypes but no genotypes: ",
                  paste(pheno.labels[is.na(match.pheno_geno)], collapse = ", ")))
  }
  
  match.pheno_geno.nonNA <- match.pheno_geno[!is.na(match.pheno_geno)]
  n.pheno.nonNA <- length(match.pheno_geno.nonNA)
  
  Z <- as.matrix(Matrix::sparseMatrix(i = 1:n.pheno.nonNA,
                                      j = match.pheno_geno.nonNA,
                                      x = rep(1, n.pheno.nonNA),
                                      dims = c(n.pheno.nonNA, n.geno)))
  rownames(Z) <- pheno.labels[!is.na(match.pheno_geno)]
  colnames(Z) <- geno.names
  
  return(Z)
}





#' Function to modify genotype and phenotype data to match
#'
#'
#' @param pheno.mat A \eqn{n _ 1 \times p} matrix of phenotype data. rownames(pheno.mat) should be genotype (line; accesion; variety) names.
#' @param geno.mat A \eqn{n _ 2 \times m} matrix of marker genotype data. rownames(geno.mat) should be genotype (line; accesion; variety) names.
#' @param pheno.labels A vector of genotype (line; accesion; variety) names which correpond to phenotypic values.
#' @param geno.names  A vector of genotype (line; accesion; variety) names for marker genotype data (duplication is not recommended).
#' @param map Data frame with the marker names in the first column. The second and third columns contain the chromosome and map position.
#' @param return.ZETA If this argument is TRUE, the list for mixed model equation (ZETA) will be returned.
#' @param return.GWAS.format If this argument is TRUE, phenotype and genotype data for GWAS will be returned.
#' @return
#' \describe{
#' \item{$geno.modi}{The modified marker genotype data.}
#' \item{$pheno.modi}{The modified phenotype data.}
#' \item{$ZETA}{The list for mixed model equation (ZETA).}
#' \item{$pheno.GWAS}{GWAS formatted phenotype data.}
#' \item{$geno.GWAS}{GWAS formatted marker genotype data.}
#' }
#'
#'
#'
modify.data <- function(pheno.mat, geno.mat, pheno.labels = NULL, geno.names = NULL, map = NULL,
                        return.ZETA = TRUE, return.GWAS.format = FALSE) {
  pheno.mat <- as.matrix(pheno.mat)
  
  if(is.null(pheno.labels)){
    pheno.labels <- as.character(rownames(pheno.mat))
  }else{
    pheno.labels <- as.character(pheno.labels)
    rownames(pheno.mat) <- pheno.labels
  }
  
  if(is.null(geno.names)){
    geno.names <- as.character(rownames(geno.mat))
  }else{
    geno.names <- as.character(geno.names)
    rownames(geno.mat) <- geno.names
  }
  both.names <- Reduce(intersect, list(pheno.labels, geno.names))
  
  
  match.pheno <- match(pheno.labels, both.names)
  match.geno <- match(geno.names, both.names)
  
  
  if(any(is.na(match.pheno))){
    warning(paste("The following lines have phenotypes but no genotypes: ",
                  paste(pheno.labels[is.na(match.pheno)], collapse = ", ")))
  }
  
  
  pheno.mat.match <- pheno.mat[!is.na(match.pheno), , drop = FALSE]
  pheno.mat.modi <- pheno.mat.match[order(match.pheno[!is.na(match.pheno)]), , drop = FALSE]
  pheno.labels.modi <- rownames(pheno.mat.modi)
  
  geno.mat.match <- geno.mat[!is.na(match.geno), , drop = FALSE]
  geno.mat.modi <- geno.mat.match[order(match.geno[!is.na(match.geno)]), , drop = FALSE]
  geno.names.modi <- rownames(geno.mat.modi)
  
  if(return.ZETA){
    K.A <- calcGRM(genoMat = geno.mat.modi)
    Z.A <- design.Z(pheno.labels = pheno.labels.modi, geno.names = geno.names.modi)
    
    ZETA <- list(A = list(Z = Z.A, K = K.A))
  }else{
    ZETA <- NULL
  }
  
  if(return.GWAS.format & (!is.null(map))){
    pheno.GWAS <- data.frame(Sample_names = pheno.labels.modi, pheno.mat.modi)
    
    geno.GWAS <- data.frame(map, t(geno.mat.modi))
    rownames(geno.GWAS) <- 1:ncol(geno.mat.modi)
    colnames(geno.GWAS) <- c("marker", "chrom", "pos", geno.names.modi)
  }else{
    pheno.GWAS <- geno.GWAS <- NULL
  }
  
  return(list(geno.modi = geno.mat.modi, pheno.modi = pheno.mat.modi,
              ZETA = ZETA, pheno.GWAS = pheno.GWAS, geno.GWAS = geno.GWAS))
}






#' Function to calculate cumulative position (beyond chromosome)
#'
#'
#' @param map Data frame with the marker names in the first column. The second and third columns contain the chromosome and map position.
#' @return Cumulative position (beyond chromosome) will be returned.
#'
#'
#'
cumsumPos <- function(map) {
  marker <- as.character(map[, 1])
  chr <- map[, 2]
  pos <- map[, 3]
  
  chr.tab <- table(chr)
  chr.max <- max(chr)
  chr.cum <- cumsum(chr.tab)
  cum.pos <- pos
  if (length(chr.tab) != 1) {
    for (i in 1:(chr.max - 1)) {
      cum.pos[(chr.cum[i] + 1):(chr.cum[i + 1])] <-
        pos[(chr.cum[i] + 1):(chr.cum[i + 1])] + cum.pos[chr.cum[i]]
    }
  }
  
  return(cum.pos)
}






#' Function to generate map for gene set
#'
#'
#' @param map Data frame with the marker names in the first column. The second and third columns contain the chromosome and map position.
#' @param gene.set Gene information with the format of a "data.frame" (whose dimension is (the number of gene) x 2).
#'            In the first column, you should assign the gene name. And in the second column, you should assign the names of each marker,
#'            which correspond to the marker names of "map" argument.
#' @param cumulative If this argument is TRUE, cumulative position will be returned.
#'
#' @return Map for gene set.
#'
#'
#'
genesetmap <- function(map, gene.set, cumulative = FALSE) {
  marker <- as.character(map[, 1])
  chr <- map[, 2]
  pos <- map[, 3]
  
  chr.tab <- table(chr)
  chr.max <- max(chr)
  chr.cum <- cumsum(chr.tab)
  cum.pos <- cumsumPos(map)
  
  gene.names <- as.character(gene.set[, 1])
  mark.id <- as.character(gene.set[, 2])
  gene.name <- as.character(unique(gene.names))
  n.scores <- length(unique(gene.set[, 1]))
  chr.set.mean <- pos.set.mean <- cum.pos.set.mean <- rep(NA, n.scores)
  ids <- as.data.frame(matrix(rep(NA, n.scores * 2), ncol = 2))
  
  marker.now <- gene.name
  for (k in 1:n.scores) {
    id <- mark.id[gene.names == gene.name[k]]
    ids[k, ] <- c(as.character(id[1]), as.character(id[length(id)]))
    num.sel <- match(id, map[, 1])
    chr.sel <- map[num.sel, 2]
    pos.sel <- map[num.sel, 3]
    cum.pos.sel <- cum.pos[num.sel]
    chr.set.mean[k] <- chr.sel[1]
    pos.set.mean[k] <- mean(pos.sel)
    cum.pos.set.mean[k] <- mean(cum.pos.sel)
  }
  
  if(!cumulative){
    map2 <- data.frame(marker = marker.now,
                       chr = chr.set.mean,
                       pos = pos.set.mean)
  }else{
    map2 <- data.frame(marker = marker.now,
                       chr = chr.set.mean,
                       pos = pos.set.mean,
                       cum.pos = cum.pos.set.mean)
  }
  
  return(map2)
}






#' Draw manhattan plot
#'
#' @param input Data frame of GWAS results where the first column is the marker names,
#' the second and third column is the chromosome amd map position, and the forth column is -log10(p) for each marker.
#' @param sig.level Significance level for the threshold. The default is 0.05.
#' @param method.thres Method for detemining threshold of significance. "BH" and "Bonferroni are offered.
#' @param y.max The maximum value for the vertical axis of manhattan plot. If NULL, automatically determined.
#' @param cex A numerical value giving the amount by which plotting text and symbols should be magnified relative to the default.
#' @param cex.lab The font size of the labels.
#' @param lwd.thres The line width for the threshold.
#' @param plot.col1 This argument determines the color of the manhattan plot.
#'  You should substitute this argument as color vector whose length is 2.
#'  plot.col1[1] for odd chromosomes and plot.col1[2] for even chromosomes.
#' @param cex.axis.x The font size of the x axis.
#' @param cex.axis.y The font size of the y axis.
#' @param plot.type  This argument determines the type of the manhattan plot. See the help page of "plot".
#' @param plot.pch This argument determines the shape of the dot of the manhattan plot. See the help page of "plot".
#'
#' @return Draw manhttan plot
#'
#'
#'
manhattan <- function(input, sig.level = 0.05, method.thres = "BH",
                      y.max = NULL, cex = 1, cex.lab = 1, lwd.thres = 1,
                      plot.col1 = c("dark blue", "cornflowerblue"),
                      cex.axis.x = 1, cex.axis.y = 1, plot.type = "p", plot.pch = 16){
  input <- input[!is.na(input[, 4]), , drop = FALSE]
  input <- input[order(input[, 2], input[, 3]), ]
  chroms <- unique(input[, 2])
  n.chrom <- length(chroms)
  chrom.start <- rep(0, n.chrom)
  chrom.mid <- rep(0, n.chrom)
  if (n.chrom > 1) {
    for (i in 1:(n.chrom - 1)) {
      chrom.start[i + 1] <- chrom.start[i] + max(input[which(input[, 2] == chroms[i]), 3]) + 1
    }
  }
  x.max <- chrom.start[n.chrom] + max(input[which(input[, 2] == chroms[n.chrom]), 3])
  if(is.null(y.max)){
    y.max <- max(input[, 4]) + 1
  }
  plot(0, 0, type = "n", xlim = c(0, x.max), ylim = c(0, y.max), ylab = expression(-log[10](italic(p))),
       xlab = "Chromosome", xaxt = "n", cex = cex, cex.lab = cex.lab, cex.axis = cex.axis.y)
  for (i in seq(1, n.chrom, by = 2)) {
    ix <- which(input[, 2] == chroms[i])
    chrom.mid[i] <- median(chrom.start[i] + input[ix, 3])
    points(chrom.start[i] + input[ix, 3], input[ix, 4], cex = cex,
           col = plot.col1[1], type = plot.type, pch = plot.pch)
  }
  if (n.chrom > 1) {
    for (i in seq(2, n.chrom, by = 2)) {
      ix <- which(input[, 2] == chroms[i])
      chrom.mid[i] <- median(chrom.start[i] + input[ix, 3])
      points(chrom.start[i] + input[ix, 3], input[ix, 4], cex = cex,
             col = plot.col1[2], type = plot.type, pch = plot.pch)
    }
  }
  
  
  
  threshold <- try(CalcThreshold(input, sig.level = sig.level, method = method.thres), silent = TRUE)
  if((!("try-error" %in% class(threshold))) & (!is.na(threshold))){
    lines(x = c(0, x.max), y = rep(threshold, 2), lty = 2, lwd = lwd.thres)
  }
  axis(side = 1, at = chrom.mid, labels = chroms, cex.axis = cex.axis.x)
}



#' Add points of -log10(p) corrected by kernel methods to manhattan plot
#'
#' @param input Data frame of GWAS results where the first column is the marker names,
#' the second and third column is the chromosome amd map position, and the forth column is -log10(p) for each marker.
#' @param checks The marker numbers whose -log10(p)s are corrected by kernel methods.
#' @param cex A numerical value giving the amount by which plotting text and symbols should be magnified relative to the default.
#' @param plot.col1 This argument determines the color of the manhattan plot.
#'  You should substitute this argument as a color vector whose length is 2.
#'  plot.col1[1] for odd chromosomes and plot.col1[2] for even chromosomes.
#' @param plot.col3 Color of -log10(p) corrected by kernel methods. plot.col3[1] for odd chromosomes and plot.col3[2] for even chromosomes
#' @param plot.type  This argument determines the type of the manhattan plot. See the help page of "plot".
#' @param plot.pch This argument determines the shape of the dot of the manhattan plot. See the help page of "plot".
#'
#' @return Draw manhttan plot
#'
#'
#'

manhattan.plus <- function(input, checks, cex = 1, plot.col1 = c("dark blue", "cornflowerblue"),
                           plot.col3 = c("red3", "orange3"), plot.type = "p", plot.pch = 16) {
  input <- input[!is.na(input[, 4]), , drop = FALSE]
  input <- input[order(input[, 2], input[, 3]), ]
  chroms <- unique(input[, 2])
  n.chrom <- length(chroms)
  chrom.start <- rep(0, n.chrom)
  chrom.mid <- rep(0, n.chrom)
  if (n.chrom > 1) {
    for (i in 1:(n.chrom - 1)) {
      chrom.start[i + 1] <- chrom.start[i] + max(input[which(input[, 2] == chroms[i]), 3]) + 1
    }
  }
  x.max <- chrom.start[n.chrom] + max(input[which(input[, 2] == chroms[n.chrom]), 3])
  for (i in seq(1, n.chrom, by = 2)) {
    ix <- checks[which(input[checks, 2] == chroms[i])]
    chrom.mid[i] <- median(chrom.start[i] + input[ix, 3])
    points(chrom.start[i] + input[ix, 3], input[ix, 4], cex = cex,
           col = plot.col3[1], type = plot.type, pch = plot.pch)
  }
  if (n.chrom > 1) {
    for (i in seq(2, n.chrom, by = 2)) {
      ix <- checks[which(input[checks, 2] == chroms[i])]
      chrom.mid[i] <- median(chrom.start[i] + input[ix, 3])
      points(chrom.start[i] + input[ix, 3], input[ix, 4], cex = cex, 
             col = plot.col3[2], type = plot.type, pch = plot.pch)
    }
  }
}



#' Draw manhattan plot (another method)
#'
#' @param input Data frame of GWAS results where the first column is the marker names,
#' the second and third column is the chromosome amd map position, and the forth column is -log10(p) for each marker.
#' @param sig.level Siginifincance level for the threshold. The default is 0.05.
#' @param method.thres Method for detemining threshold of significance. "BH" and "Bonferroni are offered.
#' @param cex A numerical value giving the amount by which plotting text and symbols should be magnified relative to the default.
#' @param plot.col2 Color of the manhattan plot. color changes with chromosome and it starts from plot.col2 + 1
#' (so plot.col2 = 1 means color starts from red.)
#' @param plot.type  This argument determines the type of the manhattan plot. See the help page of "plot".
#' @param plot.pch This argument determines the shape of the dot of the manhattan plot. See the help page of "plot".
#' @param cum.pos Cumulative position (over chromosomes) of each marker
#' @param lwd.thres The line width for the threshold.
#' @param cex.lab The font size of the labels.
#' @param cex.axis The font size of the axes.
#'
#' @return Draw manhttan plot
#'
#'
#'
manhattan2 <- function(input, sig.level = 0.05, method.thres = "BH", cex = 1, plot.col2 = 1,
                       plot.type = "p", plot.pch = 16, cum.pos = NULL, lwd.thres = 1,
                       cex.lab = 1, cex.axis = 1) {
  input <- input[!is.na(input[, 4]), , drop = FALSE]
  chr <- input[, 2]
  pos <- input[, 3]
  input <- input[order(chr, pos), ]
  chroms <- unique(chr)
  n.chrom <- length(chroms)
  chrom.start <- rep(0, n.chrom)
  chrom.mid <- rep(0, n.chrom)
  chr.tab <- table(chr)
  chr.max <- max(chr)
  chr.cum <- cumsum(chr.tab)
  if(is.null(cum.pos)){
    cum.pos <- pos
    if(length(chr.tab) != 1){
      for(i in 1:(chr.max - 1)){
        cum.pos[(chr.cum[i] + 1):(chr.cum[i + 1])] <- pos[(chr.cum[i] + 1):(chr.cum[i + 1])] + cum.pos[chr.cum[i]]
      }
    }
  }
  plot(cum.pos, input[, 4], col = chr + plot.col2, type = plot.type,
       pch = plot.pch, xlab = "Position (bp)", ylab = "-log10(p)",
       cex.lab = cex.lab, cex.axis = cex.axis, cex = cex)
  
  threshold <- try(CalcThreshold(input, sig.level = sig.level, method = method.thres), silent = TRUE)
  if((!("try-error" %in% class(threshold))) & (!is.na(threshold))){
    lines(x = c(0, max(cum.pos)), y = rep(threshold, 2), lty = 2, lwd = lwd.thres)
  }
}

#' Draw the effects of epistasis (3d plot and 2d plot)
#'
#'
#' @param input Data frame of GWAS results where the first column is the marker names,
#' the second and third column is the chromosome amd map position, and the forth column is -log10(p) for each marker.
#' @param cum.pos Cumulative position (over chromosomes) of each marker
#' @param plot.epi.3d If TRUE, draw 3d plot
#' @param plot.epi.2d If TRUE, draw 2d plot
#' @param main.epi.3d The title of 3d plot. If this argument is NULL, trait name is set as the title.
#' @param main.epi.2d The title of 2d plot. If this argument is NULL, trait name is set as the title.
#' @param saveName When drawing any plot, you can save plots in png format. In saveName, you should substitute the name you want to save.
#' When saveAt = NULL, the plot is not saved.
#'
#' @return Draw 3d plot and 2d plot to show epistatic effects
#'
#'
#'

manhattan3 <- function(input, cum.pos, plot.epi.3d = TRUE,
                       plot.epi.2d = TRUE, main.epi.3d = NULL,
                       main.epi.2d = NULL, saveName = NULL){
  x <- input[[2]]
  y <- input[[3]]
  z <- input[[4]]
  col.id <- rainbow(7)
  quan <- seq(0, max(z, na.rm = TRUE), length = 8)
  col.num <- rep(NA, length(z))
  for(j in 1:length(z)){
    if(z[j] != 0){
      col.num[j] <- max(which(z[j] - quan > 0))
    }else{
      col.num[j] <- 1
    }
  }
  
  if(plot.epi.3d){
    rgl::rgl.open()
    rgl::par3d(cex = 0.6)
    rgl::plot3d(x, y, z, col = col.id[col.num], type = "h", lwd = 3,
                xlim = c(0, max(cum.pos)), ylim = c(0, max(cum.pos)),
                xlab = "Position (bp)", ylab = "Position (bp)", zlab = "-log10(p)")
    rgl::legend3d("topright", legend = paste(round(rev(quan)[-1], 1), "~", round(rev(quan[-1]), 1)),
                  pch = 16, col = rev(rainbow(7)), cex = 0.6, inset = c(0.02))
    
    rgl::title3d(main = main.epi.3d)
  }
  
  pl.size <- 10 * z / max(z)
  if(!is.null(saveName)){
    if(plot.epi.3d){
      rgl.snapshot(paste(saveName, "_epistasis_3d_snapshot.png", sep = ""))
      if(length(grep("/", saveName)) != 0){
        spr <- strsplit(saveName, "/")[[1]]
        dir <- paste(spr[-length(spr)], collapse = "/")
        file_name <- spr[length(spr)]
      }else{
        dir <- ""
        file_name <- saveName
      }
      writeWebGL(dir = dir, filename = file.path(dir, paste(file_name, "_epistasis_3d_webplot.html")),
                 width = 1000, height = 1000)
      rgl.close()
    }
    
    if(plot.epi.2d){
      png(paste(saveName, "_epistasis_2d_plot.png", ""), width = 600, height = 500)
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
      par(mar = c(3, 3, 3, 6), xpd = T)
      plot(x, y, cex = pl.size, xlim = c(0, max(cum.pos)), ylim = c(0, max(cum.pos)), col = col.id[col.num], pch = 1)
      segments(0, 0, max(cum.pos), max(cum.pos), lty = "dotted")
      legend(oldpar$usr[2], oldpar$usr[4],
             legend = paste(round(rev(quan)[-1], 1), "~", round(rev(quan[-1]), 1)),
             pch = 1, col = rev(rainbow(7)), cex = 1)
      title(main = main.epi.2d)
      dev.off()
    }
  }else{
    if(plot.epi.2d){
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
      par(mar = c(3, 3, 3, 6), xpd = T)
      plot(x, y, cex = pl.size, xlim = c(0, max(cum.pos)), ylim = c(0, max(cum.pos)), col = col.id[col.num], pch = 1)
      segments(0, 0, max(cum.pos), max(cum.pos), lty = "dotted")
      legend(oldpar$usr[2], oldpar$usr[4],
             legend = paste(round(rev(quan)[-1], 1), "~", round(rev(quan[-1]), 1)),
             pch = 1, col = rev(rainbow(7)), cex = 1)
      title(main = main.epi.2d)
    }
  }
}



#' Draw qq plot
#'
#' @param scores A vector of -log10(p) for each marker
#'
#' @return Draw qq plot
#'
#'
#'
qq <- function(scores) {
  remove <- which(scores == 0)
  if (length(remove) > 0) {
    x <- sort(scores[-remove], decreasing = TRUE)
  }
  else {
    x <- sort(scores, decreasing = TRUE)
  }
  n <- length(x)
  unif.p <- -log10(ppoints(n))
  plot(unif.p, x, pch = 16, xlab = "Expected -log(p)",
       ylab = "Observed -log(p)")
  lines(c(0, max(unif.p)), c(0, max(unif.p)), lty = 2)
}



#' Calculate -log10(p) for single-SNP GWAS
#'
#' @description Calculate -log10(p) of each SNP by the Wald test.
#'
#'
#' @param M.now A n.sample x n.mark genotype matrix where n.sample is sample size and n.mark is the number of markers.
#' @param ZETA.now A list of variance (relationship) matrix (K; \eqn{m \times m}) and its design matrix (Z; \eqn{n \times m}) of random effects. You can use only one kernel matrix.
#' For example, ZETA = list(A = list(Z = Z, K = K))
#' Please set names of list "Z" and "K"!
#' @param y A \eqn{n \times 1} vector. A vector of phenotypic values should be used. NA is allowed.
#' @param X.now A \eqn{n \times p} matrix. You should assign mean vector (rep(1, n)) and covariates. NA is not allowed.
#' @param Hinv The inverse of \eqn{H = ZKZ' + \lambda I} where \eqn{\lambda = \sigma^2_e / \sigma^2_u}.
#' @param P3D When P3D = TRUE, variance components are estimated by REML only once, without any markers in the model.
#' When P3D = FALSE, variance components are estimated by REML for each marker separately.
#' @param optimizer The function used in the optimization process. We offer "optim", "optimx", and "nlminb" functions.
#' @param eigen.G A list with
#' \describe{
#' \item{$values}{Eigen values}
#' \item{$vectors}{Eigen vectors}
#' }
#' The result of the eigen decompsition of \eqn{G = ZKZ'}. You can use "spectralG.cpp" function in RAINBOWR.
#' If this argument is NULL, the eigen decomposition will be performed in this function.
#' We recommend you assign the result of the eigen decomposition beforehand for time saving.
#' @param min.MAF Specifies the minimum minor allele frequency (MAF).
#' If a marker has a MAF less than min.MAF, it is assigned a zero score.
#' @param count When count is TRUE, you can know how far RGWAS has ended with percent display.
#'
#' @return -log10(p) for each marker
#'
#' @references Kennedy, B.W., Quinton, M. and van Arendonk, J.A. (1992)
#' Estimation of effects of single genes on quantitative traits. J Anim Sci. 70(7): 2000-2012.
#'
#' Kang, H.M. et al. (2008) Efficient Control of Population Structure
#'  in Model Organism Association Mapping. Genetics. 178(3): 1709-1723.
#'
#' Kang, H.M. et al. (2010) Variance component model to account for sample
#'   structure in genome-wide association studies. Nat Genet. 42(4): 348-354.
#'
#' Zhang, Z. et al. (2010) Mixed linear model approach adapted for genome-wide
#'  association studies. Nat Genet. 42(4): 355-360.
#'
#'
#'
score.calc <- function(M.now, ZETA.now, y, X.now, Hinv, P3D = TRUE, optimizer = "nlminb",
                       eigen.G = NULL,  min.MAF = 0.02, count = TRUE) {
  n.mark <- ncol(M.now)
  scores <- array(NA, n.mark)
  
  lz <- length(ZETA.now)
  ZKZt <- matrix(0, nrow = length(y), ncol = length(y))
  for(list.no in lz){
    ZKZt.now <- tcrossprod(ZETA.now[[list.no]]$Z %*% ZETA.now[[list.no]]$K, ZETA.now[[list.no]]$Z)
    ZKZt <- ZKZt + ZKZt.now
  }
  rank.ZKZt <- Matrix::rankMatrix(ZKZt)[1]
  
  pb <- txtProgressBar(min = 1, max = n.mark, style = 3)
  n.mark2 <- n.mark - n.mark %% 100
  start.scorecalc <- Sys.time()
  for (i in 1:n.mark) {
    if(count){
      setTxtProgressBar(pb, i)
      if(n.mark > 100){
        if(i == (n.mark2 / 100 + 1) | i == (n.mark2 / 10 + 1) | i == (n.mark2 / 2 + 1)){
          cat("\n")
          end.scorecalc <- Sys.time()
          jikan.scorecalc <- (end.scorecalc - start.scorecalc) * (n.mark - i + 1) / (i - 1)
          print(paste0((i - 1) * 100 / n.mark2, "%...Done. ",
                       round(jikan.scorecalc, 2), " ", attr(jikan.scorecalc, "units"),
                       " to end.  Scheduled end time : ", end.scorecalc + jikan.scorecalc))
        }
      }
    }
    
    Mi <- M.now[, i]
    freq <- mean(Mi + 1, na.rm = TRUE) / 2
    MAF <- min(freq, 1 - freq)
    if (MAF >= min.MAF) {
      not.NA.geno <- which(!is.na(Mi))
      
      ni <- as.integer(min(length(not.NA.geno), rank.ZKZt))
      yi <- as.matrix(y[not.NA.geno])
      Xi <- cbind(X.now[not.NA.geno, ], Mi[not.NA.geno])
      p <- ncol(Xi)
      v1 <- 1
      v2 <- ni - p
      
      if (!P3D) {
        Xi <- make.full(Xi)
        if(length(ZETA.now) > 1){
          EMM.res <- EM3.cpp(y = yi, X0 = Xi, ZETA = ZETA.now, eigen.G = eigen.G, optimizer = optimizer,
                             tol = NULL, n.thres = 450, REML = TRUE, pred = FALSE)
        }else{
          EMM.res <- EMM.cpp(y = yi, X = Xi, ZETA = ZETA.now, eigen.G = eigen.G, optimizer = optimizer,
                             tol = NULL, n.thres = 450, REML = TRUE)
        }
        H2inv <- EMM.res$Hinv
      }else {
        H2inv <- Hinv[not.NA.geno, not.NA.geno]
      }
      
      beta.stat <- try(GWAS_F_test(y = yi, x = Xi, hinv = H2inv,
                                   v1 = v1, v2 = v2, p = p), silent = TRUE)
      
      if(!("try-error" %in% class(beta.stat))){
        scores[i] <- -log10(pbeta(beta.stat, v2 / 2, v1 / 2))
      }
    }
  }
  if (count) {
    cat("\n")
  }
  return(scores)
}




#' Calculate -log10(p) for single-SNP GWAS (multi-cores)
#'
#' @description Calculate -log10(p) of each SNP by the Wald test.
#'
#'
#' @param M.now A n.sample x n.mark genotype matrix where n.sample is sample size and n.mark is the number of markers.
#' @param ZETA.now A list of variance (relationship) matrix (K; \eqn{m \times m}) and its design matrix (Z; \eqn{n \times m}) of random effects. You can use only one kernel matrix.
#' For example, ZETA = list(A = list(Z = Z, K = K))
#' Please set names of list "Z" and "K"!
#' @param y A \eqn{n \times 1} vector. A vector of phenotypic values should be used. NA is allowed.
#' @param X.now A \eqn{n \times p} matrix. You should assign mean vector (rep(1, n)) and covariates. NA is not allowed.
#' @param Hinv The inverse of \eqn{H = ZKZ' + \lambda I} where \eqn{\lambda = \sigma^2_e / \sigma^2_u}.
#' @param n.core Setting n.core > 1 will enable parallel execution on a machine with multiple cores.
#' @param P3D When P3D = TRUE, variance components are estimated by REML only once, without any markers in the model.
#' When P3D = FALSE, variance components are estimated by REML for each marker separately.
#' @param optimizer The function used in the optimization process. We offer "optim", "optimx", and "nlminb" functions.
#' @param eigen.G A list with
#' \describe{
#' \item{$values}{Eigen values}
#' \item{$vectors}{Eigen vectors}
#' }
#' The result of the eigen decompsition of \eqn{G = ZKZ'}. You can use "spectralG.cpp" function in RAINBOWR.
#' If this argument is NULL, the eigen decomposition will be performed in this function.
#' We recommend you assign the result of the eigen decomposition beforehand for time saving.
#' @param min.MAF Specifies the minimum minor allele frequency (MAF).
#' If a marker has a MAF less than min.MAF, it is assigned a zero score.
#' @param count When count is TRUE, you can know how far RGWAS has ended with percent display.
#'
#' @return -log10(p) for each marker
#'
#' @references Kennedy, B.W., Quinton, M. and van Arendonk, J.A. (1992)
#' Estimation of effects of single genes on quantitative traits. J Anim Sci. 70(7): 2000-2012.
#'
#' Kang, H.M. et al. (2008) Efficient Control of Population Structure
#'  in Model Organism Association Mapping. Genetics. 178(3): 1709-1723.
#'
#' Kang, H.M. et al. (2010) Variance component model to account for sample
#'   structure in genome-wide association studies. Nat Genet. 42(4): 348-354.
#'
#' Zhang, Z. et al. (2010) Mixed linear model approach adapted for genome-wide
#'  association studies. Nat Genet. 42(4): 355-360.
#'
#'
#'
score.calc.MC <- function(M.now, ZETA.now, y, X.now, Hinv, n.core = 2, P3D = TRUE, optimizer = "nlminb",
                          eigen.G = NULL,  min.MAF = 0.02, count = TRUE) {
  n.mark <- ncol(M.now)
  
  lz <- length(ZETA.now)
  ZKZt <- matrix(0, nrow = length(y), ncol = length(y))
  for(list.no in lz){
    ZKZt.now <- tcrossprod(ZETA.now[[list.no]]$Z %*% ZETA.now[[list.no]]$K, ZETA.now[[list.no]]$Z)
    ZKZt <- ZKZt + ZKZt.now
  }
  rank.ZKZt <- Matrix::rankMatrix(ZKZt)[1]
  
  score.calc.MC.oneSNP <- function(markNo) {
    Mi <- M.now[, markNo]
    freq <- mean(Mi + 1, na.rm = TRUE) / 2
    MAF <- min(freq, 1 - freq)
    if (MAF >= min.MAF) {
      not.NA.geno <- which(!is.na(Mi))
      
      ni <- as.integer(min(length(not.NA.geno), rank.ZKZt))
      yi <- as.matrix(y[not.NA.geno])
      Xi <- cbind(X.now[not.NA.geno, ], Mi[not.NA.geno])
      p <- ncol(Xi)
      v1 <- 1
      v2 <- ni - p
      
      if (!P3D) {
        Xi <- make.full(Xi)
        if(length(ZETA.now) > 1){
          EMM.res <- EM3.cpp(y = yi, X0 = Xi, ZETA = ZETA.now, eigen.G = eigen.G, optimizer = optimizer,
                             tol = NULL, n.thres = 450, REML = TRUE, pred = FALSE)
        }else{
          EMM.res <- EMM.cpp(y = yi, X = Xi, ZETA = ZETA.now, eigen.G = eigen.G, optimizer = optimizer,
                             tol = NULL, n.thres = 450, REML = TRUE)
        }
        H2inv <- EMM.res$Hinv
      }else {
        H2inv <- Hinv[not.NA.geno, not.NA.geno]
      }
      
      beta.stat <- try(GWAS_F_test(y = yi, x = Xi, hinv = H2inv,
                                   v1 = v1, v2 = v2, p = p), silent = TRUE)
      
      if(!("try-error" %in% class(beta.stat))){
        scores.now <- -log10(pbeta(beta.stat, v2 / 2, v1 / 2))
      } else {
        scores.now <- NA
      }
    } else {
      scores.now <- NA
    }
    return(scores.now)
  }
  if (count) {
    cat("\n")
  }
  
  scores <- unlist(pbmcapply::pbmclapply(X = 1:n.mark, FUN = score.calc.MC.oneSNP, mc.cores = n.core))
  
  return(scores)
}






#' Change a matrix to full-rank matrix
#'
#' @param X A \eqn{n \times p} matrix which you want to change into full-rank matrix.
#'
#' @return A full-rank matrix
#'
#'
#'
make.full <- function(X) {
  svd.X <- svd(X)
  r <- max(which(svd.X$d > 1e-08))
  if(r < ncol(X)){
    newX <- svd.X$u[, 1:r, drop = FALSE]
    rownames(newX) <- rownames(X)
    colnames(newX) <- paste0("svd", 1:r)
  }else{
    newX <- X
  }
  return(newX)
}




#' Calculate -log10(p) of each SNP-set by the LR test
#'
#' @description This function calculates -log10(p) of each SNP-set by the LR (likelihood-ratio) test.
#' First, the function solves the multi-kernel mixed model and calaculates the maximum restricted log likelihood.
#' Then it performs the LR test by using the fact that the deviance
#'
#' \deqn{D = 2 \times (LL _ {alt} - LL _ {null})}
#'
#' follows the chi-square distribution.
#'
#'
#' @param M.now A n.sample x n.mark genotype matrix where n.sample is sample size and n.mark is the number of markers.
#' @param ZETA.now A list of variance (relationship) matrix (K; \eqn{m \times m}) and its design matrix (Z; \eqn{n \times m}) of random effects. You can use only one kernel matrix.
#' For example, ZETA = list(A = list(Z = Z, K = K))
#' Please set names of list "Z" and "K"!
#' @param y A \eqn{n \times 1} vector. A vector of phenotypic values should be used. NA is allowed.
#' @param X.now A \eqn{n \times p} matrix. You should assign mean vector (rep(1, n)) and covariates. NA is not allowed.
#' @param LL0 The log-likelihood for the null model.
#' @param eigen.SGS A list with
#' \describe{
#' \item{$values}{Eigen values}
#' \item{$vectors}{Eeigen vectors}
#' }
#' The result of the eigen decompsition of \eqn{SGS}, where \eqn{S = I - X(X'X)^{-1}X'}, \eqn{G = ZKZ'}.
#' You can use "spectralG.cpp" function in RAINBOWR.
#' If this argument is NULL, the eigen decomposition will be performed in this function.
#' We recommend you assign the result of the eigen decomposition beforehand for time saving.
#' @param eigen.G A list with
#' \describe{
#' \item{$values}{Eigen values}
#' \item{$vectors}{Eigen vectors}
#' }
#' The result of the eigen decompsition of \eqn{G = ZKZ'}. You can use "spectralG.cpp" function in RAINBOWR.
#' If this argument is NULL, the eigen decomposition will be performed in this function.
#' We recommend you assign the result of the eigen decomposition beforehand for time saving.
#' @param map Data frame of map information where the first column is the marker names,
#' the second and third column is the chromosome amd map position, and the forth column is -log10(p) for each marker.
#' @param kernel.method It determines how to calculate kernel. There are three methods.
#' \describe{
#' \item{"gaussian"}{It is the default method. Gaussian kernel is calculated by distance matrix.}
#' \item{"exponential"}{When this method is selected, exponential kernel is calculated by distance matrix.}
#' \item{"linear"}{When this method is selected, linear kernel is calculated by NOIA methods for additive GRM.}
#'}
#' @param kernel.h The hyper-parameter for gaussian or exponential kernel.
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
#' @param optimizer The function used in the optimization process. We offer "optim", "optimx", and "nlminb" functions.
#' @param chi0.mixture RAINBOWR assumes the deviance is considered to follow a x chisq(df = 0) + (1 - a) x chisq(df = r).
#' where r is the degree of freedom.
#' The argument chi0.mixture is a (0 <= a < 1), and default is 0.5.
#' @param weighting.center In kernel-based GWAS, weights according to the Gaussian distribution (centered on the tested SNP) are taken into account when calculating the kernel if Rainbow = TRUE.
#'           If weighting.center = FALSE, weights are not taken into account.
#' @param weighting.other You can set other weights in addition to weighting.center. The length of this argument should be equal to the number of SNPs.
#'           For example, you can assign SNP effects from the information of gene annotation.
#' @param gene.set If you have information of gene, you can use it to perform kernel-based GWAS.
#'            You should assign your gene information to gene.set in the form of a "data.frame" (whose dimension is (the number of gene) x 2).
#'            In the first column, you should assign the gene name. And in the second column, you should assign the names of each marker,
#'            which correspond to the marker names of "geno" argument.
#' @param min.MAF Specifies the minimum minor allele frequency (MAF).
#' If a marker has a MAF less than min.MAF, it is assigned a zero score.
#' @param count When count is TRUE, you can know how far RGWAS has ended with percent display.
#'
#' @return -log10(p) for each SNP-set
#'
#' @references Listgarten, J. et al. (2013) A powerful and efficient set test
#'  for genetic markers that handles confounders. Bioinformatics. 29(12): 1526-1533.
#'
#' Lippert, C. et al. (2014) Greater power and computational efficiency for kernel-based
#'  association testing of sets of genetic variants. Bioinformatics. 30(22): 3206-3214.
#'
#'
#'
score.calc.LR <- function(M.now, y, X.now, ZETA.now, LL0, eigen.SGS = NULL, eigen.G = NULL, optimizer = "nlminb",
                          map, kernel.method = "linear", kernel.h = "tuned", haplotype = TRUE, num.hap = NULL,
                          test.effect = "additive", window.size.half = 5, window.slide = 1,
                          chi0.mixture = 0.5, weighting.center = TRUE, weighting.other = NULL,
                          gene.set = NULL, min.MAF = 0.02, count = TRUE){
  n <- length(y)
  
  chr <- map[, 2]
  chr.tab <- table(chr)
  chr.max <- max(chr)
  chr.cum <- cumsum(chr.tab)
  n.scores.each <- (chr.tab + (window.slide - 1)) %/% window.slide
  cum.n.scores <- cumsum(n.scores.each)
  if(is.null(gene.set)){
    n.scores <- sum(n.scores.each)
  }else{
    gene.names <- as.character(gene.set[, 1])
    mark.id <- as.character(gene.set[, 2])
    gene.name <- as.character(unique(gene.names))
    n.scores <- length(unique(gene.set[, 1]))
  }
  
  if(kernel.method == "linear"){
    scores <- matrix(NA, nrow = n.scores, ncol = length(test.effect))
  }else{
    scores <- matrix(NA, nrow = n.scores, ncol =  1)
  }
  window.centers <- rep(NA, n.scores)
  
  probaa <- apply(M.now == -1, 2, mean)
  probAa <- apply(M.now == 0, 2, mean)
  freq <- probaa + probAa / 2
  MAF <- pmin(freq, 1 - freq)
  MAF.D <- pmin(probAa, 1 - probAa)
  
  
  pb <- txtProgressBar(min = 1, max = n.scores, style = 3)
  n.scores2 <- n.scores - n.scores %% 100
  start.scorecalc <- Sys.time()
  for (i in 1:n.scores) {
    if(count){
      setTxtProgressBar(pb, i)
      if(n.scores > 100){
        if(i == (n.scores2 / 100 + 1) | i == (n.scores2 / 10 + 1) | i == (n.scores2 / 2 + 1)){
          cat("\n")
          end.scorecalc <- Sys.time()
          jikan.scorecalc <- (end.scorecalc - start.scorecalc) * (n.scores - i + 1) / (i - 1)
          print(paste0((i - 1) * 100 / n.scores2, "%...Done. ",
                       round(jikan.scorecalc, 2), " ", attr(jikan.scorecalc, "units"),
                       " to end.  Scheduled end time : ", end.scorecalc + jikan.scorecalc))
        }
      }
    }
    
    
    
    if(is.null(gene.set)){
      i.chr <- min(which(i - cum.n.scores <= 0))
      if(i.chr >= 2){
        window.center <- window.slide * (i - cum.n.scores[i.chr - 1] - 1) + chr.cum[i.chr - 1] + 1
      }else{
        window.center <- window.slide * (i - 1) + 1
      }
      names(window.center) <- i.chr
      window.centers[i] <- window.center
      Theories1 <-  window.center < window.size.half + 1
      for(r in 1:(chr.max - 1)){
        Theory1 <- chr.cum[r] < window.center & window.center < window.size.half + 1 + chr.cum[r]
        Theories1 <- c(Theories1, Theory1)
      }
      rule1 <- sum(Theories1) != 0
      
      
      Theories2  <- NULL
      for(r in 1:chr.max){
        Theory2 <- chr.cum[r] - (window.size.half + 1) < window.center & window.center <= chr.cum[r]
        Theories2 <- c(Theories2, Theory2)
      }
      rule2 <- sum(Theories2) != 0
      
      
      
      if(rule1 & rule2){
        Mis.range.0 <- which(chr == i.chr)
        Mis.range.02 <- which(chr == i.chr) - window.center + 1 + window.size.half
      }else{
        if(rule1){
          near.min <- c(0, chr.cum)[which.min(abs(window.center - c(0, chr.cum)))]
          Mis.range.0 <- (near.min + 1):(window.center + window.size.half)
          Mis.range.02 <- (2 * window.size.half + 2 - length(Mis.range.0)):(2 * window.size.half + 1)
        }else{
          if(rule2){
            near.max <- c(0, chr.cum)[which.min(abs(window.center - c(0, chr.cum)))]
            Mis.range.0 <- (window.center - window.size.half):near.max
            Mis.range.02 <- 1:length(Mis.range.0)
          }else{
            Mis.range.0 <- (window.center - window.size.half):(window.center + window.size.half)
            Mis.range.02 <- 1:(2 * window.size.half + 1)
          }
        }
      }
    }else{
      mark.name.now <- mark.id[gene.names == gene.name[i]]
      Mis.range.0 <- match(mark.name.now, map[, 1])
      Mis.range.02 <- 1:length(Mis.range.0)
      weighting.center <- FALSE
    }
    
    Mis.0 <- M.now[, Mis.range.0, drop = FALSE]
    MAF.cut <- MAF[Mis.range.0] >= min.MAF
    if(any(test.effect %in% c("dominance", "additive+dominance"))){
      Mis.D.0 <- M.now[, Mis.range.0, drop = FALSE]
      MAF.cut.D <- MAF.D[Mis.range.0] > 0
    }else{
      MAF.cut.D <- rep(TRUE, length(MAF.cut))
    }
    
    if(any(MAF.cut)){
      Mis.0 <- Mis.0[, MAF.cut, drop = FALSE]
      Mis.range <- Mis.range.0[MAF.cut]
      Mis.range2 <- Mis.range.02[MAF.cut]
      window.size <- ncol(Mis.0)
      if(any(MAF.cut.D)){
        if(any(test.effect %in% c("dominance", "additive+dominance"))){
          Mis.D.0 <- Mis.D.0[, MAF.cut.D, drop = FALSE]
          Mis.range.D <- Mis.range.0[MAF.cut.D]
          Mis.range2.D <- Mis.range.02[MAF.cut.D]
          window.size.D <- ncol(Mis.D.0)
        }
      }
      
      if(haplotype){
        if(is.null(num.hap)){
          Mis.fac <- factor(apply(Mis.0, 1, function(x) paste(x, collapse = "")))
          Mis <- Mis.0[!duplicated(as.numeric(Mis.fac)), , drop = FALSE]
          bango <- as.factor(as.numeric(Mis.fac))
          levels(bango) <- order(unique(bango))
          bango <- as.numeric(as.character(bango))
          if(any(MAF.cut.D)){
            if(any(test.effect %in% c("dominance", "additive+dominance"))){
              Mis.D.fac <- factor(apply(Mis.D.0, 1, function(x) paste(x, collapse = "")))
              Mis.D <- Mis.D.0[!duplicated(as.numeric(Mis.D.fac)), , drop = FALSE]
              
              bango.D <- as.factor(as.numeric(Mis.D.fac))
              levels(bango.D) <- order(unique(bango.D))
              bango.D <- as.numeric(as.character(bango.D))
            }
          }
        }else{
          kmed.res <- cluster::pam(Mis.0, k = num.hap)
          Mis <- kmed.res$medoids
          bango <- kmed.res$clustering
          if(any(MAF.cut.D)){
            if(any(test.effect %in% c("dominance", "additive+dominance"))){
              kmed.res.D <- cluster::pam(Mis.D.0, k = num.hap)
              Mis.D <- kmed.res.D$medoids
              bango.D <- kmed.res.D$clustering
            }
          }
        }
        Z.part <- as.matrix(Matrix::sparseMatrix(i = 1:nrow(M.now), j = bango, x = rep(1, nrow(M.now)),
                                                 dims = c(nrow(M.now), nrow(Mis))))
        if(any(MAF.cut.D)){
          if(any(test.effect %in% c("dominance", "additive+dominance"))){
            Z.part.D <- as.matrix(Matrix::sparseMatrix(i = 1:nrow(M.now), j = bango.D, x = rep(1, nrow(M.now)),
                                                       dims = c(nrow(M.now), nrow(Mis.D))))
          }
        }
      }else{
        Mis <- Mis.0
        Mis.D <- Mis.D.0
        Z.part <- Z.part.D <- diag(nrow(M.now))
      }
      
      if(window.size != 1){
        if(weighting.center){
          weight.Mis <- dnorm((-window.size.half):(window.size.half), 0, window.size.half / 2)[Mis.range2]
          weight.Mis <- weight.Mis / apply(Mis, 2, sd)
          if(!is.null(weighting.other)){
            weight.Mis <- weight.Mis * weighting.other[Mis.range]
          }
          weight.Mis <- weight.Mis * window.size / sum(weight.Mis)
        }else{
          weight.Mis <- rep(1, window.size)
          weight.Mis <- weight.Mis / apply(Mis, 2, sd)
          if(!is.null(weighting.other)){
            weight.Mis <- weight.Mis * weighting.other[Mis.range]
          }
          weight.Mis <- weight.Mis * window.size / sum(weight.Mis)
        }
      }else{
        weight.Mis <- 1
      }
      
      if(any(MAF.cut.D)){
        if(window.size != 1){
          if(any(test.effect %in% c("dominance", "additive+dominance"))){
            if(weighting.center){
              weight.Mis.D <- dnorm((-window.size.half):(window.size.half), 0, window.size.half / 2)[Mis.range2.D]
              weight.Mis.D <- weight.Mis.D / apply(Mis.D, 2, sd)
              if(!is.null(weighting.other)){
                weight.Mis.D <- weight.Mis.D * weighting.other[Mis.range.D]
              }
              weight.Mis.D <- weight.Mis.D * window.size.D / sum(weight.Mis.D)
            }else{
              weight.Mis.D <- rep(1, window.size.D)
              weight.Mis.D <- weight.Mis.D / apply(Mis.D, 2, sd)
              if(!is.null(weighting.other)){
                weight.Mis.D <- weight.Mis.D * weighting.other[Mis.range.D]
              }
              weight.Mis.D <- weight.Mis.D * window.size.D / sum(weight.Mis.D)
            }
          }
        }else{
          weight.Mis.D <- 1
        }
      }
      
      if(kernel.method != "linear"){
        if(ncol(Mis) != 1){
          Mis.weighted <- t(apply(Mis, 1, function(x) x * weight.Mis))
        }else{
          Mis.weighted <- as.matrix(apply(Mis, 1, function(x) x * weight.Mis))
        }
        
        K.SNP <- calcGRM(genoMat = Mis.weighted,
                         methodGRM = kernel.method,
                         kernel.h = kernel.h,
                         returnWMat = FALSE)
        
        if(length(ZETA.now) == 1){
          Gammas0 <- list(K = K.SNP)
          Ws0 <- list(W = Z.part)
          Zs0 <- list(Z = diag(nrow(Mis.0)))
          EMM.res2 <- try(EM3.linker.cpp(y0 = y, X0 = X.now, ZETA = ZETA.now, optimizer = optimizer,
                                         Zs0 = Zs0, Ws0 = Ws0, Gammas0 = Gammas0,
                                         gammas.diag = FALSE, X.fix = TRUE, tol = NULL,
                                         eigen.SGS = eigen.SGS, eigen.G = eigen.G,
                                         REML = TRUE, pred = FALSE), silent = TRUE)
          if("try-error" %in% class(EMM.res2)){
            ZETA.now2 <- c(ZETA.now, list(part = list(Z = Z.part, K = K.SNP)))
            EMM.res2 <- try(EM3.cpp(y = y, X0 = X.now, ZETA = ZETA.now2, tol = NULL, optimizer = optimizer,
                                    REML = TRUE, pred = FALSE), silent = TRUE)
          }
        }else{
          ZETA.now2 <- c(ZETA.now, list(part = list(Z = Z.part, K = K.SNP)))
          EMM.res2 <- try(EM3.cpp(y = y, X0 = X.now, ZETA = ZETA.now2, tol = NULL, optimizer = optimizer,
                                  REML = TRUE, pred = FALSE), silent = TRUE)
        }
        
        if(!("try-error" %in% class(EMM.res2))){
          LL2s <- EMM.res2$LL
        }else{
          LL2s <- LL0
        }
        df <- 1
      }else{
        test.no <- match(test.effect, c("additive", "dominance", "additive+dominance"))
        if(length(test.no) == 0){
          stop("The effect to test should be 'additive', 'dominance' or 'additive+dominance'!")
        }
        
        if(any(test.effect %in% c("additive", "additive+dominance"))){
          W.A <- calcGRM(genoMat = Mis,
                         methodGRM = "addNOIA",
                         returnWMat = TRUE,
                         probaa = probaa[Mis.range],
                         probAa = probAa[Mis.range])
        }
        
        if(any(MAF.cut.D)){
          if(any(test.effect %in% c("dominance", "additive+dominance"))){
            W.D <- calcGRM(genoMat = Mis.D,
                           methodGRM = "domNOIA",
                           returnWMat = TRUE,
                           probaa = probaa[Mis.range.D],
                           probAa = probAa[Mis.range.D])
          }
        }
        
        if(length(ZETA.now) == 1){
          if(1 %in% test.no){
            Ws0.A <- list(W.A = W.A)
            Zs0.A <- list(W.A = Z.part)
            Gammas0.A <- list(W.A = diag(weight.Mis ^ 2))
          }
          
          if(any(MAF.cut.D)){
            if(2 %in% test.no){
              Ws0.D <- list(W.D = W.D)
              Zs0.D <- list(W.D = Z.part.D)
              Gammas0.D <- list(W.D = diag(weight.Mis.D ^ 2))
            }
            
            if(3 %in% test.no){
              Ws0.AD <- list(W.A = W.A, W.D = W.D)
              Zs0.AD <- list(W.A = Z.part, W.D = Z.part.D)
              Gammas0.AD <- list(W.A = diag(weight.Mis ^ 2), W.D = diag(weight.Mis.D ^ 2))
            }
          }
        }else{
          if(1 %in% test.no){
            K.A.part <- W.A %*% (t(W.A) * weight.Mis)
            ZETA.now2.A <- c(ZETA.now, list(part.A = list(Z = Z.part, K = K.A.part)))
          }
          
          if(any(MAF.cut.D)){
            if(2 %in% test.no){
              K.D.part <- W.D %*% (t(W.D) * weight.Mis.D)
              ZETA.now2.D <- c(ZETA.now, list(part.D = list(Z = Z.part.D, K = K.D.part)))
            }
            
            if(3 %in% test.no){
              K.A.part <- W.A %*% (t(W.A) * weight.Mis)
              K.D.part <- W.D %*% (t(W.D) * weight.Mis.D)
              ZETA.now2.AD <- c(ZETA.now, list(part.A = list(Z = Z.part, K = K.A.part)),
                                list(part.D = list(Z = Z.part.D, K = K.D.part)))
            }
          }
        }
        
        LL2s <- df <- rep(NA, length(test.no))
        for(j in 1:length(test.no)){
          test.no.now <- test.no[j]
          if(length(ZETA.now) == 1){
            if(test.no.now == 1){
              EMM.res2 <- try(EM3.linker.cpp(y0 = y, X0 = X.now, ZETA = ZETA.now, optimizer = optimizer,
                                             Zs0 = Zs0.A, Ws0 = Ws0.A, Gammas0 = Gammas0.A,
                                             gammas.diag = TRUE, X.fix = TRUE, tol = NULL,
                                             eigen.SGS = eigen.SGS, eigen.G = eigen.G,
                                             REML = TRUE, pred = FALSE), silent = TRUE)
            }
            
            if(test.no.now == 2){
              EMM.res2 <- try(EM3.linker.cpp(y0 = y, X0 = X.now, ZETA = ZETA.now, optimizer = optimizer,
                                             Zs0 = Zs0.D, Ws0 = Ws0.D, Gammas0 = Gammas0.D,
                                             gammas.diag = TRUE, X.fix = TRUE, tol = NULL,
                                             eigen.SGS = eigen.SGS, eigen.G = eigen.G,
                                             REML = TRUE, pred = FALSE), silent = TRUE)
            }
            
            if(test.no.now == 3){
              EMM.res2 <- try(EM3.linker.cpp(y0 = y, X0 = X.now, ZETA = ZETA.now, optimizer = optimizer,
                                             Zs0 = Zs0.AD, Ws0 = Ws0.AD, Gammas0 = Gammas0.AD,
                                             gammas.diag = TRUE, X.fix = TRUE, tol = NULL,
                                             eigen.SGS = eigen.SGS, eigen.G = eigen.G,
                                             REML = TRUE, pred = FALSE), silent = TRUE)
            }
            
            if("try-error" %in% class(EMM.res2)){
              if(1 %in% test.no){
                K.A.part <- W.A %*% (t(W.A) * weight.Mis)
                ZETA.now2.A <- c(ZETA.now, list(part.A = list(Z = Z.part, K = K.A.part)))
              }
              
              if(any(MAF.cut.D)){
                if(2 %in% test.no){
                  K.D.part <- W.D %*% (t(W.D) * weight.Mis.D)
                  ZETA.now2.D <- c(ZETA.now, list(part.D = list(Z = Z.part.D, K = K.D.part)))
                }
                
                if(3 %in% test.no){
                  K.A.part <- W.A %*% (t(W.A) * weight.Mis)
                  K.D.part <- W.D %*% (t(W.D) * weight.Mis.D)
                  ZETA.now2.AD <- c(ZETA.now, list(part.A = list(Z = Z.part, K = K.A.part)),
                                    list(part.D = list(Z = Z.part.D, K = K.D.part)))
                }
              }
              
              if(test.no.now == 1){
                EMM.res2 <- try(EM3.cpp(y = y, X0 = X.now, ZETA = ZETA.now2.A, tol = NULL, optimizer = optimizer,
                                        REML = TRUE, pred = FALSE), silent = TRUE)
              }
              
              if(test.no.now == 2){
                EMM.res2 <- try(EM3.cpp(y = y, X0 = X.now, ZETA = ZETA.now2.D, tol = NULL, optimizer = optimizer,
                                        REML = TRUE, pred = FALSE), silent = TRUE)
              }
              
              if(test.no.now == 3){
                EMM.res2 <- try(EM3.cpp(y = y, X0 = X.now, ZETA = ZETA.now2.AD, tol = NULL, optimizer = optimizer,
                                        REML = TRUE, pred = FALSE), silent = TRUE)
              }
            }
          }else{
            if(test.no.now == 1){
              EMM.res2 <- try(EM3.cpp(y = y, X0 = X.now, ZETA = ZETA.now2.A, tol = NULL, optimizer = optimizer,
                                      REML = TRUE, pred = FALSE), silent = TRUE)
            }
            
            if(test.no.now == 2){
              EMM.res2 <- try(EM3.cpp(y = y, X0 = X.now, ZETA = ZETA.now2.D, tol = NULL, optimizer = optimizer,
                                      REML = TRUE, pred = FALSE), silent = TRUE)
            }
            
            if(test.no.now == 3){
              EMM.res2 <- try(EM3.cpp(y = y, X0 = X.now, ZETA = ZETA.now2.AD, tol = NULL, optimizer = optimizer,
                                      REML = TRUE, pred = FALSE), silent = TRUE)
            }
          }
          
          if(!("try-error" %in% class(EMM.res2))){
            LL2 <- EMM.res2$LL
          }else{
            LL2 <- LL0
          }
          LL2s[j] <- LL2
        }
        df[test.no == 1] <- 1
        df[test.no == 2] <- 1
        df[test.no == 3] <- 2
      }
      
      
      deviances <- 2 * (LL2s - LL0)
      scores.now <- ifelse(deviances <= 0, 0, -log10((1 - chi0.mixture) *
                                                       pchisq(q = deviances, df = df, lower.tail = FALSE)))
      scores[i, ] <- scores.now
    }
  }
  
  if(is.null(gene.set)){
    rownames(scores) <- window.centers
  }else{
    rownames(scores) <- gene.name
  }
  
  if(kernel.method == "linear"){
    colnames(scores) <- test.effect
  }else{
    colnames(scores) <- kernel.method
  }
  
  if (count) {
    cat("\n")
  }
  return(scores)
}






#' Calculate -log10(p) of each SNP-set by the LR test (multi-cores)
#'
#' @description This function calculates -log10(p) of each SNP-set by the LR (likelihood-ratio) test.
#' First, the function solves the multi-kernel mixed model and calaculates the maximum restricted log likelihood.
#' Then it performs the LR test by using the fact that the deviance
#'
#' \deqn{D = 2 \times (LL _ {alt} - LL _ {null})}
#'
#' follows the chi-square distribution.
#'
#'
#' @param M.now A n.sample x n.mark genotype matrix where n.sample is sample size and n.mark is the number of markers.
#' @param ZETA.now A list of variance (relationship) matrix (K; \eqn{m \times m}) and its design matrix (Z; \eqn{n \times m}) of random effects. You can use only one kernel matrix.
#' For example, ZETA = list(A = list(Z = Z, K = K))
#' Please set names of list "Z" and "K"!
#' @param y A \eqn{n \times 1} vector. A vector of phenotypic values should be used. NA is allowed.
#' @param X.now A \eqn{n \times p} matrix. You should assign mean vector (rep(1, n)) and covariates. NA is not allowed.
#' @param LL0 The log-likelihood for the null model.
#' @param eigen.SGS A list with
#' \describe{
#' \item{$values}{Eigen values}
#' \item{$vectors}{Eigen vectors}
#' }
#' The result of the eigen decompsition of \eqn{SGS}, where \eqn{S = I - X(X'X)^{-1}X'}, \eqn{G = ZKZ'}.
#' You can use "spectralG.cpp" function in RAINBOWR.
#' If this argument is NULL, the eigen decomposition will be performed in this function.
#' We recommend you assign the result of the eigen decomposition beforehand for time saving.
#' @param eigen.G A list with
#' \describe{
#' \item{$values}{Eigen values}
#' \item{$vectors}{Eigen vectors}
#' }
#' The result of the eigen decompsition of \eqn{G = ZKZ'}. You can use "spectralG.cpp" function in RAINBOWR.
#' If this argument is NULL, the eigen decomposition will be performed in this function.
#' We recommend you assign the result of the eigen decomposition beforehand for time saving.
#' @param n.core Setting n.core > 1 will enable parallel execution on a machine with multiple cores.
#' @param map Data frame of map information where the first column is the marker names,
#' the second and third column is the chromosome amd map position, and the forth column is -log10(p) for each marker.
#' @param kernel.method It determines how to calculate kernel. There are three methods.
#' \describe{
#' \item{"gaussian"}{It is the default method. Gaussian kernel is calculated by distance matrix.}
#' \item{"exponential"}{When this method is selected, exponential kernel is calculated by distance matrix.}
#' \item{"linear"}{When this method is selected, linear kernel is calculated by NOIA methods for additive GRM.}
#'}
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
#' @param optimizer The function used in the optimization process. We offer "optim", "optimx", and "nlminb" functions.
#' @param chi0.mixture RAINBOWR assumes the deviance is considered to follow a x chisq(df = 0) + (1 - a) x chisq(df = r).
#' where r is the degree of freedom.
#' The argument chi0.mixture is a (0 <= a < 1), and default is 0.5.
#' @param weighting.center In kernel-based GWAS, weights according to the Gaussian distribution (centered on the tested SNP) are taken into account when calculating the kernel if Rainbow = TRUE.
#'           If weighting.center = FALSE, weights are not taken into account.
#' @param weighting.other You can set other weights in addition to weighting.center. The length of this argument should be equal to the number of SNPs.
#'           For example, you can assign SNP effects from the information of gene annotation.
#' @param gene.set If you have information of gene, you can use it to perform kernel-based GWAS.
#'            You should assign your gene information to gene.set in the form of a "data.frame" (whose dimension is (the number of gene) x 2).
#'            In the first column, you should assign the gene name. And in the second column, you should assign the names of each marker,
#'            which correspond to the marker names of "geno" argument.
#' @param min.MAF Specifies the minimum minor allele frequency (MAF).
#' If a marker has a MAF less than min.MAF, it is assigned a zero score.
#' @param count When count is TRUE, you can know how far RGWAS has ended with percent display.
#'
#' @return -log10(p) for each SNP-set
#'
#' @references Listgarten, J. et al. (2013) A powerful and efficient set test
#'  for genetic markers that handles confounders. Bioinformatics. 29(12): 1526-1533.
#'
#' Lippert, C. et al. (2014) Greater power and computational efficiency for kernel-based
#'  association testing of sets of genetic variants. Bioinformatics. 30(22): 3206-3214.
#'
#'
#'
score.calc.LR.MC <- function(M.now, y, X.now, ZETA.now, LL0, eigen.SGS = NULL, eigen.G = NULL, n.core = 2,
                             map, kernel.method = "linear", kernel.h = "tuned", haplotype = TRUE, num.hap = NULL,
                             test.effect = "additive", window.size.half = 5, window.slide = 1, optimizer = "nlminb",
                             chi0.mixture = 0.5, weighting.center = TRUE, weighting.other = NULL,
                             gene.set = NULL, min.MAF = 0.02, count = TRUE){
  n <- length(y)
  
  chr <- map[, 2]
  chr.tab <- table(chr)
  chr.max <- max(chr)
  chr.cum <- cumsum(chr.tab)
  n.scores.each <- (chr.tab + (window.slide - 1)) %/% window.slide
  cum.n.scores <- cumsum(n.scores.each)
  if(is.null(gene.set)){
    n.scores <- sum(n.scores.each)
  }else{
    gene.names <- as.character(gene.set[, 1])
    mark.id <- as.character(gene.set[, 2])
    gene.name <- as.character(unique(gene.names))
    n.scores <- length(unique(gene.set[, 1]))
  }
  
  if(kernel.method == "linear"){
    ncol.scores <- length(test.effect)
  }else{
    ncol.scores <- 1
  }
  window.centers <- rep(NA, n.scores)
  probaa <- apply(M.now == -1, 2, mean)
  probAa <- apply(M.now == 0, 2, mean)
  freq <- probaa + probAa / 2
  MAF <- pmin(freq, 1 - freq)
  MAF.D <- pmin(probAa, 1 - probAa)
  
  
  score.calc.LR.MC.oneSNP <- function(markNo) {
    if(is.null(gene.set)){
      markNo.chr <- min(which(markNo - cum.n.scores <= 0))
      if(markNo.chr >= 2){
        window.center <- window.slide * (markNo - cum.n.scores[markNo.chr - 1] - 1) + chr.cum[markNo.chr - 1] + 1
      }else{
        window.center <- window.slide * (markNo - 1) + 1
      }
      names(window.center) <- markNo.chr
      Theories1 <-  window.center < window.size.half + 1
      for(r in 1:(chr.max - 1)){
        Theory1 <- chr.cum[r] < window.center & window.center < window.size.half + 1 + chr.cum[r]
        Theories1 <- c(Theories1, Theory1)
      }
      rule1 <- sum(Theories1) != 0
      
      
      Theories2  <- NULL
      for(r in 1:chr.max){
        Theory2 <- chr.cum[r] - (window.size.half + 1) < window.center & window.center <= chr.cum[r]
        Theories2 <- c(Theories2, Theory2)
      }
      rule2 <- sum(Theories2) != 0
      
      
      
      if(rule1 & rule2){
        Mis.range.0 <- which(chr == markNo.chr)
        Mis.range.02 <- which(chr == markNo.chr) - window.center + 1 + window.size.half
      }else{
        if(rule1){
          near.min <- c(0, chr.cum)[which.min(abs(window.center - c(0, chr.cum)))]
          Mis.range.0 <- (near.min + 1):(window.center + window.size.half)
          Mis.range.02 <- (2 * window.size.half + 2 - length(Mis.range.0)):(2 * window.size.half + 1)
        }else{
          if(rule2){
            near.max <- c(0, chr.cum)[which.min(abs(window.center - c(0, chr.cum)))]
            Mis.range.0 <- (window.center - window.size.half):near.max
            Mis.range.02 <- 1:length(Mis.range.0)
          }else{
            Mis.range.0 <- (window.center - window.size.half):(window.center + window.size.half)
            Mis.range.02 <- 1:(2 * window.size.half + 1)
          }
        }
      }
    }else{
      mark.name.now <- mark.id[gene.names == gene.name[markNo]]
      Mis.range.0 <- match(mark.name.now, map[, 1])
      Mis.range.02 <- 1:length(Mis.range.0)
      weighting.center <- FALSE
    }
    
    Mis.0 <- M.now[, Mis.range.0, drop = FALSE]
    MAF.cut <- MAF[Mis.range.0] >= min.MAF
    if(any(test.effect %in% c("dominance", "additive+dominance"))){
      Mis.D.0 <- M.now[, Mis.range.0, drop = FALSE]
      MAF.cut.D <- MAF.D[Mis.range.0] > 0
    }else{
      MAF.cut.D <- rep(TRUE, length(MAF.cut))
    }
    
    if(any(MAF.cut)){
      Mis.0 <- Mis.0[, MAF.cut, drop = FALSE]
      Mis.range <- Mis.range.0[MAF.cut]
      Mis.range2 <- Mis.range.02[MAF.cut]
      window.size <- ncol(Mis.0)
      if(any(MAF.cut.D)){
        if(any(test.effect %in% c("dominance", "additive+dominance"))){
          Mis.D.0 <- Mis.D.0[, MAF.cut.D, drop = FALSE]
          Mis.range.D <- Mis.range.0[MAF.cut.D]
          Mis.range2.D <- Mis.range.02[MAF.cut.D]
          window.size.D <- ncol(Mis.D.0)
        }
      }
      
      if(haplotype){
        if(is.null(num.hap)){
          Mis.fac <- factor(apply(Mis.0, 1, function(x) paste(x, collapse = "")))
          Mis <- Mis.0[!duplicated(as.numeric(Mis.fac)), , drop = FALSE]
          bango <- as.factor(as.numeric(Mis.fac))
          levels(bango) <- order(unique(bango))
          bango <- as.numeric(as.character(bango))
          if(any(MAF.cut.D)){
            if(any(test.effect %in% c("dominance", "additive+dominance"))){
              Mis.D.fac <- factor(apply(Mis.D.0, 1, function(x) paste(x, collapse = "")))
              Mis.D <- Mis.D.0[!duplicated(as.numeric(Mis.D.fac)), , drop = FALSE]
              
              bango.D <- as.factor(as.numeric(Mis.D.fac))
              levels(bango.D) <- order(unique(bango.D))
              bango.D <- as.numeric(as.character(bango.D))
            }
          }
        }else{
          kmed.res <- cluster::pam(Mis.0, k = num.hap)
          Mis <- kmed.res$medoids
          bango <- kmed.res$clustering
          if(any(MAF.cut.D)){
            if(any(test.effect %in% c("dominance", "additive+dominance"))){
              kmed.res.D <- cluster::pam(Mis.D.0, k = num.hap)
              Mis.D <- kmed.res.D$medoids
              bango.D <- kmed.res.D$clustering
            }
          }
        }
        Z.part <- as.matrix(Matrix::sparseMatrix(i = 1:nrow(M.now), j = bango, x = rep(1, nrow(M.now)),
                                                 dims = c(nrow(M.now), nrow(Mis))))
        if(any(MAF.cut.D)){
          if(any(test.effect %in% c("dominance", "additive+dominance"))){
            Z.part.D <- as.matrix(Matrix::sparseMatrix(i = 1:nrow(M.now), j = bango.D, x = rep(1, nrow(M.now)),
                                                       dims = c(nrow(M.now), nrow(Mis.D))))
          }
        }
      }else{
        Mis <- Mis.0
        Mis.D <- Mis.D.0
        Z.part <- Z.part.D <- diag(nrow(M.now))
      }
      
      if(window.size != 1){
        if(weighting.center){
          weight.Mis <- dnorm((-window.size.half):(window.size.half), 0, window.size.half / 2)[Mis.range2]
          weight.Mis <- weight.Mis / apply(Mis, 2, sd)
          if(!is.null(weighting.other)){
            weight.Mis <- weight.Mis * weighting.other[Mis.range]
          }
          weight.Mis <- weight.Mis * window.size / sum(weight.Mis)
        }else{
          weight.Mis <- rep(1, window.size)
          weight.Mis <- weight.Mis / apply(Mis, 2, sd)
          if(!is.null(weighting.other)){
            weight.Mis <- weight.Mis * weighting.other[Mis.range]
          }
          weight.Mis <- weight.Mis * window.size / sum(weight.Mis)
        }
      }else{
        weight.Mis <- 1
      }
      
      if(any(MAF.cut.D)){
        if(window.size != 1){
          if(any(test.effect %in% c("dominance", "additive+dominance"))){
            if(weighting.center){
              weight.Mis.D <- dnorm((-window.size.half):(window.size.half), 0, window.size.half / 2)[Mis.range2.D]
              weight.Mis.D <- weight.Mis.D / apply(Mis.D, 2, sd)
              if(!is.null(weighting.other)){
                weight.Mis.D <- weight.Mis.D * weighting.other[Mis.range.D]
              }
              weight.Mis.D <- weight.Mis.D * window.size.D / sum(weight.Mis.D)
            }else{
              weight.Mis.D <- rep(1, window.size.D)
              weight.Mis.D <- weight.Mis.D / apply(Mis.D, 2, sd)
              if(!is.null(weighting.other)){
                weight.Mis.D <- weight.Mis.D * weighting.other[Mis.range.D]
              }
              weight.Mis.D <- weight.Mis.D * window.size.D / sum(weight.Mis.D)
            }
          }
        }else{
          weight.Mis.D <- 1
        }
      }
      
      if(kernel.method != "linear"){
        if(ncol(Mis) != 1){
          Mis.weighted <- t(apply(Mis, 1, function(x) x * weight.Mis))
        }else{
          Mis.weighted <- as.matrix(apply(Mis, 1, function(x) x * weight.Mis))
        }
        
        K.SNP <- calcGRM(genoMat = Mis.weighted,
                         methodGRM = kernel.method,
                         kernel.h = kernel.h,
                         returnWMat = FALSE)
        
        
        if(length(ZETA.now) == 1){
          Gammas0 <- list(K = K.SNP)
          Ws0 <- list(W = Z.part)
          Zs0 <- list(Z = diag(nrow(Mis.0)))
          EMM.res2 <- try(EM3.linker.cpp(y0 = y, X0 = X.now, ZETA = ZETA.now, optimizer = optimizer,
                                         Zs0 = Zs0, Ws0 = Ws0, Gammas0 = Gammas0,
                                         gammas.diag = FALSE, X.fix = TRUE, tol = NULL,
                                         eigen.SGS = eigen.SGS, eigen.G = eigen.G,
                                         REML = TRUE, pred = FALSE), silent = TRUE)
          if("try-error" %in% class(EMM.res2)){
            ZETA.now2 <- c(ZETA.now, list(part = list(Z = Z.part, K = K.SNP)))
            EMM.res2 <- try(EM3.cpp(y = y, X0 = X.now, ZETA = ZETA.now2, tol = NULL, optimizer = optimizer,
                                    REML = TRUE, pred = FALSE), silent = TRUE)
          }
        }else{
          ZETA.now2 <- c(ZETA.now, list(part = list(Z = Z.part, K = K.SNP)))
          EMM.res2 <- try(EM3.cpp(y = y, X0 = X.now, ZETA = ZETA.now2, tol = NULL, optimizer = optimizer,
                                  REML = TRUE, pred = FALSE), silent = TRUE)
        }
        
        if(!("try-error" %in% class(EMM.res2))){
          LL2s <- EMM.res2$LL
        }else{
          LL2s <- LL0
        }
        df <- 1
      }else{
        test.no <- match(test.effect, c("additive", "dominance", "additive+dominance"))
        if(length(test.no) == 0){
          stop("The effect to test should be 'additive', 'dominance' or 'additive+dominance'!")
        }
        
        if(any(test.effect %in% c("additive", "additive+dominance"))){
          W.A <- calcGRM(genoMat = Mis,
                         methodGRM = "addNOIA",
                         returnWMat = TRUE,
                         probaa = probaa[Mis.range],
                         probAa = probAa[Mis.range])
        }
        
        if(any(MAF.cut.D)){
          if(any(test.effect %in% c("dominance", "additive+dominance"))){
            W.D <- calcGRM(genoMat = Mis.D,
                           methodGRM = "domNOIA",
                           returnWMat = TRUE,
                           probaa = probaa[Mis.range.D],
                           probAa = probAa[Mis.range.D])
          }
        }
        
        if(length(ZETA.now) == 1){
          if(1 %in% test.no){
            Ws0.A <- list(W.A = W.A)
            Zs0.A <- list(W.A = Z.part)
            Gammas0.A <- list(W.A = diag(weight.Mis ^ 2))
          }
          
          if(any(MAF.cut.D)){
            if(2 %in% test.no){
              Ws0.D <- list(W.D = W.D)
              Zs0.D <- list(W.D = Z.part.D)
              Gammas0.D <- list(W.D = diag(weight.Mis.D ^ 2))
            }
            
            if(3 %in% test.no){
              Ws0.AD <- list(W.A = W.A, W.D = W.D)
              Zs0.AD <- list(W.A = Z.part, W.D = Z.part.D)
              Gammas0.AD <- list(W.A = diag(weight.Mis ^ 2), W.D = diag(weight.Mis.D ^ 2))
            }
          }
        }else{
          if(1 %in% test.no){
            K.A.part <- W.A %*% (t(W.A) * weight.Mis)
            ZETA.now2.A <- c(ZETA.now, list(part.A = list(Z = Z.part, K = K.A.part)))
          }
          
          if(any(MAF.cut.D)){
            if(2 %in% test.no){
              K.D.part <- W.D %*% (t(W.D) * weight.Mis.D)
              ZETA.now2.D <- c(ZETA.now, list(part.D = list(Z = Z.part.D, K = K.D.part)))
            }
            
            if(3 %in% test.no){
              K.A.part <- W.A %*% (t(W.A) * weight.Mis)
              K.D.part <- W.D %*% (t(W.D) * weight.Mis.D)
              ZETA.now2.AD <- c(ZETA.now, list(part.A = list(Z = Z.part, K = K.A.part)),
                                list(part.D = list(Z = Z.part.D, K = K.D.part)))
            }
          }
        }
        
        LL2s <- df <- rep(NA, length(test.no))
        for(j in 1:length(test.no)){
          test.no.now <- test.no[j]
          if(length(ZETA.now) == 1){
            if(test.no.now == 1){
              EMM.res2 <- try(EM3.linker.cpp(y0 = y, X0 = X.now, ZETA = ZETA.now, optimizer = optimizer,
                                             Zs0 = Zs0.A, Ws0 = Ws0.A, Gammas0 = Gammas0.A,
                                             gammas.diag = TRUE, X.fix = TRUE, tol = NULL,
                                             eigen.SGS = eigen.SGS, eigen.G = eigen.G,
                                             REML = TRUE, pred = FALSE), silent = TRUE)
            }
            
            if(test.no.now == 2){
              EMM.res2 <- try(EM3.linker.cpp(y0 = y, X0 = X.now, ZETA = ZETA.now, optimizer = optimizer,
                                             Zs0 = Zs0.D, Ws0 = Ws0.D, Gammas0 = Gammas0.D,
                                             gammas.diag = TRUE, X.fix = TRUE, tol = NULL,
                                             eigen.SGS = eigen.SGS, eigen.G = eigen.G,
                                             REML = TRUE, pred = FALSE), silent = TRUE)
            }
            
            if(test.no.now == 3){
              EMM.res2 <- try(EM3.linker.cpp(y0 = y, X0 = X.now, ZETA = ZETA.now, optimizer = optimizer,
                                             Zs0 = Zs0.AD, Ws0 = Ws0.AD, Gammas0 = Gammas0.AD,
                                             gammas.diag = TRUE, X.fix = TRUE, tol = NULL,
                                             eigen.SGS = eigen.SGS, eigen.G = eigen.G,
                                             REML = TRUE, pred = FALSE), silent = TRUE)
            }
            
            if("try-error" %in% class(EMM.res2)){
              if(1 %in% test.no){
                K.A.part <- W.A %*% (t(W.A) * weight.Mis)
                ZETA.now2.A <- c(ZETA.now, list(part.A = list(Z = Z.part, K = K.A.part)))
              }
              
              if(any(MAF.cut.D)){
                if(2 %in% test.no){
                  K.D.part <- W.D %*% (t(W.D) * weight.Mis.D)
                  ZETA.now2.D <- c(ZETA.now, list(part.D = list(Z = Z.part.D, K = K.D.part)))
                }
                
                if(3 %in% test.no){
                  K.A.part <- W.A %*% (t(W.A) * weight.Mis)
                  K.D.part <- W.D %*% (t(W.D) * weight.Mis.D)
                  ZETA.now2.AD <- c(ZETA.now, list(part.A = list(Z = Z.part, K = K.A.part)),
                                    list(part.D = list(Z = Z.part.D, K = K.D.part)))
                }
              }
              
              if(test.no.now == 1){
                EMM.res2 <- try(EM3.cpp(y = y, X0 = X.now, ZETA = ZETA.now2.A, tol = NULL, optimizer = optimizer,
                                        REML = TRUE, pred = FALSE), silent = TRUE)
              }
              
              if(test.no.now == 2){
                EMM.res2 <- try(EM3.cpp(y = y, X0 = X.now, ZETA = ZETA.now2.D, tol = NULL, optimizer = optimizer,
                                        REML = TRUE, pred = FALSE), silent = TRUE)
              }
              
              if(test.no.now == 3){
                EMM.res2 <- try(EM3.cpp(y = y, X0 = X.now, ZETA = ZETA.now2.AD, tol = NULL, optimizer = optimizer,
                                        REML = TRUE, pred = FALSE), silent = TRUE)
              }
            }
          }else{
            if(test.no.now == 1){
              EMM.res2 <- try(EM3.cpp(y = y, X0 = X.now, ZETA = ZETA.now2.A, tol = NULL, optimizer = optimizer,
                                      REML = TRUE, pred = FALSE), silent = TRUE)
            }
            
            if(test.no.now == 2){
              EMM.res2 <- try(EM3.cpp(y = y, X0 = X.now, ZETA = ZETA.now2.D, tol = NULL, optimizer = optimizer,
                                      REML = TRUE, pred = FALSE), silent = TRUE)
            }
            
            if(test.no.now == 3){
              EMM.res2 <- try(EM3.cpp(y = y, X0 = X.now, ZETA = ZETA.now2.AD, tol = NULL, optimizer = optimizer,
                                      REML = TRUE, pred = FALSE), silent = TRUE)
            }
          }
          
          if(!("try-error" %in% class(EMM.res2))){
            LL2 <- EMM.res2$LL
          }else{
            LL2 <- LL0
          }
          LL2s[j] <- LL2
        }
        df[test.no == 1] <- 1
        df[test.no == 2] <- 1
        df[test.no == 3] <- 2
      }
      
      
      deviances <- 2 * (LL2s - LL0)
      scores.now <- ifelse(deviances <= 0, 0, -log10((1 - chi0.mixture) *
                                                       pchisq(q = deviances, df = df, lower.tail = FALSE)))
    } else {
      scores.now <- rep(NA, ncol.scores)
    }
    
    if(is.null(gene.set)){
      return(list(scores = scores.now, window.center = window.center))
    } else {
      return(list(scores = scores.now))
    }
  }
  
  all.res <- pbmcapply::pbmclapply(1:n.scores, score.calc.LR.MC.oneSNP, mc.cores = n.core)
  scores <- unlist(lapply(all.res, function(x) x$scores))
  scores <- matrix(scores, nrow = n.scores, ncol = ncol.scores, byrow = TRUE)
  
  if(is.null(gene.set)){
    window.centers <- unlist(lapply(all.res, function(x) x$window.center))
    rownames(scores) <- window.centers
  }else{
    rownames(scores) <- gene.name
  }
  
  if(kernel.method == "linear"){
    colnames(scores) <- test.effect
  }else{
    colnames(scores) <- kernel.method
  }
  
  if (count) {
    cat("\n")
  }
  return(scores)
}










#' Calculate -log10(p) of each SNP-set by the score test
#'
#' @description This function calculates -log10(p) of each SNP-set by the score test.
#' First, the function calculates the score statistic
#' without solving the multi-kernel mixed model for each SNP-set.
#' Then it performs the score test by using the fact that the score statistic follows the chi-square distribution.
#'
#'
#' @param M.now A n.sample x n.mark genotype matrix where n.sample is sample size and n.mark is the number of markers.
#' @param ZETA.now A list of variance (relationship) matrix (K; \eqn{m \times m}) and its design matrix (Z; \eqn{n \times m}) of random effects. You can use only one kernel matrix.
#' For example, ZETA = list(A = list(Z = Z, K = K))
#' Please set names of list "Z" and "K"!
#' @param y A \eqn{n \times 1} vector. A vector of phenotypic values should be used. NA is allowed.
#' @param X.now A \eqn{n \times p} matrix. You should assign mean vector (rep(1, n)) and covariates. NA is not allowed.
#' @param LL0 The log-likelihood for the null model.
#' @param Gu A \eqn{n \times n} matrix. You should assign \eqn{ZKZ'}, where K is covariance (relationship) matrix and Z is its design matrix.
#' @param Ge A \eqn{n \times n} matrix. You should assign identity matrix I (diag(n)).
#' @param P0 \eqn{n \times n} matrix. The Moore-Penrose generalized inverse of \eqn{SV0S}, where \eqn{S = X(X'X)^{-1}X'} and
#' \eqn{V0 = \sigma^2_u Gu + \sigma^2_e Ge}. \eqn{\sigma^2_u} and \eqn{\sigma^2_e} are estimators of the null model.
#' @param map Data frame of map information where the first column is the marker names,
#' the second and third column is the chromosome amd map position, and the forth column is -log10(p) for each marker.
#' @param kernel.method It determines how to calculate kernel. There are three methods.
#' \describe{
#' \item{"gaussian"}{It is the default method. Gaussian kernel is calculated by distance matrix.}
#' \item{"exponential"}{When this method is selected, exponential kernel is calculated by distance matrix.}
#' \item{"linear"}{When this method is selected, linear kernel is calculated by NOIA methods for additive GRM.}
#'}
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
#' @param chi0.mixture RAINBOWR assumes the test statistic \eqn{l1' F l1} is considered to follow a x chisq(df = 0) + (1 - a) x chisq(df = r).
#' where l1 is the first derivative of the log-likelihood and F is the Fisher information. And r is the degree of freedom.
#' The argument chi0.mixture is a (0 <= a < 1), and default is 0.5.
#' @param weighting.center In kernel-based GWAS, weights according to the Gaussian distribution (centered on the tested SNP) are taken into account when calculating the kernel if Rainbow = TRUE.
#'           If weighting.center = FALSE, weights are not taken into account.
#' @param weighting.other You can set other weights in addition to weighting.center. The length of this argument should be equal to the number of SNPs.
#'           For example, you can assign SNP effects from the information of gene annotation.
#' @param gene.set If you have information of gene, you can use it to perform kernel-based GWAS.
#'            You should assign your gene information to gene.set in the form of a "data.frame" (whose dimension is (the number of gene) x 2).
#'            In the first column, you should assign the gene name. And in the second column, you should assign the names of each marker,
#'            which correspond to the marker names of "geno" argument.
#' @param min.MAF Specifies the minimum minor allele frequency (MAF).
#' If a marker has a MAF less than min.MAF, it is assigned a zero score.
#' @param count When count is TRUE, you can know how far RGWAS has ended with percent display.
#'
#' @return -log10(p) for each SNP-set
#'
#' @references Listgarten, J. et al. (2013) A powerful and efficient set test
#'  for genetic markers that handles confounders. Bioinformatics. 29(12): 1526-1533.
#'
#' Lippert, C. et al. (2014) Greater power and computational efficiency for kernel-based
#'  association testing of sets of genetic variants. Bioinformatics. 30(22): 3206-3214.
#'
#'
#'
#'
score.calc.score <- function(M.now, y, X.now, ZETA.now, LL0, Gu, Ge, P0,
                             map, kernel.method = "linear", kernel.h = "tuned", haplotype = TRUE, num.hap = NULL,
                             test.effect = "additive", window.size.half = 5, window.slide = 1,
                             chi0.mixture = 0.5, weighting.center = TRUE, weighting.other = NULL,
                             gene.set = NULL, min.MAF = 0.02, count = TRUE){
  chr <- map[, 2]
  chr.tab <- table(chr)
  chr.max <- max(chr)
  chr.cum <- cumsum(chr.tab)
  n.scores.each <- (chr.tab + (window.slide - 1)) %/% window.slide
  cum.n.scores <- cumsum(n.scores.each)
  if(is.null(gene.set)){
    n.scores <- sum(n.scores.each)
  }else{
    gene.names <- as.character(gene.set[, 1])
    mark.id <- as.character(gene.set[, 2])
    gene.name <- as.character(unique(gene.names))
    n.scores <- length(unique(gene.set[, 1]))
  }
  
  if(kernel.method == "linear"){
    scores <- matrix(NA, nrow = n.scores, ncol = length(test.effect))
  }else{
    scores <- matrix(NA, nrow = n.scores, ncol =  1)
  }
  
  window.centers <- rep(NA, n.scores)
  probaa <- apply(M.now == -1, 2, mean)
  probAa <- apply(M.now == 0, 2, mean)
  freq <- probaa + probAa / 2
  MAF <- pmin(freq, 1 - freq)
  MAF.D <- pmin(probAa, 1 - probAa)
  
  
  pb <- txtProgressBar(min = 1, max = n.scores, style = 3)
  n.scores2 <- n.scores - n.scores %% 100
  start.scorecalc <- Sys.time()
  for (i in 1:n.scores) {
    if(count){
      setTxtProgressBar(pb, i)
      if(n.scores > 100){
        if(i == (n.scores2 / 100 + 1) | i == (n.scores2 / 10 + 1) | i == (n.scores2 / 2 + 1)){
          cat("\n")
          end.scorecalc <- Sys.time()
          jikan.scorecalc <- (end.scorecalc - start.scorecalc) * (n.scores - i + 1) / (i - 1)
          print(paste0((i - 1) * 100 / n.scores2, "%...Done. ",
                       round(jikan.scorecalc, 2), " ", attr(jikan.scorecalc, "units"),
                       " to end.  Scheduled end time : ", end.scorecalc + jikan.scorecalc))
        }
      }
    }
    
    if(is.null(gene.set)){
      i.chr <- min(which(i - cum.n.scores <= 0))
      if(i.chr >= 2){
        window.center <- window.slide * (i - cum.n.scores[i.chr - 1] - 1) + chr.cum[i.chr - 1] + 1
      }else{
        window.center <- window.slide * (i - 1) + 1
      }
      names(window.center) <- i.chr
      window.centers[i] <- window.center
      Theories1 <-  window.center < window.size.half + 1
      for(r in 1:(chr.max - 1)){
        Theory1 <- chr.cum[r] < window.center & window.center < window.size.half + 1 + chr.cum[r]
        Theories1 <- c(Theories1, Theory1)
      }
      rule1 <- sum(Theories1) != 0
      
      
      Theories2  <- NULL
      for(r in 1:chr.max){
        Theory2 <- chr.cum[r] - (window.size.half + 1) < window.center & window.center <= chr.cum[r]
        Theories2 <- c(Theories2, Theory2)
      }
      rule2 <- sum(Theories2) != 0
      
      
      if(rule1 & rule2){
        Mis.range.0 <- which(chr == i.chr)
        Mis.range.02 <- which(chr == i.chr) - window.center + 1 + window.size.half
      }else{
        if(rule1){
          near.min <- c(0, chr.cum)[which.min(abs(window.center - c(0, chr.cum)))]
          Mis.range.0 <- (near.min + 1):(window.center + window.size.half)
          Mis.range.02 <- (2 * window.size.half + 2 - length(Mis.range.0)):(2 * window.size.half + 1)
        }else{
          if(rule2){
            near.max <- c(0, chr.cum)[which.min(abs(window.center - c(0, chr.cum)))]
            Mis.range.0 <- (window.center - window.size.half):near.max
            Mis.range.02 <- 1:length(Mis.range.0)
          }else{
            Mis.range.0 <- (window.center - window.size.half):(window.center + window.size.half)
            Mis.range.02 <- 1:(2 * window.size.half + 1)
          }
        }
      }
    }else{
      mark.name.now <- mark.id[gene.names == gene.name[i]]
      Mis.range.0 <- match(mark.name.now, map[, 1])
      Mis.range.02 <- 1:length(Mis.range.0)
      weighting.center <- FALSE
    }
    
    Mis.0 <- M.now[, Mis.range.0, drop = FALSE]
    MAF.cut <- MAF[Mis.range.0] >= min.MAF
    if(any(test.effect %in% c("dominance", "additive+dominance"))){
      Mis.D.0 <- M.now[, Mis.range.0, drop = FALSE]
      MAF.cut.D <- MAF.D[Mis.range.0] > 0
    }else{
      MAF.cut.D <- rep(TRUE, length(MAF.cut))
    }
    
    if(any(MAF.cut)){
      Mis.0 <- Mis.0[, MAF.cut, drop = FALSE]
      Mis.range <- Mis.range.0[MAF.cut]
      Mis.range2 <- Mis.range.02[MAF.cut]
      window.size <- ncol(Mis.0)
      if(any(MAF.cut.D)){
        if(any(test.effect %in% c("dominance", "additive+dominance"))){
          Mis.D.0 <- Mis.D.0[, MAF.cut.D, drop = FALSE]
          Mis.range.D <- Mis.range.0[MAF.cut.D]
          Mis.range2.D <- Mis.range.02[MAF.cut.D]
          window.size.D <- ncol(Mis.D.0)
        }
      }
      
      if(haplotype){
        if(is.null(num.hap)){
          Mis.fac <- factor(apply(Mis.0, 1, function(x) paste(x, collapse = "")))
          Mis <- Mis.0[!duplicated(as.numeric(Mis.fac)), , drop = FALSE]
          
          bango <- as.factor(as.numeric(Mis.fac))
          levels(bango) <- order(unique(bango))
          bango <- as.numeric(as.character(bango))
          if(any(MAF.cut.D)){
            if(any(test.effect %in% c("dominance", "additive+dominance"))){
              Mis.D.fac <- factor(apply(Mis.D.0, 1, function(x) paste(x, collapse = "")))
              Mis.D <- Mis.D.0[!duplicated(as.numeric(Mis.D.fac)), , drop = FALSE]
              
              
              bango.D <- as.factor(as.numeric(Mis.D.fac))
              levels(bango.D) <- order(unique(bango.D))
              bango.D <- as.numeric(as.character(bango.D))
            }
          }
        }else{
          kmed.res <- cluster::pam(Mis.0, k = num.hap)
          Mis <- kmed.res$medoids
          bango <- kmed.res$clustering
          if(any(MAF.cut.D)){
            if(any(test.effect %in% c("dominance", "additive+dominance"))){
              kmed.res.D <- cluster::pam(Mis.D.0, k = num.hap)
              Mis.D <- kmed.res.D$medoids
              bango.D <- kmed.res.D$clustering
            }
          }
        }
        Z.part <- as.matrix(Matrix::sparseMatrix(i = 1:nrow(M.now), j = bango, x = rep(1, nrow(M.now)),
                                                 dims = c(nrow(M.now), nrow(Mis))))
        if(any(MAF.cut.D)){
          if(any(test.effect %in% c("dominance", "additive+dominance"))){
            Z.part.D <- as.matrix(Matrix::sparseMatrix(i = 1:nrow(M.now), j = bango.D, x = rep(1, nrow(M.now)),
                                                       dims = c(nrow(M.now), nrow(Mis.D))))
          }
        }
      }else{
        Mis <- Mis.0
        Mis.D <- Mis.D.0
        Z.part <- Z.part.D <- diag(nrow(M.now))
      }
      
      if(window.size != 1){
        if(weighting.center){
          weight.Mis <- dnorm((-window.size.half):(window.size.half), 0, window.size.half / 2)[Mis.range2]
          weight.Mis <- weight.Mis / apply(Mis, 2, sd)
          if(!is.null(weighting.other)){
            weight.Mis <- weight.Mis * weighting.other[Mis.range]
          }
          weight.Mis <- weight.Mis * window.size / sum(weight.Mis)
        }else{
          weight.Mis <- rep(1, window.size)
          weight.Mis <- weight.Mis / apply(Mis, 2, sd)
          if(!is.null(weighting.other)){
            weight.Mis <- weight.Mis * weighting.other[Mis.range]
          }
          weight.Mis <- weight.Mis * window.size / sum(weight.Mis)
        }
      }else{
        weight.Mis <- 1
      }
      
      if(any(MAF.cut.D)){
        if(window.size != 1){
          if(any(test.effect %in% c("dominance", "additive+dominance"))){
            if(weighting.center){
              weight.Mis.D <- dnorm((-window.size.half):(window.size.half), 0, window.size.half / 2)[Mis.range2.D]
              weight.Mis.D <- weight.Mis.D / apply(Mis.D, 2, sd)
              if(!is.null(weighting.other)){
                weight.Mis.D <- weight.Mis.D * weighting.other[Mis.range.D]
              }
              weight.Mis.D <- weight.Mis.D * window.size.D / sum(weight.Mis.D)
            }else{
              weight.Mis.D <- rep(1, window.size.D)
              weight.Mis.D <- weight.Mis.D / apply(Mis.D, 2, sd)
              if(!is.null(weighting.other)){
                weight.Mis.D <- weight.Mis.D * weighting.other[Mis.range.D]
              }
              weight.Mis.D <- weight.Mis.D * window.size.D / sum(weight.Mis.D)
            }
          }
        }else{
          weight.Mis.D <- 1
        }
      }
      
      if(kernel.method != "linear"){
        if(ncol(Mis) != 1){
          Mis.weighted <- t(apply(Mis, 1, function(x) x * weight.Mis))
        }else{
          Mis.weighted <- as.matrix(apply(Mis, 1, function(x) x * weight.Mis))
        }
        
        K.SNP <- calcGRM(genoMat = Mis.weighted,
                         methodGRM = kernel.method,
                         kernel.h = kernel.h,
                         returnWMat = FALSE)        
        
        Ws <- list(W = Z.part)
        Gammas <- list(Gamma = K.SNP)
        scores.now <- score.linker.cpp(y, Ws = Ws, Gammas = Gammas,
                                       gammas.diag = FALSE, Gu = Gu, Ge = Ge,
                                       P0 = P0, chi0.mixture = chi0.mixture)
      }else{
        test.no <- match(test.effect, c("additive", "dominance", "additive+dominance"))
        if(length(test.no) == 0){
          stop("The effect to test should be 'additive', 'dominance' or 'additive+dominance'!")
        }
        
        if(any(test.effect %in% c("additive", "additive+dominance"))){
          W.A <- calcGRM(genoMat = Mis,
                         methodGRM = "addNOIA",
                         returnWMat = TRUE,
                         probaa = probaa[Mis.range],
                         probAa = probAa[Mis.range])
        }
        if(any(MAF.cut.D)){
          if(any(test.effect %in% c("dominance", "additive+dominance"))){
            W.D <- calcGRM(genoMat = Mis.D,
                           methodGRM = "domNOIA",
                           returnWMat = TRUE,
                           probaa = probaa[Mis.range.D],
                           probAa = probAa[Mis.range.D])
          }
        }
        
        if(1 %in% test.no){
          Ws.A <- list(W.A = Z.part %*% W.A)
          Gammas.A <- list(W.A = diag(weight.Mis ^ 2))
        }
        
        if(any(MAF.cut.D)){
          if(2 %in% test.no){
            Ws.D <- list(W.D = Z.part.D %*% W.D)
            Gammas.D <- list(W.D = diag(weight.Mis.D ^ 2))
          }
          
          if(3 %in% test.no){
            Ws.AD <- list(W.A = Z.part %*% W.A, Z.part.D %*% W.D)
            Gammas.AD <- list(W.A = diag(weight.Mis ^ 2), W.D = diag(weight.Mis.D ^ 2))
          }
        }
        
        scores.now <- rep(NA, length(test.no))
        for(j in 1:length(test.no)){
          test.no.now <- test.no[j]
          if(test.no.now == 1){
            score.now <- score.linker.cpp(y, Ws = Ws.A, Gammas = Gammas.A,
                                          gammas.diag = TRUE, Gu = Gu, Ge = Ge,
                                          P0 = P0, chi0.mixture = chi0.mixture)
          }
          
          if(test.no.now == 2){
            if(any(MAF.cut.D)){
              score.now <- score.linker.cpp(y, Ws = Ws.D, Gammas = Gammas.D,
                                            gammas.diag = TRUE, Gu = Gu, Ge = Ge,
                                            P0 = P0, chi0.mixture = chi0.mixture)
            }else{
              score.now <- 0
            }
          }
          
          if(test.no.now == 3){
            if(any(MAF.cut.D)){
              score.now <- score.linker.cpp(y, Ws = Ws.AD, Gammas = Gammas.AD,
                                            gammas.diag = TRUE, Gu = Gu, Ge = Ge,
                                            P0 = P0, chi0.mixture = chi0.mixture)
            }else{
              score.now <- 0
            }
          }
          
          scores.now[j] <- score.now
        }
      }
      scores[i, ] <- scores.now
    }
  }
  
  if(is.null(gene.set)){
    rownames(scores) <- window.centers
  }else{
    rownames(scores) <- gene.name
  }
  
  if(kernel.method == "linear"){
    colnames(scores) <- test.effect
  }else{
    colnames(scores) <- kernel.method
  }
  
  if (count) {
    cat("\n")
  }
  return(scores)
}





#' Calculate -log10(p) of each SNP-set by the score test (multi-cores)
#'
#' @description This function calculates -log10(p) of each SNP-set by the score test.
#' First, the function calculates the score statistic
#' without solving the multi-kernel mixed model for each SNP-set.
#' Then it performs the score test by using the fact that the score statistic follows the chi-square distribution.
#'
#'
#' @param M.now A n.sample x n.mark genotype matrix where n.sample is sample size and n.mark is the number of markers.
#' @param ZETA.now A list of variance (relationship) matrix (K; \eqn{m \times m}) and its design matrix (Z; \eqn{n \times m}) of random effects. You can use only one kernel matrix.
#' For example, ZETA = list(A = list(Z = Z, K = K))
#' Please set names of list "Z" and "K"!
#' @param y A \eqn{n \times 1} vector. A vector of phenotypic values should be used. NA is allowed.
#' @param X.now A \eqn{n \times p} matrix. You should assign mean vector (rep(1, n)) and covariates. NA is not allowed.
#' @param LL0 The log-likelihood for the null model.
#' @param Gu A \eqn{n \times n} matrix. You should assign \eqn{ZKZ'}, where K is covariance (relationship) matrix and Z is its design matrix.
#' @param Ge A \eqn{n \times n} matrix. You should assign identity matrix I (diag(n)).
#' @param P0 A \eqn{n \times n} matrix. The Moore-Penrose generalized inverse of \eqn{SV0S}, where \eqn{S = X(X'X)^{-1}X'} and
#' \eqn{V0 = \sigma^2_u Gu + \sigma^2_e Ge}. \eqn{\sigma^2_u} and \eqn{\sigma^2_e} are estimators of the null model.
#' @param n.core Setting n.core > 1 will enable parallel execution on a machine with multiple cores.
#' @param map Data frame of map information where the first column is the marker names,
#' the second and third column is the chromosome amd map position, and the forth column is -log10(p) for each marker.
#' @param kernel.method It determines how to calculate kernel. There are three methods.
#' \describe{
#' \item{"gaussian"}{It is the default method. Gaussian kernel is calculated by distance matrix.}
#' \item{"exponential"}{When this method is selected, exponential kernel is calculated by distance matrix.}
#' \item{"linear"}{When this method is selected, linear kernel is calculated by NOIA methods for additive GRM.}
#'}
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
#' @param chi0.mixture RAINBOWR assumes the test statistic \eqn{l1' F l1} is considered to follow a x chisq(df = 0) + (1 - a) x chisq(df = r).
#' where l1 is the first derivative of the log-likelihood and F is the Fisher information. And r is the degree of freedom.
#' The argument chi0.mixture is a (0 <= a < 1), and default is 0.5.
#' @param weighting.center In kernel-based GWAS, weights according to the Gaussian distribution (centered on the tested SNP) are taken into account when calculating the kernel if Rainbow = TRUE.
#'           If weighting.center = FALSE, weights are not taken into account.
#' @param weighting.other You can set other weights in addition to weighting.center. The length of this argument should be equal to the number of SNPs.
#'           For example, you can assign SNP effects from the information of gene annotation.
#' @param gene.set If you have information of gene, you can use it to perform kernel-based GWAS.
#'            You should assign your gene information to gene.set in the form of a "data.frame" (whose dimension is (the number of gene) x 2).
#'            In the first column, you should assign the gene name. And in the second column, you should assign the names of each marker,
#'            which correspond to the marker names of "geno" argument.
#' @param min.MAF Specifies the minimum minor allele frequency (MAF).
#' If a marker has a MAF less than min.MAF, it is assigned a zero score.
#' @param count When count is TRUE, you can know how far RGWAS has ended with percent display.
#'
#' @return -log10(p) for each SNP-set
#'
#' @references Listgarten, J. et al. (2013) A powerful and efficient set test
#'  for genetic markers that handles confounders. Bioinformatics. 29(12): 1526-1533.
#'
#' Lippert, C. et al. (2014) Greater power and computational efficiency for kernel-based
#'  association testing of sets of genetic variants. Bioinformatics. 30(22): 3206-3214.
#'
#'
#'
#'
score.calc.score.MC <- function(M.now, y, X.now, ZETA.now, LL0, Gu, Ge, P0, n.core = 2,
                                map, kernel.method = "linear", kernel.h = "tuned", haplotype = TRUE, num.hap = NULL,
                                test.effect = "additive", window.size.half = 5, window.slide = 1,
                                chi0.mixture = 0.5, weighting.center = TRUE, weighting.other = NULL,
                                gene.set = NULL, min.MAF = 0.02, count = TRUE){
  chr <- map[, 2]
  chr.tab <- table(chr)
  chr.max <- max(chr)
  chr.cum <- cumsum(chr.tab)
  n.scores.each <- (chr.tab + (window.slide - 1)) %/% window.slide
  cum.n.scores <- cumsum(n.scores.each)
  if(is.null(gene.set)){
    n.scores <- sum(n.scores.each)
  }else{
    gene.names <- as.character(gene.set[, 1])
    mark.id <- as.character(gene.set[, 2])
    gene.name <- as.character(unique(gene.names))
    n.scores <- length(unique(gene.set[, 1]))
  }
  
  if(kernel.method == "linear"){
    ncol.scores <- length(test.effect)
  }else{
    ncol.scores <- 1
  }
  
  window.centers <- rep(NA, n.scores)
  probaa <- apply(M.now == -1, 2, mean)
  probAa <- apply(M.now == 0, 2, mean)
  freq <- probaa + probAa / 2
  MAF <- pmin(freq, 1 - freq)
  MAF.D <- pmin(probAa, 1 - probAa)
  
  score.calc.score.MC.oneSNP <- function(markNo) {
    if(is.null(gene.set)){
      markNo.chr <- min(which(markNo - cum.n.scores <= 0))
      if(markNo.chr >= 2){
        window.center <- window.slide * (markNo - cum.n.scores[markNo.chr - 1] - 1) + chr.cum[markNo.chr - 1] + 1
      }else{
        window.center <- window.slide * (markNo - 1) + 1
      }
      names(window.center) <- markNo.chr
      Theories1 <-  window.center < window.size.half + 1
      for(r in 1:(chr.max - 1)){
        Theory1 <- chr.cum[r] < window.center & window.center < window.size.half + 1 + chr.cum[r]
        Theories1 <- c(Theories1, Theory1)
      }
      rule1 <- sum(Theories1) != 0
      
      
      Theories2  <- NULL
      for(r in 1:chr.max){
        Theory2 <- chr.cum[r] - (window.size.half + 1) < window.center & window.center <= chr.cum[r]
        Theories2 <- c(Theories2, Theory2)
      }
      rule2 <- sum(Theories2) != 0
      
      
      if(rule1 & rule2){
        Mis.range.0 <- which(chr == markNo.chr)
        Mis.range.02 <- which(chr == markNo.chr) - window.center + 1 + window.size.half
      }else{
        if(rule1){
          near.min <- c(0, chr.cum)[which.min(abs(window.center - c(0, chr.cum)))]
          Mis.range.0 <- (near.min + 1):(window.center + window.size.half)
          Mis.range.02 <- (2 * window.size.half + 2 - length(Mis.range.0)):(2 * window.size.half + 1)
        }else{
          if(rule2){
            near.max <- c(0, chr.cum)[which.min(abs(window.center - c(0, chr.cum)))]
            Mis.range.0 <- (window.center - window.size.half):near.max
            Mis.range.02 <- 1:length(Mis.range.0)
          }else{
            Mis.range.0 <- (window.center - window.size.half):(window.center + window.size.half)
            Mis.range.02 <- 1:(2 * window.size.half + 1)
          }
        }
      }
    }else{
      mark.name.now <- mark.id[gene.names == gene.name[markNo]]
      Mis.range.0 <- match(mark.name.now, map[, 1])
      Mis.range.02 <- 1:length(Mis.range.0)
      weighting.center <- FALSE
    }
    
    Mis.0 <- M.now[, Mis.range.0, drop = FALSE]
    MAF.cut <- MAF[Mis.range.0] >= min.MAF
    if(any(test.effect %in% c("dominance", "additive+dominance"))){
      Mis.D.0 <- M.now[, Mis.range.0, drop = FALSE]
      MAF.cut.D <- MAF.D[Mis.range.0] > 0
    }else{
      MAF.cut.D <- rep(TRUE, length(MAF.cut))
    }
    
    if(any(MAF.cut)){
      Mis.0 <- Mis.0[, MAF.cut, drop = FALSE]
      Mis.range <- Mis.range.0[MAF.cut]
      Mis.range2 <- Mis.range.02[MAF.cut]
      window.size <- ncol(Mis.0)
      if(any(MAF.cut.D)){
        if(any(test.effect %in% c("dominance", "additive+dominance"))){
          Mis.D.0 <- Mis.D.0[, MAF.cut.D, drop = FALSE]
          Mis.range.D <- Mis.range.0[MAF.cut.D]
          Mis.range2.D <- Mis.range.02[MAF.cut.D]
          window.size.D <- ncol(Mis.D.0)
        }
      }
      
      if(haplotype){
        if(is.null(num.hap)){
          Mis.fac <- factor(apply(Mis.0, 1, function(x) paste(x, collapse = "")))
          Mis <- Mis.0[!duplicated(as.numeric(Mis.fac)), , drop = FALSE]
          
          bango <- as.factor(as.numeric(Mis.fac))
          levels(bango) <- order(unique(bango))
          bango <- as.numeric(as.character(bango))
          if(any(MAF.cut.D)){
            if(any(test.effect %in% c("dominance", "additive+dominance"))){
              Mis.D.fac <- factor(apply(Mis.D.0, 1, function(x) paste(x, collapse = "")))
              Mis.D <- Mis.D.0[!duplicated(as.numeric(Mis.D.fac)), , drop = FALSE]
              
              bango.D <- as.factor(as.numeric(Mis.D.fac))
              levels(bango.D) <- order(unique(bango.D))
              bango.D <- as.numeric(as.character(bango.D))
            }
          }
        }else{
          kmed.res <- cluster::pam(Mis.0, k = num.hap)
          Mis <- kmed.res$medoids
          bango <- kmed.res$clustering
          if(any(MAF.cut.D)){
            if(any(test.effect %in% c("dominance", "additive+dominance"))){
              kmed.res.D <- cluster::pam(Mis.D.0, k = num.hap)
              Mis.D <- kmed.res.D$medoids
              bango.D <- kmed.res.D$clustering
            }
          }
        }
        Z.part <- as.matrix(Matrix::sparseMatrix(i = 1:nrow(M.now), j = bango, x = rep(1, nrow(M.now)),
                                                 dims = c(nrow(M.now), nrow(Mis))))
        if(any(MAF.cut.D)){
          if(any(test.effect %in% c("dominance", "additive+dominance"))){
            Z.part.D <- as.matrix(Matrix::sparseMatrix(i = 1:nrow(M.now), j = bango.D, x = rep(1, nrow(M.now)),
                                                       dims = c(nrow(M.now), nrow(Mis.D))))
          }
        }
      }else{
        Mis <- Mis.0
        Mis.D <- Mis.D.0
        Z.part <- Z.part.D <- diag(nrow(M.now))
      }
      
      if(window.size != 1){
        if(weighting.center){
          weight.Mis <- dnorm((-window.size.half):(window.size.half), 0, window.size.half / 2)[Mis.range2]
          weight.Mis <- weight.Mis / apply(Mis, 2, sd)
          if(!is.null(weighting.other)){
            weight.Mis <- weight.Mis * weighting.other[Mis.range]
          }
          weight.Mis <- weight.Mis * window.size / sum(weight.Mis)
        }else{
          weight.Mis <- rep(1, window.size)
          weight.Mis <- weight.Mis / apply(Mis, 2, sd)
          if(!is.null(weighting.other)){
            weight.Mis <- weight.Mis * weighting.other[Mis.range]
          }
          weight.Mis <- weight.Mis * window.size / sum(weight.Mis)
        }
      }else{
        weight.Mis <- 1
      }
      
      if(any(MAF.cut.D)){
        if(window.size != 1){
          if(any(test.effect %in% c("dominance", "additive+dominance"))){
            if(weighting.center){
              weight.Mis.D <- dnorm((-window.size.half):(window.size.half), 0, window.size.half / 2)[Mis.range2.D]
              weight.Mis.D <- weight.Mis.D / apply(Mis.D, 2, sd)
              if(!is.null(weighting.other)){
                weight.Mis.D <- weight.Mis.D * weighting.other[Mis.range.D]
              }
              weight.Mis.D <- weight.Mis.D * window.size.D / sum(weight.Mis.D)
            }else{
              weight.Mis.D <- rep(1, window.size.D)
              weight.Mis.D <- weight.Mis.D / apply(Mis.D, 2, sd)
              if(!is.null(weighting.other)){
                weight.Mis.D <- weight.Mis.D * weighting.other[Mis.range.D]
              }
              weight.Mis.D <- weight.Mis.D * window.size.D / sum(weight.Mis.D)
            }
          }
        }else{
          weight.Mis.D <- 1
        }
      }
      
      if(kernel.method != "linear"){
        if(ncol(Mis) != 1){
          Mis.weighted <- t(apply(Mis, 1, function(x) x * weight.Mis))
        }else{
          Mis.weighted <- as.matrix(apply(Mis, 1, function(x) x * weight.Mis))
        }
        
        
        K.SNP <- calcGRM(genoMat = Mis.weighted,
                         methodGRM = kernel.method,
                         kernel.h = kernel.h,
                         returnWMat = FALSE)
        
        Ws <- list(W = Z.part)
        Gammas <- list(Gamma = K.SNP)
        scores.now <- score.linker.cpp(y, Ws = Ws, Gammas = Gammas,
                                       gammas.diag = FALSE, Gu = Gu, Ge = Ge,
                                       P0 = P0, chi0.mixture = chi0.mixture)
      }else{
        test.no <- match(test.effect, c("additive", "dominance", "additive+dominance"))
        if(length(test.no) == 0){
          stop("The effect to test should be 'additive', 'dominance' or 'additive+dominance'!")
        }
        
        if(any(test.effect %in% c("additive", "additive+dominance"))){
          W.A <- calcGRM(genoMat = Mis,
                         methodGRM = "addNOIA",
                         returnWMat = TRUE,
                         probaa = probaa[Mis.range],
                         probAa = probAa[Mis.range])
        }
        if(any(MAF.cut.D)){
          if(any(test.effect %in% c("dominance", "additive+dominance"))){
            W.D <- calcGRM(genoMat = Mis.D,
                           methodGRM = "domNOIA",
                           returnWMat = TRUE,
                           probaa = probaa[Mis.range.D],
                           probAa = probAa[Mis.range.D])
          }
        }
        
        if(1 %in% test.no){
          Ws.A <- list(W.A = Z.part %*% W.A)
          Gammas.A <- list(W.A = diag(weight.Mis ^ 2))
        }
        
        if(any(MAF.cut.D)){
          if(2 %in% test.no){
            Ws.D <- list(W.D = Z.part.D %*% W.D)
            Gammas.D <- list(W.D = diag(weight.Mis.D ^ 2))
          }
          
          if(3 %in% test.no){
            Ws.AD <- list(W.A = Z.part %*% W.A, Z.part.D %*% W.D)
            Gammas.AD <- list(W.A = diag(weight.Mis ^ 2), W.D = diag(weight.Mis.D ^ 2))
          }
        }
        
        scores.now <- rep(NA, length(test.no))
        for(j in 1:length(test.no)){
          test.no.now <- test.no[j]
          if(test.no.now == 1){
            score.now <- score.linker.cpp(y, Ws = Ws.A, Gammas = Gammas.A,
                                          gammas.diag = TRUE, Gu = Gu, Ge = Ge,
                                          P0 = P0, chi0.mixture = chi0.mixture)
          }
          
          if(test.no.now == 2){
            if(any(MAF.cut.D)){
              score.now <- score.linker.cpp(y, Ws = Ws.D, Gammas = Gammas.D,
                                            gammas.diag = TRUE, Gu = Gu, Ge = Ge,
                                            P0 = P0, chi0.mixture = chi0.mixture)
            }else{
              score.now <- 0
            }
          }
          
          if(test.no.now == 3){
            if(any(MAF.cut.D)){
              score.now <- score.linker.cpp(y, Ws = Ws.AD, Gammas = Gammas.AD,
                                            gammas.diag = TRUE, Gu = Gu, Ge = Ge,
                                            P0 = P0, chi0.mixture = chi0.mixture)
            }else{
              score.now <- 0
            }
          }
          
          scores.now[j] <- score.now
        }
      }
    } else {
      scores.now <- rep(NA, ncol.scores)
    }
    
    if(is.null(gene.set)){
      return(list(scores = scores.now, window.center = window.center))
    } else {
      return(list(scores = scores.now))
    }
  }
  
  all.res <- pbmcapply::pbmclapply(1:n.scores, score.calc.score.MC.oneSNP, mc.cores = n.core)
  scores <- unlist(lapply(all.res, function(x) x$scores))
  scores <- matrix(scores, nrow = n.scores, ncol = ncol.scores, byrow = TRUE)
  
  
  if(is.null(gene.set)){
    window.centers <- unlist(lapply(all.res, function(x) x$window.center))
    rownames(scores) <- window.centers
  }else{
    rownames(scores) <- gene.name
  }
  
  if(kernel.method == "linear"){
    colnames(scores) <- test.effect
  }else{
    colnames(scores) <- kernel.method
  }
  
  if (count) {
    cat("\n")
  }
  return(scores)
}






#' Calculate -log10(p) of epistatic effects by LR test
#'
#'
#' @param M.now A n.sample x n.mark genotype matrix where n.sample is sample size and n.mark is the number of markers.
#' @param ZETA.now A list of variance (relationship) matrix (K; \eqn{m \times m}) and its design matrix (Z; \eqn{n \times m}) of random effects. You can use only one kernel matrix.
#' For example, ZETA = list(A = list(Z = Z, K = K))
#' Please set names of list "Z" and "K"!
#' @param y A \eqn{n \times 1} vector. A vector of phenotypic values should be used. NA is allowed.
#' @param X.now A \eqn{n \times p} matrix. You should assign mean vector (rep(1, n)) and covariates. NA is not allowed.
#' @param eigen.SGS A list with
#' \describe{
#' \item{$values}{Eigen values}
#' \item{$vectors}{Eigen vectors}
#' }
#' The result of the eigen decompsition of \eqn{SGS}, where \eqn{S = I - X(X'X)^{-1}X'}, \eqn{G = ZKZ'}.
#' You can use "spectralG.cpp" function in RAINBOWR.
#' If this argument is NULL, the eigen decomposition will be performed in this function.
#' We recommend you assign the result of the eigen decomposition beforehand for time saving.
#' @param eigen.G A list with
#' \describe{
#' \item{$values}{Eigen values}
#' \item{$vectors}{Eigen vectors}
#' }
#' The result of the eigen decompsition of \eqn{G = ZKZ'}. You can use "spectralG.cpp" function in RAINBOWR.
#' If this argument is NULL, the eigen decomposition will be performed in this function.
#' We recommend you assign the result of the eigen decomposition beforehand for time saving.
#' @param map Data frame of map information where the first column is the marker names,
#' the second and third column is the chromosome amd map position, and the forth column is -log10(p) for each marker.
#' @param haplotype If the number of lines of your data is large (maybe > 100), you should set haplotype = TRUE.
#'             When haplotype = TRUE, haplotype-based kernel will be used for calculating -log10(p).
#'             (So the dimension of this gram matrix will be smaller.)
#'             The result won't be changed, but the time for the calculation will be shorter.
#' @param num.hap When haplotype = TRUE, you can set the number of haplotypes which you expect.
#'           Then similar arrays are considered as the same haplotype, and then make kernel(K.SNP) whose dimension is num.hap x num.hap.
#'           When num.hap = NULL (default), num.hap will be set as the maximum number which reflects the difference between lines.
#' @param optimizer The function used in the optimization process. We offer "optim", "optimx", and "nlminb" functions.
#' @param window.size.half This argument decides how many SNPs (around the SNP you want to test) are used to calculated K.SNP.
#' More precisely, the number of SNPs will be 2 * window.size.half + 1.
#' @param window.slide This argument determines how often you test markers. If window.slide = 1, every marker will be tested.
#' If you want to perform SNP set by bins, please set window.slide = 2 * window.size.half + 1.
#' @param chi0.mixture RAINBOWR assumes the tdeviance is considered to follow a x chisq(df = 0) + (1 - a) x chisq(df = r).
#' where r is the degree of freedom.
#' The argument chi0.mixture is a (0 <= a < 1), and default is 0.5.
#' @param gene.set If you have information of gene, you can use it to perform kernel-based GWAS.
#'            You should assign your gene information to gene.set in the form of a "data.frame" (whose dimension is (the number of gene) x 2).
#'            In the first column, you should assign the gene name. And in the second column, you should assign the names of each marker,
#'            which correspond to the marker names of "geno" argument.
#' @param dominance.eff If this argument is TRUE, dominance effect is included in the model,
#' and additive x dominance and dominance x dominance are also tested as epistatic effects.
#' When you use inbred lines, please set this argument FALSE.
#' @param min.MAF Specifies the minimum minor allele frequency (MAF).
#' If a marker has a MAF less than min.MAF, it is assigned a zero score.
#' @param count When count is TRUE, you can know how far RGWAS has ended with percent display.
#'
#' @return -log10(p) of epistatic effects for each SNP-set
#'
#' @references Listgarten, J. et al. (2013) A powerful and efficient set test
#'  for genetic markers that handles confounders. Bioinformatics. 29(12): 1526-1533.
#'
#' Lippert, C. et al. (2014) Greater power and computational efficiency for kernel-based
#'  association testing of sets of genetic variants. Bioinformatics. 30(22): 3206-3214.
#'
#' Jiang, Y. and Reif, J.C. (2015) Modeling epistasis in genomic selection. Genetics. 201(2): 759-768.
#'
#'
#'
score.calc.epistasis.LR <- function(M.now, y, X.now, ZETA.now, eigen.SGS = NULL, eigen.G = NULL, optimizer = "nlminb",
                                    map, haplotype = TRUE, num.hap = NULL, window.size.half = 5, window.slide = 1,
                                    chi0.mixture = 0.5, gene.set = NULL, dominance.eff = TRUE,
                                    min.MAF = 0.02, count = TRUE){
  chr <- map[, 2]
  chr.tab <- table(chr)
  chr.max <- max(chr)
  chr.cum <- cumsum(chr.tab)
  n.scores.each <- (chr.tab + (window.slide - 1)) %/% window.slide
  cum.n.scores <- cumsum(n.scores.each)
  if(is.null(gene.set)){
    n.scores <- sum(n.scores.each)
  }else{
    gene.names <- as.character(gene.set[, 1])
    mark.id <- as.character(gene.set[, 2])
    gene.name <- as.character(unique(gene.names))
    n.scores <- length(unique(gene.set[, 1]))
  }
  
  scores <- matrix(0, nrow = n.scores, ncol = n.scores)
  window.centers <- rep(NA, n.scores)
  probaa <- apply(M.now == -1, 2, mean)
  probAa <- apply(M.now == 0, 2, mean)
  freq <- probaa + probAa / 2
  MAF <- pmin(freq, 1 - freq)
  MAF.D <- pmin(probAa, 1 - probAa)
  
  
  n.sample.now <- nrow(M.now)
  Z.normal <- diag(n.sample.now)
  W.A.list <- W.A.0.list  <- W.D.list <- W.D.0.list <-
    Z.A.part.list <- Z.D.part.list <- rep(list(NA), n.scores)
  pb <- txtProgressBar(min = 1, max = n.scores, style = 3)
  n.scores2 <- n.scores - n.scores %% 100
  start.scorecalc <- Sys.time()
  
  for(i in 1:n.scores){
    if(is.null(gene.set)){
      i.chr <- min(which(i - cum.n.scores <= 0))
      if(i.chr >= 2){
        window.center <- window.slide * (i - cum.n.scores[i.chr - 1] - 1) + chr.cum[i.chr - 1] + 1
      }else{
        window.center <- window.slide * (i - 1) + 1
      }
      names(window.center) <- i.chr
      window.centers[i] <- window.center
      Theories1 <-  window.center < window.size.half + 1
      for(r in 1:(chr.max - 1)){
        Theory1 <- chr.cum[r] < window.center & window.center < window.size.half + 1 + chr.cum[r]
        Theories1 <- c(Theories1, Theory1)
      }
      rule1 <- sum(Theories1) != 0
      
      
      Theories2  <- NULL
      for(r in 1:chr.max){
        Theory2 <- chr.cum[r] - (window.size.half + 1) < window.center & window.center <= chr.cum[r]
        Theories2 <- c(Theories2, Theory2)
      }
      rule2 <- sum(Theories2) != 0
      
      
      
      if(rule1 & rule2){
        Mis.range.0 <- which(chr == i.chr)
        Mis.range.02 <- which(chr == i.chr) - window.center + 1 + window.size.half
      }else{
        if(rule1){
          near.min <- c(0, chr.cum)[which.min(abs(window.center - c(0, chr.cum)))]
          Mis.range.0 <- (near.min + 1):(window.center + window.size.half)
          Mis.range.02 <- (2 * window.size.half + 2 - length(Mis.range.0)):(2 * window.size.half + 1)
        }else{
          if(rule2){
            near.max <- c(0, chr.cum)[which.min(abs(window.center - c(0, chr.cum)))]
            Mis.range.0 <- (window.center - window.size.half):near.max
            Mis.range.02 <- 1:length(Mis.range.0)
          }else{
            Mis.range.0 <- (window.center - window.size.half):(window.center + window.size.half)
            Mis.range.02 <- 1:(2 * window.size.half + 1)
          }
        }
      }
    }else{
      mark.name.now <- mark.id[gene.names == gene.name[i]]
      Mis.range.0 <- match(mark.name.now, map[, 1])
      Mis.range.02 <- 1:length(Mis.range.0)
      weighting.center <- FALSE
    }
    
    Mis.0.0 <- M.now[, Mis.range.0, drop = FALSE]
    MAF.cut <- MAF[Mis.range.0] >= min.MAF
    if(dominance.eff){
      Mis.D.0.0 <- M.now[, Mis.range.0, drop = FALSE]
      MAF.cut.D <- MAF.D[Mis.range.0] > 0
    }else{
      MAF.cut.D <- rep(TRUE, length(MAF.cut))
    }
    
    
    if(any(MAF.cut)){
      Mis.0 <- Mis.0.0[, MAF.cut, drop = FALSE]
      Mis.range <- Mis.range.0[MAF.cut]
      Mis.range2 <- Mis.range.02[MAF.cut]
      window.size <- ncol(Mis.0)
      if(any(MAF.cut.D)){
        if(dominance.eff){
          Mis.D.0 <- Mis.D.0.0[, MAF.cut.D, drop = FALSE]
          Mis.range.D <- Mis.range.0[MAF.cut.D]
          Mis.range2.D <- Mis.range.02[MAF.cut.D]
          window.size.D <- ncol(Mis.D.0)
        }
      }
      
      if(haplotype){
        if(is.null(num.hap)){
          Mis.fac <- factor(apply(Mis.0, 1, function(x) paste(x, collapse = "")))
          Mis <- Mis.0[!duplicated(as.numeric(Mis.fac)), , drop = FALSE]
          
          bango <- as.factor(as.numeric(Mis.fac))
          levels(bango) <- order(unique(bango))
          bango <- as.numeric(as.character(bango))
          if(any(MAF.cut.D)){
            if(dominance.eff){
              Mis.D.fac <- factor(apply(Mis.D.0, 1, function(x) paste(x, collapse = "")))
              Mis.D <- Mis.D.0[!duplicated(as.numeric(Mis.D.fac)), , drop = FALSE]
              
              bango.D <- as.factor(as.numeric(Mis.D.fac))
              levels(bango.D) <- order(unique(bango.D))
              bango.D <- as.numeric(as.character(bango.D))
            }
          }
        }else{
          kmed.res <- cluster::pam(Mis.0, k = num.hap)
          Mis <- kmed.res$medoids
          bango <- kmed.res$clustering
          if(any(MAF.cut.D)){
            if(dominance.eff){
              kmed.res.D <- cluster::pam(Mis.D.0, k = num.hap)
              Mis.D <- kmed.res.D$medoids
              bango.D <- kmed.res.D$clustering
            }
          }
        }
        Z.part <- as.matrix(Matrix::sparseMatrix(i = 1:nrow(M.now), j = bango, x = rep(1, nrow(M.now)),
                                                 dims = c(nrow(M.now), nrow(Mis))))
        if(any(MAF.cut.D)){
          if(dominance.eff){
            Z.part.D <- as.matrix(Matrix::sparseMatrix(i = 1:nrow(M.now), j = bango.D, x = rep(1, nrow(M.now)),
                                                       dims = c(nrow(M.now), nrow(Mis.D))))
          }
        }
      }else{
        Mis <- Mis.0
        Mis.D <- Mis.D.0        
        Z.part <- Z.part.D <- diag(nrow(M.now))
      }
      
      
      W.A <- calcGRM(genoMat = Mis,
                     methodGRM = "addNOIA",
                     returnWMat = TRUE,
                     probaa = probaa[Mis.range],
                     probAa = probAa[Mis.range])
      
      W.A.0 <- calcGRM(genoMat = Mis.0.0,
                       methodGRM = "addNOIA",
                       returnWMat = TRUE,
                       probaa = probaa[Mis.range.0],
                       probAa = probAa[Mis.range.0])
      
      W.A.list[[i]] <- W.A
      Z.A.part.list[[i]] <- Z.part
      W.A.0.list[[i]] <- W.A.0
      
      if(any(MAF.cut.D)){
        if(dominance.eff){
          W.D <- calcGRM(genoMat = Mis.D,
                         methodGRM = "domNOIA",
                         returnWMat = TRUE,
                         probaa = probaa[Mis.range.D],
                         probAa = probAa[Mis.range.D])
          
          W.D.0 <- calcGRM(genoMat = Mis.0.0,
                           methodGRM = "domNOIA",
                           returnWMat = TRUE,
                           probaa = probaa[Mis.range.0],
                           probAa = probAa[Mis.range.0])
          
          W.D.list[[i]] <- W.D
          Z.D.part.list[[i]] <- Z.part.D
          W.D.0.list[[i]] <- W.D.0
        }
      }
    }
  }
  
  n.calc <- n.scores * (n.scores + 1) / 2
  pb <- txtProgressBar(min = 1, max = n.calc, style = 3)
  n.calc2 <- n.calc - n.calc %% 100
  start.scorecalc <- Sys.time()
  
  
  for(i in 1:n.scores){
    W.A.1 <- W.A.list[[i]]
    Z.A.1.part <- Z.A.part.list[[i]]
    W.A.0.1 <- W.A.0.list[[i]]
    m.A.1 <- ncol(W.A.0.1)
    
    if(dominance.eff){
      W.D.1 <- W.D.list[[i]]
      Z.D.1.part <- Z.D.part.list[[i]]
      W.D.0.1 <- W.D.0.list[[i]]
      isna.1 <- any(is.na(W.D.0.1))
      m.D.1 <- ncol(W.D.0.1)
    }else{
      isna.1 <- TRUE
    }
    for(j in i:n.scores){
      num.now <- (i - 1) * n.scores + (j - i + 1)
      if(count){
        setTxtProgressBar(pb, num.now)
        if(n.calc > 100){
          if(num.now == (n.calc2 / 100 + 1) | num.now == (n.calc2 / 10 + 1) | num.now == (n.calc2 / 2 + 1)){
            cat("\n")
            end.scorecalc <- Sys.time()
            jikan.scorecalc <- (end.scorecalc - start.scorecalc) * (n.calc - num.now + 1) / (num.now - 1)
            print(paste0((num.now - 1) * 100 / n.calc2, "%...Done. ",
                         round(jikan.scorecalc, 2), " ", attr(jikan.scorecalc, "units"),
                         " to end.  Scheduled end time : ", end.scorecalc + jikan.scorecalc))
          }
        }
      }
      
      if(i == j){
        if((!dominance.eff) | isna.1){
          W.AA <- W.A.0.1 ^ 2
          W.AA <- W.AA * sqrt(nrow(W.AA) / sum(W.AA * W.AA))
          
          Ws0.null <- list(W.A = W.A.1)
          Ws0.alt <- list(W.A = W.A.1, W.AA = W.AA)
          
          Zs0.null <- list(W.A = Z.A.1.part)
          Zs0.alt <- list(W.A = Z.A.1.part, W.AA = Z.normal)
          
          lin.method <- TRUE
          df <- 1
        }else{
          if(m.A.1 == m.D.1){
            W.AA <- W.A.0.1 ^ 2
            W.AA <- W.AA * sqrt(nrow(W.AA) / sum(W.AA * W.AA))
            W.AD <- W.A.0.1 * W.D.0.1
            W.AD <- W.AD * sqrt(nrow(W.AD) / sum(W.AD * W.AD))
            W.DD <- W.D.0.1 ^ 2
            W.DD <- W.DD * sqrt(nrow(W.DD) / sum(W.DD * W.DD))
            
            Ws0.null <- list(W.A = W.A.1, W.D = W.D.1)
            Ws0.alt <- list(W.A = W.A.1, W.D = W.D.1, W.AA = W.AA, W.AD = W.AD, W.DD = W.DD)
            Zs0.null <- list(W.A = Z.A.1.part, W.D = Z.D.1.part)
            Zs0.alt <- list(W.A = Z.A.1.part, W.D = Z.D.1.part, W.AA = Z.normal,
                            W.AD = Z.normal, W.DD = Z.normal)
            
            lin.method <- TRUE
          }else{
            K.A.1.part <- tcrossprod(W.A.1)
            K.D.1.part <- tcrossprod(W.D.1)
            
            K.A.0.1.part <- tcrossprod(W.A.0.1)
            K.D.0.1.part <- tcrossprod(W.D.0.1)
            K.AA.part <- K.A.0.1.part ^ 2
            K.AA.part <- K.AA.part * sqrt(nrow(K.AA.part) / sum(K.AA.part * K.AA.part))
            K.AD.part <- K.A.0.1.part * K.D.0.1.part
            K.AD.part <- K.AD.part * sqrt(nrow(K.AD.part) / sum(K.AD.part * K.AD.part))
            K.DD.part <- K.D.0.1.part ^ 2
            K.DD.part <- K.DD.part * sqrt(nrow(K.DD.part) / sum(K.DD.part * K.DD.part))
            
            ZETA.now2.null <- c(ZETA.now, list(A.part = list(Z = Z.A.1.part, K = K.A.1.part)),
                                list(D.part = list(Z = Z.D.1.part, K = K.D.1.part)))
            
            ZETA.now2.alt <- c(ZETA.now, list(A.part = list(Z = Z.A.1.part, K = K.A.1.part)),
                               list(D.part = list(Z = Z.D.1.part, K = K.D.1.part)),
                               list(AA.part = list(Z = Z.normal, K = K.AA.part)),
                               list(AD.part = list(Z = Z.normal, K = K.AD.part)),
                               list(DD.part = list(Z = Z.normal, K = K.DD.part)))
            
            lin.method <- FALSE
          }
          df <- 3
        }
      }else{
        W.A.2 <- W.A.list[[j]]
        Z.A.2.part <- Z.A.part.list[[j]]
        W.A.0.2 <- W.A.0.list[[j]]
        m.A.2 <- ncol(W.A.0.2)
        
        if(dominance.eff){
          W.D.2 <- W.D.list[[j]]
          Z.D.2.part <- Z.D.part.list[[j]]
          W.D.0.2 <- W.D.0.list[[j]]
          isna.2 <- any(is.na(W.D.0.2))
          m.D.2 <- ncol(W.D.0.2)
          
          isnas <- c(isna.1, isna.2)
        }else{
          isna.2 <- TRUE
          isnas <- c(isna.1, isna.2)
        }
        
        if((!dominance.eff) | all(isnas)){
          if(m.A.1 == m.A.2){
            W.AA <- W.A.0.1 * W.A.0.2
            W.AA <- W.AA * sqrt(nrow(W.AA) / sum(W.AA * W.AA))
            
            Ws0.null <- list(W.A.1 = W.A.1, W.A.2 = W.A.2)
            Ws0.alt <- list(W.A.1 = W.A.1, W.A.2 = W.A.2, W.AA = W.AA)
            
            Zs0.null <- list(W.A.1 = Z.A.1.part, W.A.2 = Z.A.2.part)
            Zs0.alt <- list(W.A.1 = Z.A.1.part, W.A.2 = Z.A.2.part, W.AA = Z.normal)
            
            lin.method <- TRUE
          }else{
            K.A.1.part <- tcrossprod(W.A.1)
            K.A.2.part <- tcrossprod(W.A.2)
            
            K.A.0.1.part <- tcrossprod(W.A.0.1)
            K.A.0.2.part <- tcrossprod(W.A.0.2)
            K.AA.part <- K.A.0.1.part * K.A.0.2.part
            K.AA.part <- K.AA.part * sqrt(nrow(K.AA.part) / sum(K.AA.part * K.AA.part))
            
            ZETA.now2.null <- c(ZETA.now, list(A.1.part = list(Z = Z.A.1.part, K = K.A.1.part)),
                                list(A.2.part = list(Z = Z.A.2.part, K = K.A.2.part)))
            
            ZETA.now2.alt <- c(ZETA.now, list(A.1.part = list(Z = Z.A.1.part, K = K.A.1.part)),
                               list(A.2.part = list(Z = Z.A.2.part, K = K.A.2.part)),
                               list(AA.part = list(Z = Z.normal, K = K.AA.part)))
            
            lin.method <- FALSE
          }
          df <- 1
        }else{
          if(isna.1){
            if(all(c(m.A.2, m.D.2) == m.A.1)){
              W.AA <- W.A.0.1 * W.A.0.2
              W.AD <- W.A.0.1 * W.D.0.2
              W.AA <- W.AA * sqrt(nrow(W.AA) / sum(W.AA * W.AA))
              W.AD <- W.AD * sqrt(nrow(W.AD) / sum(W.AD * W.AD))
              
              Ws0.null <- list(W.A.1 = W.A.1, W.A.2 = W.A.2, W.D.2 = W.D.2)
              Ws0.alt <- list(W.A.1 = W.A.1, W.A.2 = W.A.2, W.D.2 = W.D.2,
                              W.AA = W.AA, W.AD = W.AD)
              
              Zs0.null <- list(W.A.1 = Z.A.1.part, W.A.2 = Z.A.2.part, W.D.2 = Z.D.2.part)
              Zs0.alt <- list(W.A.1 = Z.A.1.part, W.A.2 = Z.A.2.part, W.D.2 = Z.D.2.part,
                              W.AA = Z.normal, W.AD = Z.normal)
              
              lin.method <- TRUE
            }else{
              K.A.1.part <- tcrossprod(W.A.1)
              K.A.2.part <- tcrossprod(W.A.2)
              K.D.2.part <- tcrossprod(W.D.2)
              
              K.A.0.1.part <- tcrossprod(W.A.0.1)
              K.A.0.2.part <- tcrossprod(W.A.0.2)
              K.D.0.2.part <- tcrossprod(W.D.0.2)
              K.AA.part <- K.A.0.1.part * K.A.0.2.part
              K.AD.part <- K.A.0.1.part * K.D.0.2.part
              K.AA.part <- K.AA.part * sqrt(nrow(K.AA.part) / sum(K.AA.part * K.AA.part))
              K.AD.part <- K.AD.part * sqrt(nrow(K.AD.part) / sum(K.AD.part * K.AD.part))
              
              ZETA.now2.null <- c(ZETA.now, list(A.1.part = list(Z = Z.A.1.part, K = K.A.1.part)),
                                  list(A.2.part = list(Z = Z.A.2.part, K = K.A.2.part)),
                                  list(D.2.part = list(Z = Z.D.2.part, K = K.D.2.part)))
              
              ZETA.now2.alt <- c(ZETA.now, list(A.1.part = list(Z = Z.A.1.part, K = K.A.1.part)),
                                 list(A.2.part = list(Z = Z.A.2.part, K = K.A.2.part)),
                                 list(D.2.part = list(Z = Z.D.2.part, K = K.D.2.part)),
                                 list(AA.part = list(Z = Z.normal, K = K.AA.part)),
                                 list(AD.part = list(Z = Z.normal, K = K.AD.part)))
              
              lin.method <- FALSE
            }
            df <- 2
          }else{
            if(isna.2){
              if(all(c(m.A.1, m.D.1) == m.A.2)){
                W.AA <- W.A.0.1 * W.A.0.2
                W.DA <- W.D.0.1 * W.A.0.2
                W.AA <- W.AA * sqrt(nrow(W.AA) / sum(W.AA * W.AA))
                W.DA <- W.DA * sqrt(nrow(W.DA) / sum(W.DA * W.DA))
                
                Ws0.null <- list(W.A.1 = W.A.1, W.A.2 = W.A.2, W.D.1 = W.D.1)
                Ws0.alt <- list(W.A.1 = W.A.1, W.A.2 = W.A.2, W.D.1 = W.D.1,
                                W.AA = W.AA, W.DA = W.DA)
                
                Zs0.null <- list(W.A.1 = Z.A.1.part, W.A.2 = Z.A.2.part, W.D.1 = Z.D.1.part)
                Zs0.alt <- list(W.A.1 = Z.A.1.part, W.A.2 = Z.A.2.part, W.D.1 = Z.D.1.part,
                                W.AA = Z.normal, W.DA = Z.normal)
                
                lin.method <- TRUE
              }else{
                K.A.1.part <- tcrossprod(W.A.1)
                K.A.2.part <- tcrossprod(W.A.2)
                K.D.1.part <- tcrossprod(W.D.1)
                
                K.A.0.1.part <- tcrossprod(W.A.0.1)
                K.A.0.2.part <- tcrossprod(W.A.0.2)
                K.D.0.1.part <- tcrossprod(W.D.0.1)
                K.AA.part <- K.A.0.1.part * K.A.0.2.part
                K.DA.part <- K.D.0.1.part * K.A.0.2.part
                K.AA.part <- K.AA.part * sqrt(nrow(K.AA.part) / sum(K.AA.part * K.AA.part))
                K.DA.part <- K.DA.part * sqrt(nrow(K.DA.part) / sum(K.DA.part * K.DA.part))
                
                ZETA.now2.null <- c(ZETA.now, list(A.1.part = list(Z = Z.A.1.part, K = K.A.1.part)),
                                    list(A.2.part = list(Z = Z.A.2.part, K = K.A.2.part)),
                                    list(D.1.part = list(Z = Z.D.1.part, K = K.D.1.part)))
                
                ZETA.now2.alt <- c(ZETA.now, list(A.1.part = list(Z = Z.A.1.part, K = K.A.1.part)),
                                   list(A.2.part = list(Z = Z.A.2.part, K = K.A.2.part)),
                                   list(D.1.part = list(Z = Z.D.1.part, K = K.D.1.part)),
                                   list(AA.part = list(Z = Z.normal, K = K.AA.part)),
                                   list(DA.part = list(Z = Z.normal, K = K.DA.part)))
                
                lin.method <- FALSE
              }
              df <- 2
            }else{
              if(all(c(m.D.1, m.A.2, m.D.2) == m.A.1)){
                W.AA <- W.A.0.1 * W.A.0.2
                W.AD <- W.A.0.1 * W.D.0.2
                W.DA <- W.D.0.1 * W.A.0.2
                W.DD <- W.D.0.1 * W.D.0.2
                W.AA <- W.AA * sqrt(nrow(W.AA) / sum(W.AA * W.AA))
                W.AD <- W.AD * sqrt(nrow(W.AD) / sum(W.AD * W.AD))
                W.DA <- W.DA * sqrt(nrow(W.DA) / sum(W.DA * W.DA))
                W.DD <- W.DD * sqrt(nrow(W.DD) / sum(W.DD * W.DD))
                
                Ws0.null <- list(W.A.1 = W.A.1, W.A.2 = W.A.2, W.D.1 = W.D.1, W.D.2 = W.D.2)
                Ws0.alt <- list(W.A.1 = W.A.1, W.A.2 = W.A.2, W.D.1 = W.D.1, W.D.2 = W.D.2,
                                W.AA = W.AA, W.AD = W.AD, W.DA = W.DA, W.DD = W.DD)
                
                Zs0.null <- list(W.A.1 = Z.A.1.part, W.A.2 = Z.A.2.part,
                                 W.D.1 = Z.D.1.part, W.D.2 = Z.D.2.part)
                Zs0.alt <- list(W.A.1 = Z.A.1.part, W.A.2 = Z.A.2.part,
                                W.D.1 = Z.D.1.part, W.D.2 = Z.D.2.part,
                                W.AA = Z.normal, W.AD = Z.normal,
                                W.DA = Z.normal, W.DD = Z.normal)
                
                lin.method <- TRUE
              }else{
                K.A.1.part <- tcrossprod(W.A.1)
                K.A.2.part <- tcrossprod(W.A.2)
                K.D.1.part <- tcrossprod(W.D.1)
                K.D.2.part <- tcrossprod(W.D.2)
                
                K.A.0.1.part <- tcrossprod(W.A.0.1)
                K.A.0.2.part <- tcrossprod(W.A.0.2)
                K.D.0.1.part <- tcrossprod(W.D.0.1)
                K.D.0.2.part <- tcrossprod(W.D.0.2)
                K.AA.part <- K.A.0.1.part * K.A.0.2.part
                K.AD.part <- K.A.0.1.part * K.D.0.2.part
                K.DA.part <- K.D.0.1.part * K.A.0.2.part
                K.DD.part <- K.D.0.1.part * K.D.0.2.part
                K.AA.part <- K.AA.part * sqrt(nrow(K.AA.part) / sum(K.AA.part * K.AA.part))
                K.AD.part <- K.AD.part * sqrt(nrow(K.AD.part) / sum(K.AD.part * K.AD.part))
                K.DA.part <- K.DA.part * sqrt(nrow(K.DA.part) / sum(K.DA.part * K.DA.part))
                K.DD.part <- K.DD.part * sqrt(nrow(K.DD.part) / sum(K.DD.part * K.DD.part))
                
                
                ZETA.now2.null <- c(ZETA.now, list(A.1.part = list(Z = Z.A.1.part, K = K.A.1.part)),
                                    list(A.2.part = list(Z = Z.A.2.part, K = K.A.2.part)),
                                    list(D.1.part = list(Z = Z.D.1.part, K = K.D.1.part)),
                                    list(D.2.part = list(Z = Z.D.2.part, K = K.D.2.part)))
                
                ZETA.now2.alt <- c(ZETA.now, list(A.1.part = list(Z = Z.A.1.part, K = K.A.1.part)),
                                   list(A.2.part = list(Z = Z.A.2.part, K = K.A.2.part)),
                                   list(D.1.part = list(Z = Z.D.1.part, K = K.D.1.part)),
                                   list(D.2.part = list(Z = Z.D.2.part, K = K.D.2.part)),
                                   list(AA.part = list(Z = Z.normal, K = K.AA.part)),
                                   list(AD.part = list(Z = Z.normal, K = K.AD.part)),
                                   list(DA.part = list(Z = Z.normal, K = K.DA.part)),
                                   list(DD.part = list(Z = Z.normal, K = K.DD.part)))
                
                
                lin.method <- FALSE
              }
              df <- 4
            }
          }
        }
      }
      
      
      if(lin.method){
        Gammas0.null <- lapply(Ws0.null, function(x) diag(ncol(x)))
        EMM.res.null <-  try(EM3.linker.cpp(y0 = y, X0 = X.now, ZETA = ZETA.now, tol = NULL, optimizer = optimizer,
                                            Zs0 = Zs0.null, Ws0 = Ws0.null, Gammas0 = Gammas0.null,
                                            gammas.diag = TRUE, X.fix = TRUE,
                                            eigen.SGS = eigen.SGS, eigen.G = eigen.G,
                                            REML = TRUE, pred = FALSE), silent = TRUE)
        
        Gammas0.alt <- lapply(Ws0.alt, function(x) diag(ncol(x)))
        EMM.res.alt <-  try(EM3.linker.cpp(y0 = y, X0 = X.now, ZETA = ZETA.now, tol = NULL, optimizer = optimizer,
                                           Zs0 = Zs0.alt, Ws0 = Ws0.alt, Gammas0 = Gammas0.alt,
                                           gammas.diag = TRUE, X.fix = TRUE,
                                           eigen.SGS = eigen.SGS, eigen.G = eigen.G,
                                           REML = TRUE, pred = FALSE), silent = TRUE)
      }else{
        EMM.res.null <- try(EM3.cpp(y = y, X0 = X.now, ZETA = ZETA.now2.null, tol = NULL, optimizer = optimizer,
                                    REML = TRUE, pred = FALSE), silent = TRUE)
        
        EMM.res.alt <- try(EM3.cpp(y = y, X0 = X.now, ZETA = ZETA.now2.alt, tol = NULL, optimizer = optimizer,
                                   REML = TRUE, pred = FALSE), silent = TRUE)
      }
      
      if(!("try-error" %in% c(class(EMM.res.null), class(EMM.res.alt)))){
        LL.null <- EMM.res.null$LL
        LL.alt <- EMM.res.alt$LL
        
        deviance <- 2 * (LL.alt - LL.null)
        score.now <- ifelse(deviance <= 0, 0, -log10((1 - chi0.mixture) *
                                                       pchisq(q = deviance, df = df, lower.tail = FALSE)))
        
        scores[i, j] <- score.now
      }
    }
  }
  scores <- scores + t(scores)
  diag(scores) <- diag(scores) / 2
  
  if(is.null(gene.set)){
    rownames(scores) <- colnames(scores) <- window.centers
  }else{
    rownames(scores) <- colnames(scores) <- gene.name
  }
  
  if (count) {
    cat("\n")
  }
  return(scores)
}







#' Calculate -log10(p) of epistatic effects with score test
#'
#'
#' @param M.now A n.sample x n.mark genotype matrix where n.sample is sample size and n.mark is the number of markers.
#' @param ZETA.now A list of variance (relationship) matrix (K; \eqn{m \times m}) and its design matrix (Z; \eqn{n \times m}) of random effects. You can use only one kernel matrix.
#' For example, ZETA = list(A = list(Z = Z, K = K))
#' Please set names of list "Z" and "K"!
#' @param y A \eqn{n \times 1} vector. A vector of phenotypic values should be used. NA is allowed.
#' @param X.now A \eqn{n \times p} matrix. You should assign mean vector (rep(1, n)) and covariates. NA is not allowed.
#' @param Gu A \eqn{n \times n} matrix. You should assign \eqn{ZKZ'}, where K is covariance (relationship) matrix and Z is its design matrix.
#' @param Ge A \eqn{n \times n} matrix. You should assign identity matrix I (diag(n)).
#' @param P0 A \eqn{n \times n} matrix. The Moore-Penrose generalized inverse of \eqn{SV0S}, where \eqn{S = X(X'X)^{-1}X'} and
#' \eqn{V0 = \sigma^2_u Gu + \sigma^2_e Ge}. \eqn{\sigma^2_u} and \eqn{\sigma^2_e} are estimators of the null model.
#' @param map Data frame of map information where the first column is the marker names,
#' the second and third column is the chromosome amd map position, and the forth column is -log10(p) for each marker.
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
#' @param chi0.mixture RAINBOWR assumes the test statistic \eqn{l1' F l1} is considered to follow a x chisq(df = 0) + (1 - a) x chisq(df = r).
#' where l1 is the first derivative of the log-likelihood and F is the Fisher information. And r is the degree of freedom.
#' The argument chi0.mixture is a (0 <= a < 1), and default is 0.5.
#' @param gene.set If you have information of gene, you can use it to perform kernel-based GWAS.
#'            You should assign your gene information to gene.set in the form of a "data.frame" (whose dimension is (the number of gene) x 2).
#'            In the first column, you should assign the gene name. And in the second column, you should assign the names of each marker,
#'            which correspond to the marker names of "geno" argument.
#' @param dominance.eff If this argument is TRUE, dominance effect is included in the model,
#' and additive x dominance and dominance x dominance are also tested as epistatic effects.
#' When you use inbred lines, please set this argument FALSE.
#' @param min.MAF Specifies the minimum minor allele frequency (MAF).
#' If a marker has a MAF less than min.MAF, it is assigned a zero score.
#' @param count When count is TRUE, you can know how far RGWAS has ended with percent display.
#'
#' @return -log10(p) of epistatic effects for each SNP-set
#'
#' @references Listgarten, J. et al. (2013) A powerful and efficient set test
#'  for genetic markers that handles confounders. Bioinformatics. 29(12): 1526-1533.
#'
#' Lippert, C. et al. (2014) Greater power and computational efficiency for kernel-based
#'  association testing of sets of genetic variants. Bioinformatics. 30(22): 3206-3214.
#'
#' Jiang, Y. and Reif, J.C. (2015) Modeling epistasis in genomic selection. Genetics. 201(2): 759-768.
#'
#'
#'
score.calc.epistasis.score <- function(M.now, y, X.now, ZETA.now, Gu, Ge, P0,
                                       map, haplotype = TRUE, num.hap = NULL, window.size.half = 5, window.slide = 1,
                                       chi0.mixture = 0.5, gene.set = NULL, dominance.eff = TRUE,
                                       min.MAF = 0.02, count = TRUE){
  chr <- map[, 2]
  chr.tab <- table(chr)
  chr.max <- max(chr)
  chr.cum <- cumsum(chr.tab)
  n.scores.each <- (chr.tab + (window.slide - 1)) %/% window.slide
  cum.n.scores <- cumsum(n.scores.each)
  if(is.null(gene.set)){
    n.scores <- sum(n.scores.each)
  }else{
    gene.names <- as.character(gene.set[, 1])
    mark.id <- as.character(gene.set[, 2])
    gene.name <- as.character(unique(gene.names))
    n.scores <- length(unique(gene.set[, 1]))
  }
  
  scores <- matrix(0, nrow = n.scores, ncol = n.scores)
  window.centers <- rep(NA, n.scores)
  probaa <- apply(M.now == -1, 2, mean)
  probAa <- apply(M.now == 0, 2, mean)
  freq <- probaa + probAa / 2
  MAF <- pmin(freq, 1 - freq)
  MAF.D <- pmin(probAa, 1 - probAa)
  
  
  n.sample.now <- nrow(M.now)
  Z.normal <- diag(n.sample.now)
  W.A.list <- W.A.0.list  <- W.D.list <- W.D.0.list <-
    Z.A.part.list <- Z.D.part.list <- rep(list(NA), n.scores)
  pb <- txtProgressBar(min = 1, max = n.scores, style = 3)
  n.scores2 <- n.scores - n.scores %% 100
  start.scorecalc <- Sys.time()
  
  for(i in 1:n.scores){
    if(is.null(gene.set)){
      i.chr <- min(which(i - cum.n.scores <= 0))
      if(i.chr >= 2){
        window.center <- window.slide * (i - cum.n.scores[i.chr - 1] - 1) + chr.cum[i.chr - 1] + 1
      }else{
        window.center <- window.slide * (i - 1) + 1
      }
      names(window.center) <- i.chr
      window.centers[i] <- window.center
      Theories1 <-  window.center < window.size.half + 1
      for(r in 1:(chr.max - 1)){
        Theory1 <- chr.cum[r] < window.center & window.center < window.size.half + 1 + chr.cum[r]
        Theories1 <- c(Theories1, Theory1)
      }
      rule1 <- sum(Theories1) != 0
      
      
      Theories2  <- NULL
      for(r in 1:chr.max){
        Theory2 <- chr.cum[r] - (window.size.half + 1) < window.center & window.center <= chr.cum[r]
        Theories2 <- c(Theories2, Theory2)
      }
      rule2 <- sum(Theories2) != 0
      
      
      if(rule1 & rule2){
        Mis.range.0 <- which(chr == i.chr)
        Mis.range.02 <- which(chr == i.chr) - window.center + 1 + window.size.half
      }else{
        if(rule1){
          near.min <- c(0, chr.cum)[which.min(abs(window.center - c(0, chr.cum)))]
          Mis.range.0 <- (near.min + 1):(window.center + window.size.half)
          Mis.range.02 <- (2 * window.size.half + 2 - length(Mis.range.0)):(2 * window.size.half + 1)
        }else{
          if(rule2){
            near.max <- c(0, chr.cum)[which.min(abs(window.center - c(0, chr.cum)))]
            Mis.range.0 <- (window.center - window.size.half):near.max
            Mis.range.02 <- 1:length(Mis.range.0)
          }else{
            Mis.range.0 <- (window.center - window.size.half):(window.center + window.size.half)
            Mis.range.02 <- 1:(2 * window.size.half + 1)
          }
        }
      }
    }else{
      mark.name.now <- mark.id[gene.names == gene.name[i]]
      Mis.range.0 <- match(mark.name.now, map[, 1])
      Mis.range.02 <- 1:length(Mis.range.0)
      weighting.center <- FALSE
    }
    
    Mis.0.0 <- M.now[, Mis.range.0, drop = FALSE]
    MAF.cut <- MAF[Mis.range.0] >= min.MAF
    if(dominance.eff){
      Mis.D.0.0 <- M.now[, Mis.range.0, drop = FALSE]
      MAF.cut.D <- MAF.D[Mis.range.0] > 0
    }else{
      MAF.cut.D <- rep(TRUE, length(MAF.cut))
    }
    
    
    if(any(MAF.cut)){
      Mis.0 <- Mis.0.0[, MAF.cut, drop = FALSE]
      Mis.range <- Mis.range.0[MAF.cut]
      Mis.range2 <- Mis.range.02[MAF.cut]
      window.size <- ncol(Mis.0)
      if(any(MAF.cut.D)){
        if(dominance.eff){
          Mis.D.0 <- Mis.D.0.0[, MAF.cut.D, drop = FALSE]
          Mis.range.D <- Mis.range.0[MAF.cut.D]
          Mis.range2.D <- Mis.range.02[MAF.cut.D]
          window.size.D <- ncol(Mis.D.0)
        }
      }
      
      if(haplotype){
        if(is.null(num.hap)){
          Mis.fac <- factor(apply(Mis.0, 1, function(x) paste(x, collapse = "")))
          Mis <- Mis.0[!duplicated(as.numeric(Mis.fac)), , drop = FALSE]
          bango <- as.factor(as.numeric(Mis.fac))
          levels(bango) <- order(unique(bango))
          bango <- as.numeric(as.character(bango))
          if(any(MAF.cut.D)){
            if(dominance.eff){
              Mis.D.fac <- factor(apply(Mis.D.0, 1, function(x) paste(x, collapse = "")))
              Mis.D <- Mis.D.0[!duplicated(as.numeric(Mis.D.fac)), , drop = FALSE]
              
              bango.D <- as.factor(as.numeric(Mis.D.fac))
              levels(bango.D) <- order(unique(bango.D))
              bango.D <- as.numeric(as.character(bango.D))
            }
          }
        }else{
          kmed.res <- cluster::pam(Mis.0, k = num.hap)
          Mis <- kmed.res$medoids
          bango <- kmed.res$clustering
          if(any(MAF.cut.D)){
            if(dominance.eff){
              kmed.res.D <- cluster::pam(Mis.D.0, k = num.hap)
              Mis.D <- kmed.res.D$medoids
              bango.D <- kmed.res.D$clustering
            }
          }
        }
        Z.part <- as.matrix(Matrix::sparseMatrix(i = 1:nrow(M.now), j = bango, x = rep(1, nrow(M.now)),
                                                 dims = c(nrow(M.now), nrow(Mis))))
        if(any(MAF.cut.D)){
          if(dominance.eff){
            Z.part.D <- as.matrix(Matrix::sparseMatrix(i = 1:nrow(M.now), j = bango.D, x = rep(1, nrow(M.now)),
                                                       dims = c(nrow(M.now), nrow(Mis.D))))
          }
        }
      }else{
        Mis <- Mis.0
        Mis.D <- Mis.D.0
        Z.part <- Z.part.D <- diag(nrow(M.now))
      }
      
      
      W.A <- calcGRM(genoMat = Mis,
                     methodGRM = "addNOIA",
                     returnWMat = TRUE,
                     probaa = probaa[Mis.range],
                     probAa = probAa[Mis.range])
      
      W.A.0 <- calcGRM(genoMat = Mis.0.0,
                       methodGRM = "addNOIA",
                       returnWMat = TRUE,
                       probaa = probaa[Mis.range.0],
                       probAa = probAa[Mis.range.0])
      
      W.A.list[[i]] <- W.A
      Z.A.part.list[[i]] <- Z.part
      W.A.0.list[[i]] <- W.A.0
      
      if(any(MAF.cut.D)){
        if(dominance.eff){
          W.D <- calcGRM(genoMat = Mis.D,
                         methodGRM = "domNOIA",
                         returnWMat = TRUE,
                         probaa = probaa[Mis.range.D],
                         probAa = probAa[Mis.range.D])
          
          W.D.0 <- calcGRM(genoMat = Mis.0.0,
                           methodGRM = "domNOIA",
                           returnWMat = TRUE,
                           probaa = probaa[Mis.range.0],
                           probAa = probAa[Mis.range.0])
          
          W.D.list[[i]] <- W.D
          Z.D.part.list[[i]] <- Z.part.D
          W.D.0.list[[i]] <- W.D.0
        }
      }
    }
  }
  
  n.calc <- n.scores * (n.scores + 1) / 2
  pb <- txtProgressBar(min = 1, max = n.calc, style = 3)
  n.calc2 <- n.calc - n.calc %% 100
  start.scorecalc <- Sys.time()
  
  
  for(i in 1:n.scores){
    W.A.1 <- W.A.list[[i]]
    Z.A.1.part <- Z.A.part.list[[i]]
    W.A.0.1 <- W.A.0.list[[i]]
    m.A.1 <- ncol(W.A.0.1)
    
    if(dominance.eff){
      W.D.1 <- W.D.list[[i]]
      Z.D.1.part <- Z.D.part.list[[i]]
      W.D.0.1 <- W.D.0.list[[i]]
      isna.1 <- any(is.na(W.D.0.1))
      m.D.1 <- ncol(W.D.0.1)
    }else{
      isna.1 <- TRUE
    }
    for(j in i:n.scores){
      num.now <- (i - 1) * n.scores + (j - i + 1)
      if(count){
        setTxtProgressBar(pb, num.now)
        if(n.calc > 100){
          if(num.now == (n.calc2 / 100 + 1) | num.now == (n.calc2 / 10 + 1) | num.now == (n.calc2 / 2 + 1)){
            cat("\n")
            end.scorecalc <- Sys.time()
            jikan.scorecalc <- (end.scorecalc - start.scorecalc) * (n.calc - num.now + 1) / (num.now - 1)
            print(paste0((num.now - 1) * 100 / n.calc2, "%...Done. ",
                         round(jikan.scorecalc, 2), " ", attr(jikan.scorecalc, "units"),
                         " to end.  Scheduled end time : ", end.scorecalc + jikan.scorecalc))
          }
        }
      }
      
      if(i == j){
        if((!dominance.eff) | isna.1){
          W.AA <- W.A.0.1 ^ 2
          W.AA <- W.AA * sqrt(nrow(W.AA) / sum(W.AA * W.AA))
          
          Ws.null <- list(W.A = Z.A.1.part %*% W.A.1)
          Ws.alt <- list(W.A = Z.A.1.part %*% W.A.1, W.AA = W.AA)
          
          gammas.diag <- TRUE
          df <- 1
        }else{
          if(m.A.1 == m.D.1){
            W.AA <- W.A.0.1 ^ 2
            W.AD <- W.A.0.1 * W.D.0.1
            W.DD <- W.D.0.1 ^ 2
            W.AA <- W.AA * sqrt(nrow(W.AA) / sum(W.AA * W.AA))
            W.AD <- W.AD * sqrt(nrow(W.AD) / sum(W.AD * W.AD))
            W.DD <- W.DD * sqrt(nrow(W.DD) / sum(W.DD * W.DD))
            
            Ws.null <- list(W.A = Z.A.1.part %*% W.A.1, W.D = Z.D.1.part %*% W.D.1)
            Ws.alt <- list(W.A = Z.A.1.part %*% W.A.1, W.D = Z.D.1.part %*% W.D.1,
                           W.AA = W.AA, W.AD = W.AD, W.DD = W.DD)
            
            gammas.diag <- TRUE
          }else{
            K.A.1.part <- tcrossprod(W.A.1)
            K.D.1.part <- tcrossprod(W.D.1)
            
            K.A.0.1.part <- tcrossprod(W.A.0.1)
            K.D.0.1.part <- tcrossprod(W.D.0.1)
            K.AA.part <- K.A.0.1.part ^ 2
            K.AD.part <- K.A.0.1.part * K.D.0.1.part
            K.DD.part <- K.D.0.1.part ^ 2
            K.AA.part <- K.AA.part * sqrt(nrow(K.AA.part) / sum(K.AA.part * K.AA.part))
            K.AD.part <- K.AD.part * sqrt(nrow(K.AD.part) / sum(K.AD.part * K.AD.part))
            K.DD.part <- K.DD.part * sqrt(nrow(K.DD.part) / sum(K.DD.part * K.DD.part))
            
            Ws.null <- list(A.part = Z.A.1.part, D.part = Z.D.1.part)
            Ws.alt <- list(A.part = Z.A.1.part, D.part = Z.D.1.part,
                           AA.part = Z.normal, AD.part = Z.normal, DD.part = Z.normal)
            
            Gammas.null <- list(A.part = K.A.1.part, D.part = K.D.1.part)
            Gammas.alt <- list(A.part = K.A.1.part, D.part = K.D.1.part,
                               AA.part = K.AA.part, AD.part = K.AD.part,
                               DD.part = K.DD.part)
            
            gammas.diag <- FALSE
          }
          df <- 3
        }
      }else{
        W.A.2 <- W.A.list[[j]]
        Z.A.2.part <- Z.A.part.list[[j]]
        W.A.0.2 <- W.A.0.list[[j]]
        m.A.2 <- ncol(W.A.0.2)
        
        if(dominance.eff){
          W.D.2 <- W.D.list[[j]]
          Z.D.2.part <- Z.D.part.list[[j]]
          W.D.0.2 <- W.D.0.list[[j]]
          isna.2 <- any(is.na(W.D.0.2))
          m.D.2 <- ncol(W.D.0.2)
          
          isnas <- c(isna.1, isna.2)
        }else{
          isna.2 <- TRUE
          isnas <- c(isna.1, isna.2)
        }
        
        if((!dominance.eff) | all(isnas)){
          if(m.A.1 == m.A.2){
            W.AA <- W.A.0.1 * W.A.0.2
            W.AA <- W.AA * sqrt(nrow(W.AA) / sum(W.AA * W.AA))
            
            Ws.null <- list(W.A.1 = Z.A.1.part %*% W.A.1, W.A.2 = Z.A.2.part %*% W.A.2)
            Ws.alt <- list(W.A.1 = Z.A.1.part %*% W.A.1, W.A.2 = Z.A.2.part %*% W.A.2, W.AA = W.AA)
            
            gammas.diag <- TRUE
          }else{
            K.A.1.part <- tcrossprod(W.A.1)
            K.A.2.part <- tcrossprod(W.A.2)
            
            K.A.0.1.part <- tcrossprod(W.A.0.1)
            K.A.0.2.part <- tcrossprod(W.A.0.2)
            K.AA.part <- K.A.0.1.part * K.A.0.2.part
            K.AA.part <- K.AA.part * sqrt(nrow(K.AA.part) / sum(K.AA.part * K.AA.part))
            
            Ws.null <- list(A.part.1 = Z.A.1.part, A.part.2 = Z.A.2.part)
            Ws.alt <- list(A.part.1 = Z.A.1.part, A.part.2 = Z.A.2.part, AA.part = Z.normal)
            
            Gammas.null <- list(A.part.1 = K.A.1.part, A.part.2 = K.A.2.part)
            Gammas.alt <- list(A.part.1 = K.A.1.part, A.part.2 = K.A.2.part, AA.part = K.AA.part)
            
            gammas.diag <- FALSE
          }
          df <- 1
        }else{
          if(isna.1){
            if(all(c(m.A.2, m.D.2) == m.A.1)){
              W.AA <- W.A.0.1 * W.A.0.2
              W.AD <- W.A.0.1 * W.D.0.2
              W.AA <- W.AA * sqrt(nrow(W.AA) / sum(W.AA * W.AA))
              W.AD <- W.AD * sqrt(nrow(W.AD) / sum(W.AD * W.AD))
              
              Ws.null <- list(W.A.1 = Z.A.1.part %*% W.A.1, W.A.2 = Z.A.2.part %*% W.A.2,
                              W.D.2 = Z.D.2.part %*% W.D.2)
              Ws.alt <- list(W.A.1 = Z.A.1.part %*% W.A.1, W.A.2 = Z.A.2.part %*% W.A.2,
                             W.D.2 = Z.D.2.part %*% W.D.2, W.AA = W.AA, W.AD = W.AD)
              
              gammas.diag <- TRUE
            }else{
              K.A.1.part <- tcrossprod(W.A.1)
              K.A.2.part <- tcrossprod(W.A.2)
              K.D.2.part <- tcrossprod(W.D.2)
              
              K.A.0.1.part <- tcrossprod(W.A.0.1)
              K.A.0.2.part <- tcrossprod(W.A.0.2)
              K.D.0.2.part <- tcrossprod(W.D.0.2)
              K.AA.part <- K.A.0.1.part * K.A.0.2.part
              K.AD.part <- K.A.0.1.part * K.D.0.2.part
              K.AA.part <- K.AA.part * sqrt(nrow(K.AA.part) / sum(K.AA.part * K.AA.part))
              K.AD.part <- K.AD.part * sqrt(nrow(K.AD.part) / sum(K.AD.part * K.AD.part))
              
              Ws.null <- list(A.part.1 = Z.A.1.part, A.part.2 = Z.A.2.part,
                              D.part.2 = Z.D.2.part)
              Ws.alt <- list(A.part.1 = Z.A.1.part, A.part.2 = Z.A.2.part,
                             D.part.2 = Z.D.2.part, AA.part = Z.normal, AD.part = Z.normal)
              
              Gammas.null <- list(A.part.1 = K.A.1.part, A.part.2 = K.A.2.part,
                                  D.part.2 = K.D.2.part)
              Gammas.alt <- list(A.part.1 = K.A.1.part, A.part.2 = K.A.2.part, D.part.2 = K.D.2.part,
                                 AA.part = K.AA.part, AD.part = K.AD.part)
              
              gammas.diag <- FALSE
            }
            df <- 2
          }else{
            if(isna.2){
              if(all(c(m.A.1, m.D.1) == m.A.2)){
                W.AA <- W.A.0.1 * W.A.0.2
                W.DA <- W.D.0.1 * W.A.0.2
                W.AA <- W.AA * sqrt(nrow(W.AA) / sum(W.AA * W.AA))
                W.AD <- W.AD * sqrt(nrow(W.AD) / sum(W.AD * W.AD))
                W.DA <- W.DA * sqrt(nrow(W.DA) / sum(W.DA * W.DA))
                
                Ws.null <- list(W.A.1 = Z.A.1.part %*% W.A.1, W.A.2 = Z.A.2.part %*% W.A.2,
                                W.D.1 = Z.D.1.part %*% W.D.1)
                Ws.alt <- list(W.A.1 = Z.A.1.part %*% W.A.1, W.A.2 = Z.A.2.part %*% W.A.2,
                               W.D.1 = Z.D.1.part %*% W.D.1, W.AA = W.AA, W.DA = W.DA)
                
                gammas.diag <- TRUE
              }else{
                K.A.1.part <- tcrossprod(W.A.1)
                K.A.2.part <- tcrossprod(W.A.2)
                K.D.1.part <- tcrossprod(W.D.1)
                
                K.A.0.1.part <- tcrossprod(W.A.0.1)
                K.A.0.2.part <- tcrossprod(W.A.0.2)
                K.D.0.1.part <- tcrossprod(W.D.0.1)
                K.AA.part <- K.A.0.1.part * K.A.0.2.part
                K.DA.part <- K.D.0.1.part * K.A.0.2.part
                K.AA.part <- K.AA.part * sqrt(nrow(K.AA.part) / sum(K.AA.part * K.AA.part))
                K.DA.part <- K.DA.part * sqrt(nrow(K.DA.part) / sum(K.DA.part * K.DA.part))
                
                Ws.null <- list(A.part.1 = Z.A.1.part, A.part.2 = Z.A.2.part,
                                D.part.1 = Z.D.1.part)
                Ws.alt <- list(A.part.1 = Z.A.1.part, A.part.2 = Z.A.2.part,
                               D.part.1 = Z.D.1.part, AA.part = Z.normal, DA.part = Z.normal)
                
                Gammas.null <- list(A.part.1 = K.A.1.part, A.part.2 = K.A.2.part,
                                    D.part.1 = K.D.1.part)
                Gammas.alt <- list(A.part.1 = K.A.1.part, A.part.2 = K.A.2.part, D.part.1 = K.D.1.part,
                                   AA.part = K.AA.part, DA.part = K.DA.part)
                
                gammas.diag <- FALSE
              }
              df <- 2
            }else{
              if(all(c(m.D.1, m.A.2, m.D.2) == m.A.1)){
                W.AA <- W.A.0.1 * W.A.0.2
                W.AD <- W.A.0.1 * W.D.0.2
                W.DA <- W.D.0.1 * W.A.0.2
                W.DD <- W.D.0.1 * W.D.0.2
                W.AA <- W.AA * sqrt(nrow(W.AA) / sum(W.AA * W.AA))
                W.AD <- W.AD * sqrt(nrow(W.AD) / sum(W.AD * W.AD))
                W.DA <- W.DA * sqrt(nrow(W.DA) / sum(W.DA * W.DA))
                W.DD <- W.DD * sqrt(nrow(W.DD) / sum(W.DD * W.DD))
                
                Ws.null <- list(W.A.1 = Z.A.1.part %*% W.A.1, W.A.2 = Z.A.2.part %*% W.A.2,
                                W.D.1 = Z.D.1.part %*% W.D.1, W.D.2 = Z.D.2.part %*% W.D.2)
                Ws.alt <- list(W.A.1 = Z.A.1.part %*% W.A.1, W.A.2 = Z.A.2.part %*% W.A.2,
                               W.D.1 = Z.D.1.part %*% W.D.1, W.D.2 = Z.D.2.part %*% W.D.2,
                               W.AA = W.AA, W.AD = W.AD, W.DA = W.DA, W.DD = W.DD)
                
                gammas.diag <- TRUE
              }else{
                K.A.1.part <- tcrossprod(W.A.1)
                K.A.2.part <- tcrossprod(W.A.2)
                K.D.1.part <- tcrossprod(W.D.1)
                K.D.2.part <- tcrossprod(W.D.2)
                
                K.A.0.1.part <- tcrossprod(W.A.0.1)
                K.A.0.2.part <- tcrossprod(W.A.0.2)
                K.D.0.1.part <- tcrossprod(W.D.0.1)
                K.D.0.2.part <- tcrossprod(W.D.0.2)
                K.AA.part <- K.A.0.1.part * K.A.0.2.part
                K.AD.part <- K.A.0.1.part * K.D.0.2.part
                K.DA.part <- K.D.0.1.part * K.A.0.2.part
                K.DD.part <- K.D.0.1.part * K.D.0.2.part
                K.AA.part <- K.AA.part * sqrt(nrow(K.AA.part) / sum(K.AA.part * K.AA.part))
                K.AD.part <- K.AD.part * sqrt(nrow(K.AD.part) / sum(K.AD.part * K.AD.part))
                K.DA.part <- K.DA.part * sqrt(nrow(K.DA.part) / sum(K.DA.part * K.DA.part))
                K.DD.part <- K.DD.part * sqrt(nrow(K.DD.part) / sum(K.DD.part * K.DD.part))
                
                Ws.null <- list(A.part.1 = Z.A.1.part, A.part.2 = Z.A.2.part,
                                D.part.1 = Z.D.1.part, D.part.2 = Z.D.2.part)
                Ws.alt <- list(A.part.1 = Z.A.1.part, A.part.2 = Z.A.2.part,
                               D.part.1 = Z.D.1.part, D.part.2 = Z.D.2.part,
                               AA.part = Z.normal, AD.part = Z.normal,
                               DA.part = Z.normal, DD.part = Z.normal)
                
                Gammas.null <- list(A.part.1 = K.A.1.part, A.part.2 = K.A.2.part,
                                    D.part.1 = K.D.1.part, D.part.2 = K.D.2.part)
                Gammas.alt <- list(A.part.1 = K.A.1.part, A.part.2 = K.A.2.part,
                                   D.part.1 = K.D.1.part, D.part.2 = K.D.2.part,
                                   AA.part = K.AA.part, AD.part = K.AD.part,
                                   DA.part = K.DA.part, DD.part = K.DD.part)
                
                
                gammas.diag <- FALSE
              }
              df <- 4
            }
          }
        }
      }
      
      
      if(gammas.diag){
        Gammas.null <- lapply(Ws.null, function(x) diag(ncol(x)))
        Gammas.alt <- lapply(Ws.alt, function(x) diag(ncol(x)))
      }
      score.null <- try(score.linker.cpp(y, Ws = Ws.null, Gammas = Gammas.null,
                                         gammas.diag = gammas.diag, Gu = Gu, Ge = Ge,
                                         P0 = P0, chi0.mixture = chi0.mixture), silent = TRUE)
      
      
      score.alt <- try(score.linker.cpp(y, Ws = Ws.alt, Gammas = Gammas.alt,
                                        gammas.diag = gammas.diag, Gu = Gu, Ge = Ge,
                                        P0 = P0, chi0.mixture = chi0.mixture), silent = TRUE)
      
      if(!("try-error" %in% c(class(score.null), class(score.alt)))){
        stat.null <- qchisq(p = 10 ^ (-score.null) / (1 - chi0.mixture), df = df, lower.tail = FALSE)
        stat.alt <- qchisq(p = 10 ^ (-score.alt) / (1 - chi0.mixture), df = df, lower.tail = FALSE)
        deviance <- stat.alt - stat.null
        score.now <- ifelse(deviance <= 0, 0, -log10((1 - chi0.mixture) *
                                                       pchisq(q = deviance, df = df, lower.tail = FALSE)))
        
        scores[i, j] <- score.now
      }
    }
  }
  scores <- scores + t(scores)
  diag(scores) <- diag(scores) / 2
  
  if(is.null(gene.set)){
    rownames(scores) <- colnames(scores) <- window.centers
  }else{
    rownames(scores) <- colnames(scores) <- gene.name
  }
  
  if (count) {
    cat("\n")
  }
  return(scores)
}

