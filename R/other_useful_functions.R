#' Function to view the first part of data (like head(), tail())
#'
#' @param data Your data. 'vector', 'matrix', 'array' (whose dimensions <= 4), 'data.frame' are supported format.
#' If other formatted data is assigned, str(data) will be returned.
#' @param fh  From head. If this argument is TRUE, first part (row) of data will be shown (like head() function).
#' If FALSE, last part (row) of your data will be shown (like tail() function).
#' @param fl From left. If this argument is TRUE, first part (column) of data will be shown (like head() function).
#' If FALSE, last part (column) of your data will be shown (like tail() function).
#' @param rown  The number of rows shown in console.
#' @param coln  The number of columns shown in console.
#' @param rowst  The start point for the direction of row.
#' @param colst  The start point for the direction of column.
#' @param narray The number of dimensions other than row and column shown in console.
#' This argument is effective only your data is array (whose dimensions >= 3).
#' @param drop  When rown = 1 or coln = 1, the dimension will be reduced if this argument is TRUE.
#' @param save.variable If you want to assign the result to a variable, please set this agument TRUE.
#' @param verbose If TRUE, print the first part of data.
#'
#' @return If save.variable is FALSE, NULL. If TRUE, the first part of your data will be returned.
#'
#'
See <- function(data, fh = TRUE, fl = TRUE, rown = 6, coln = 6,
                rowst = 1, colst = 1, narray = 2, drop = FALSE,
                save.variable = FALSE, verbose = TRUE){
  islist <- is.list(data)
  isvec <- is.vector(data)
  isfac <- is.factor(data)
  
  
  
  if ((isvec | isfac) & (!islist)) {
    n.data <- length(data)
    if (fh) {
      start <- min(rowst, n.data)
      end <- min(rowst + rown - 1, n.data)
    } else {
      start <- max(n.data - rowst - rown + 2, 1)
      end <- max(n.data - rowst + 1, 1)
    }
    data.show <- data[start:end]
    dim.show <- length(data)
  } else {
    ismat <- is.matrix(data)
    isdf <- is.data.frame(data)
    
    if (ismat | isdf) {
      n.data.row <- nrow(data)
      if (fh) {
        start.row <- min(rowst, n.data.row)
        end.row <- min(rowst + rown - 1, n.data.row)
      } else {
        start.row <- max(n.data.row - rowst - rown + 2, 1)
        end.row <- max(n.data.row - rowst + 1, 1)
      }
      
      n.data.col <- ncol(data)
      if (fl) {
        start.col <- min(colst, n.data.col)
        end.col <- min(colst + coln - 1, n.data.col)
      } else {
        start.col <- max(n.data.col - colst - coln + 2, 1)
        end.col <- max(n.data.col - colst + 1, 1)
      }
      
      class.each <- rep(NA, end.col - start.col + 1)
      
      for (i in 1:(end.col - start.col + 1)) {
        class.each[i] <- class(data[, i])
      }
      class.show <- paste0("<", class.each, ">")
      data.show <-
        data[start.row:end.row, start.col:end.col, drop = drop]
      data.show <-
        as.data.frame(rbind(class.show, apply(data.show, 2, as.character)))
      data.rowname <- rownames(data)
      if (!is.null(data.rowname)) {
        rownames(data.show) <- c("class", rownames(data)[start.row:end.row])
      } else {
        rownames(data.show) <-
          c("class", paste0("NULL_", start.row:end.row))
      }
      dim.show <- dim(data)
    } else {
      isarray <- is.array(data)
      
      if (isarray) {
        n.array <- length(dim(data))
        
        if (narray <= 4) {
          n.data.row <- nrow(data)
          if (fh) {
            start.row <- min(rowst, n.data.row)
            end.row <- min(rowst + rown - 1, n.data.row)
          } else {
            start.row <- max(n.data.row - rowst - rown + 2, 1)
            end.row <- max(n.data.row - rowst + 1, 1)
          }
          
          n.data.col <- ncol(data)
          if (fl) {
            start.col <- min(colst, n.data.col)
            end.col <- min(colst + coln - 1, n.data.col)
          } else {
            start.col <- max(n.data.col - colst - coln + 2, 1)
            end.col <- max(n.data.col - colst + 1, 1)
          }
          
          start.other <- 1
          end.other <- pmin(rep(narray, n.array - 2), dim(data)[-c(1:2)])
          
          if (n.array == 1) {
            data.show <- data[start.row:end.row, drop = drop]
          }
          
          if (n.array == 2) {
            data.show <- data[start.row:end.row, start.col:end.col, drop = drop]
          }
          
          if (n.array == 3) {
            data.show <- data[start.row:end.row, start.col:end.col,
                              start.other:end.other, drop = drop]
          }
          
          if (n.array == 4) {
            data.show <- data[start.row:end.row, start.col:end.col,
                              start.other:end.other, start.other:end.other, drop = drop]
          }
          dim.show <- dim(data)
        }else{
          stop("You can only see data whose # of the dimensions <= 4!!")
        }
      }else{
        warning("We cannot offer the simple view of your data. Instead we will offer the structure of your data.")
        data.show <- str(data)
        dim.show <-  NULL
      }
    }
  }
  
  if (verbose) {
    if (!is.null(data.show)) {
      print(data.show)
    }
    print(paste0("class: ", paste(class(data), collapse = " & ")))
    print(paste0("dimension: ", paste(dim.show, collapse = " x ")))
  }
  
  if (save.variable) {
    return(data.show)
  }
}




#' Function to remove the minor alleles
#'
#' @param x.0 A \eqn{n \times m} original marker genotype matrix.
#' @param map.0  Data frame with the marker names in the first column. The second and third columns contain the chromosome and map position.
#' @param min.MAF Specifies the minimum minor allele frequency (MAF).
#' If a marker has a MAF less than min.MAF, it is removed from the original marker genotype data.
#' @param max.MS Specifies the maximum missing rate (MS).
#' If a marker has a MS more than max.MS, it is removed from the original marker genotype data.
#' @param return.MAF If TRUE, MAF will be returned.
#'
#' @return
#' \describe{
#' \item{$x}{The modified marker genotype data whose SNPs with MAF <= min.MAF were removed.}
#' \item{$map}{The modified map information whose SNPs with MAF <= min.MAF were removed.}
#' \item{$before}{Minor allele frequencies of the original marker genotype.}
#' \item{$after}{Minor allele frequencies of the modified marker genotype.}
#'}
#'
MAF.cut <-  function(x.0, map.0 = NULL, min.MAF = 0.05,
                     max.MS = 0.05, return.MAF = FALSE) {
  x.unique <- sort(unique(c(x.0)), decreasing = FALSE)
  len.x.unique <- length(x.unique)
  
  if (len.x.unique == 2) {
    is.scoring1 <- all(x.unique == c(-1, 1))
    is.scoring2 <- all(x.unique == c(0, 2))
  } else{
    if (len.x.unique == 3) {
      is.scoring1 <- all(x.unique == c(-1, 0, 1))
      is.scoring2 <- all(x.unique == c(0, 1, 2))
    } else{
      stop("Something wrong with your genotype data!!")
    }
  }
  
  if (is.scoring1) {
    freq <- apply(x.0, 2, function(x) {
      return(mean(x + 1, na.rm = TRUE) / 2)
    })
  } else{
    if (is.scoring2) {
      freq <- apply(x.0, 2, function(x) {
        return(mean(x, na.rm = TRUE) / 2)
      })
    } else{
      stop("Genotype data should be scored with (-1, 0, 1) or (0, 1, 2)!!")
    }
  }
  
  MAF.before <- pmin(freq, 1 - freq)
  mark.remain.MAF <- MAF.before >= min.MAF
  
  MS.rate <- apply(x.0, 2, function(x)
    mean(is.na(x)))
  mark.remain.MS <- MS.rate <= max.MS
  
  mark.remain <- mark.remain.MAF & mark.remain.MS
  
  
  x <- x.0[, mark.remain]
  
  if (!is.null(map.0)) {
    map <- map.0[mark.remain,]
  } else{
    map <- NULL
  }
  MAF.after <- MAF.before[mark.remain]
  
  if (return.MAF) {
    return(list(
      data = list(x = x, map = map),
      MAF = list(before = MAF.before, after = MAF.after)
    ))
  } else{
    return(list(x = x, map = map))
  }
}











#' Function to estimate & plot phylogenetic tree
#'
#' @param blockInterest A \eqn{n \times M} matrix representing the marker genotype that belongs to the haplotype block of interest.
#' If this argument is NULL, this argument will automatically be determined by `geno`, 
#' @param gwasRes You can use the results (data.frame) of haplotype-based (SNP-set) GWAS by `RGWAS.multisnp` function. 
#' @param nTopRes Haplotype blocks (or gene sets, SNP-sets) with top `nTopRes` p-values by `gwasRes` will be used.
#' @param gene.set If you have information of gene (or haplotype block), you can use it to perform kernel-based GWAS.
#'            You should assign your gene information to gene.set in the form of a "data.frame" (whose dimension is (the number of gene) x 2).
#'            In the first column, you should assign the gene name. And in the second column, you should assign the names of each marker,
#'            which correspond to the marker names of "geno" argument.
#' @param indexRegion You can specify the haplotype block (or gene set, SNP-set) of interest by the marker index in `geno`.
#' @param chrInterest You can specify the haplotype block (or gene set, SNP-set) of interest by the marker position in `geno`.
#' Please assign the chromosome number to this argument.
#' @param posRegion You can specify the haplotype block (or gene set, SNP-set) of interest by the marker position in `geno`.
#' Please assign the position in the chromosome to this argument.
#' @param blockName You can specify the haplotype block (or gene set, SNP-set) of interest by the name of haplotype block in `geno`.
#' @param pheno Data frame where the first column is the line name (gid). 
#' The remaining columns should be a phenotype to test.
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
#' @param chi2Test If TRUE, chi-square test for the relationship between haplotypes & subpopulations will be performed. 
#' @param thresChi2Test The threshold for the chi-square test.
#' @param plotTree If TRUE, the function will return the plot of phylogenetic tree.
#' @param distMat You can assign the distance matrix of the block of interest. 
#' If NULL, the distance matrix will be computed in this function.
#' @param distMethod You can choose the method to calculate distance between accessions.
#' This argument corresponds to the `method` argument in the `dist` function.
#' @param evolutionDist If TRUE, the evolution distance will be used instead of the pure distance.
#' The `distMat` will be converted to the distance matrix by the evolution distance.
#' @param subpopInfo The information of subpopulations. This argument should be a vector of factor. 
#' @param groupingMethod If `subpopInfo` argument is NULL, this function estimates subpopulation information from marker genotype.
#' You can choose the grouping method from `kmeans`, `kmedoids`, and `hclust`. 
#' @param nGrp The number of groups (or subpopulations) grouped by `groupingMethod`.
#' If this argument is 0, the subpopulation information will not be estimated.
#' @param nIterClustering If `groupingMethod` = `kmeans`, the clustering will be performed multiple times.
#' This argument specifies the number of classification performed by the function.
#' @param kernelTypes In the function, similarlity matrix between accessions will be computed from marker genotype to estimate genotypic values.
#' This argument specifies the method to compute similarity matrix: 
#' If this argument is `addNOIA` (or one of other options in `methodGRM` in `calcGRM`), 
#' then the `addNOIA` (or corresponding) option in the `calcGRM` function will be used,
#' and if this argument is `phylo`, the gaussian kernel based on phylogenetic distance will be computed from phylogenetic tree.
#' You can assign more than one kernelTypes for this argument; for example, kernelTypes = c("addNOIA", "phylo").
#' @param nCores The number of cores used for optimization.
#' @param hOpt Optimized hyper parameter for constructing kernel when estimating haplotype effects.
#'  If hOpt = "optimized", hyper parameter will be optimized in the function.
#'  If hOpt = "tuned", hyper parameter will be replaced by the median of off-diagonal of distance matrix.
#'  If hOpt is numeric, that value will be directly used in the function.
#' @param hOpt2 Optimized hyper parameter for constructing kernel when estimating haplotype effects of nodes.
#'  If hOpt2 = "optimized", hyper parameter will be optimized in the function.
#'  If hOpt2 = "tuned", hyper parameter will be replaced by the median of off-diagonal of distance matrix.
#'  If hOpt2 is numeric, that value will be directly used in the function. 
#' @param maxIter Max number of iterations for optimization algorithm.
#' @param rangeHStart The median of off-diagonal of distance matrix multiplied by rangeHStart will be used 
#' as the initial values for optimization of hyper parameters.
#' @param saveName When drawing any plot, you can save plots in png format. In saveName, you should substitute the name you want to save.
#' When saveName = NULL, the plot is not saved.
#' @param saveStyle This argument specifies how to save the plot of phylogenetic tree.
#' The function offers `png`, `pdf`, `jpg`, and `tiff`.
#' @param pchBase A vector of two integers specifying the plot types for the positive and negative genotypic values respectively.
#' @param colNodeBase A vector of two integers or chracters specifying color of nodes for the positive and negative genotypic values respectively.
#' @param colTipBase A vector of integers or chracters specifying color of tips for the positive and negative genotypic values respectively.
#' The length of the vector should equal to the number of subpopulations.
#' @param cexMax A numeric specifying the maximum point size of the plot.
#' @param cexMin A numeric specifying the minimum point size of the plot.
#' @param edgeColoring If TRUE, the edge branch of phylogenetic tree wiil be colored.
#' @param tipLabel If TRUE, lavels for tips will be shown.
#' @param ggPlotTree If TRUE, the function will return the ggplot version of phylogenetic tree. 
#' It offers the precise information on subgroups for each haplotype.
#' @param cexMaxForGG A numeric specifying the maximum point size of the plot for ggtree,
#'  relative to the range of x and y-axes (0 < cexMaxForGG <= 1).
#' @param cexMinForGG A numeric specifying the minimum point size of the plot for ggtree,
#'  relative to the range of x and y-axes (0 < cexMaxForGG <= 1).
#' @param alphaBase alpha (parameter that indicates the opacity of a geom) for tip with positive / negative effects.
#' alpha for node will be same as the alpha for tip with negative effects.
#' @param verbose If this argument is TRUE, messages for the current step_s will be shown.
#'
#' @return
#' \describe{A list / lists of 
#' \item{$haplotypeInfo}{\describe{A list of haplotype information with 
#' \item{$haploCluster}{A vector indicating each individual belongs to which haplotypes.}
#' \item{$haploBlock}{Marker genotype of haplotype block of interest for the representing haplotypes.}
#' }
#' }
#' \item{$subpopInfo}{The information of subpopulations.}
#' \item{$distMats}{\describe{A list of distance matrix: 
#' \item{$distMat}{Distance matrix between haplotypes.}
#' \item{$distMatEvol}{Evolutionary distance matrix between haplotypes.}
#' \item{$distMatNJ}{Phylogenetic distance matrix between haplotypes including nodes.}
#' }
#' }
#' \item{$pValChi2Test}{A p-value of the chi-square test for the dependency between haplotypes & subpopulations.
#' If `chi2Test = FALSE`, `NA` will be returned.}
#' \item{$njRes}{The result of phylogenetic tree by neighborhood-joining method}
#' \item{$gvTotal}{Estimated genotypic values by kernel regression for each haplotype.}
#' \item{$gvTotalForLine}{Estimated genotypic values by kernel regression for each individual.}
#' \item{$minuslog10p}{\eqn{-log_{10}(p)} for haplotype block of interest.
#'  p is the p-value for the siginifacance of the haplotype block effect.}
#' \item{$hOpts}{Optimized hyper parameters, hOpt1 & hOpt2.}
#' \item{$EMMResults}{\describe{A list of estimated results of kernel regression: 
#' \item{$EM3Res}{Estimated results of kernel regression for the estimation of haplotype effects. (1st step)}
#' \item{$EMMRes}{Estimated results of kernel regression for the estimation of haplotype effects of nodes. (2nd step)}
#' \item{$EMM0Res}{Estimated results of kernel regression for the null model.}
#' }
#' }
#' \item{$clusterNosForHaplotype}{A list of cluster Nos of individuals that belong to each haplotype.}
#'}
#'
estPhylo <- function(blockInterest = NULL, gwasRes = NULL, nTopRes = 1, gene.set = NULL,
                     indexRegion = 1:10, chrInterest = NULL, posRegion = NULL, blockName = NULL,
                     pheno = NULL, geno = NULL, ZETA = NULL, 
                     chi2Test = TRUE, thresChi2Test = 5e-2,  plotTree = TRUE,
                     distMat = NULL, distMethod = "manhattan", evolutionDist = FALSE,
                     subpopInfo = NULL, groupingMethod = "kmedoids",
                     nGrp = 3, nIterClustering = 100, kernelTypes = "addNOIA",
                     nCores = parallel::detectCores(), hOpt = "optimized",
                     hOpt2 = "optimized", maxIter = 20, rangeHStart = 10 ^ c(-1:1),
                     saveName = NULL, saveStyle = "png",
                     pchBase = c(1, 16), colNodeBase = c(2, 4), colTipBase = c(3, 5, 6),
                     cexMax = 2, cexMin = 0.7, edgeColoring = TRUE, tipLabel = TRUE,
                     ggPlotTree = FALSE, cexMaxForGG = 0.12, cexMinForGG = 0.06,
                     alphaBase = c(0.9, 0.3), verbose = TRUE) {
  if (!is.null(geno)) {
    M <- t(geno[, -c(1:3)])
    map <- geno[, 1:3]
    
    nLine <- nrow(M)
    lineNames <- rownames(M)
    
    if (is.null(ZETA)) {
      K <- calcGRM(M)
      Z <- diag(nLine)
      
      rownames(Z) <- colnames(Z) <- lineNames
      
      ZETA <- list(A = list(Z = Z, K = K))
    }
    
    
    if (is.null(subpopInfo)) {
      if (verbose) {
        print("Now clustering...")
      }
      if (nGrp > 0) {
        if (groupingMethod == "kmeans") {
          bwSSRatios <- rep(NA, nIterClustering)
          kmResList <- NULL
          
          for (iterNo in 1:nIterClustering) {
            kmResNow <- kmeans(x = M, centers = nGrp)
            bwSSRatio <- kmResNow$betweenss / (kmResNow$betweenss + kmResNow$tot.withinss)
            
            bwSSRatios[iterNo] <- bwSSRatio
            kmResList <- c(kmResList, list(kmResNow))
          }
          
          maxNo <- which(bwSSRatios == max(bwSSRatios))
          maxNoNow <- sample(maxNo, 1)
          clusterNos <- match(kmResList[[maxNoNow]]$cluster, unique(kmResList[[maxNoNow]]$cluster))
        } else if (groupingMethod == "kmedoids") {
          kmResNow <- cluster::pam(x = M, k = nGrp)
          clusterNos <- match(kmResNow$clustering, unique(kmResNow$clustering))
        } else if (groupingMethod == "hclust") {
          tre <- hclust(dist(M, method = distMethod))
          cutreeRes <- cutree(tree = tre, k = nGrp)
          clusterNos <- match(cutreeRes, unique(cutreeRes))
        } else {
          stop("We only offer 'kmeans', 'kmedoids', and 'hclust' methods for grouping methods.")
        }
        clusterRes <- factor(paste0("cluster_", clusterNos))
        names(clusterRes) <- lineNames
        
        subpopInfo <- clusterRes
      } else {
        subpopInfo <- NULL
      }
    } else {
      if (!is.factor(subpopInfo)) {
        subpopInfo <- as.factor(subpopInfo)
      }
    }
    
    nGrp <- length(levels(subpopInfo))
    
    
    if (is.null(blockInterest)) {
      if (!is.null(gwasRes)) {
        gwasResOrd <- gwasRes[order(gwasRes[, 4], decreasing = TRUE), ]
        
        if (nTopRes == 1) {
          blockName <- as.character(gwasResOrd[1, 1])
          mrkInBlock <- gene.set[gene.set[, 1] %in% blockName, 2]
          
          blockInterest <- M[, geno[, 1] %in% mrkInBlock]
        } else {
          blockInterest <- NULL
          blockNames <- rep(NA, nTopRes)
          
          for (topNo in 1:nTopRes) {
            blockName <- as.character(gwasResOrd[topNo, 1])
            mrkInBlock <- gene.set[gene.set[, 1] %in% blockName, 2]
            
            blockInterestNow <- M[, geno[, 1] %in% mrkInBlock]
            
            blockNames[topNo] <- blockName
            blockInterest <- c(blockInterest, list(blockInterestNow))
          }
          names(blockInterest) <- blockNames
        }
      } else if (!is.null(posRegion)) {
        if (is.null(chrInterest)) {
          stop("Please input the chromosome number of interest!")
        } else {
          chrCondition <- map[, 2] %in% chrInterest
          posCondition <- (posRegion[1] <= map[, 3]) & (posRegion[2] >= map[, 3])
          
          indexRegion <- which(chrCondition & posCondition)
        }
      }
      blockInterest <- M[, indexRegion]
    } 
  } else {
    blockInterestCheck <- !is.null(blockInterest)
    ZETACheck <- !is.null(ZETA)
    
    if (blockInterestCheck & ZETACheck) {
      lineNames <- rownames(blockInterest)
      nLine <- nrow(blockInterest)
    } else {
      stop("Please input 'geno', or both 'blockInterest' and 'ZETA'.")
    }
  }
  
  estPhyloEachBlock <- function(blockInterest,
                                blockName = NULL) {
    
    stringBlock <- apply(blockInterest, 1, function(x) paste0(x, collapse = ""))
    blockInterestUnique <- blockInterest[!duplicated(stringBlock), ]
    nHaplo <- length(unique(stringBlock))
    haploClusterNow <- as.numeric(factor(stringBlock))
    names(haploClusterNow) <- lineNames
    
    blockInterestUniqueSorted <- blockInterestUnique[order(unique(stringBlock)), ]
    haploNames <- paste0("haplo_", 1:nHaplo)
    rownames(blockInterestUniqueSorted) <- haploNames
    haploCluster <- haploNames[haploClusterNow]
    names(haploCluster) <- lineNames
    
    if (!is.null(subpopInfo)) {
      tableRes <- table(haploCluster = haploCluster,
                        subpop = subpopInfo)[haploNames, ]
      
      if (verbose) {
        cat("Tabular for haplotype x subpopulation: \n")
        print(head(tableRes))
        cat("\nThe number of haplotypes: ", nHaplo, "\n")
      }
      
      if (chi2Test) {
        chi2TestRes <- chisq.test(tableRes)
        pValChi2Test <- chi2TestRes$p.value
        
        if (verbose) {
          cat("\n")
          print(chi2TestRes)
          
          cat("Haplotypes & Subpopulations are: \n")
          if (pValChi2Test <= thresChi2Test) {
            cat("\tSiginificantly dependent\n")
          } else {
            cat("\tindependent\n")
          }
        }
      } else {
        pValChi2Test <- NA
      }
    }
    
    haplotypeInfo <- list(haploCluster = haploCluster,
                          haploBlock = blockInterestUniqueSorted)
    
    
    nMrkInBlock <- ncol(blockInterest)
    
    if (is.null(distMat)) {
      distMat <- as.matrix(dist(blockInterestUniqueSorted, method = distMethod))
    }
    
    if (evolutionDist) {
      inLog <- 1 - 4 * (distMat / (ncol(blockInterest) * 2)) / 3
      
      inLog[inLog <= 0] <- min(inLog[inLog > 0]) / 10
      
      dist4Nj <- - 3 * log(inLog) / 4
    } else {
      dist4Nj <- distMat
    }
    
    njRes <- ape::nj(X = as.dist(dist4Nj))
    distNodes <- ape::dist.nodes(njRes) / nMrkInBlock
    
    minuslog10ps <- c() 
    gvEstTotals <- gvEstTotalForLines <- 
      hOptsList <- EMMResultsList <- list()
    hOptBase <- hOpt
    hOptBase2 <- hOpt2
    
    
    for (kernelType in kernelTypes){
      if (verbose) {
        print(paste0("Now optimizing for kernelType: ", kernelType))
      }
      
      
      if (!is.null(pheno)) {
        rownames(distNodes)[1:nHaplo] <- colnames(distNodes)[1:nHaplo] <- haploNames
        
        nTotal <- nrow(distNodes)
        nNode <- nTotal - nHaplo
        
        ZgKernelPart <- as.matrix(Matrix::sparseMatrix(i = 1:nLine,
                                                       j = haploClusterNow,
                                                       x = rep(1, nLine),
                                                       dims = c(nLine, nHaplo),
                                                       dimnames = list(lineNames, haploNames)))
        
        hInv <- median(distNodes[upper.tri(distNodes)])
        h <- 1 / hInv
        hStarts <- h * rangeHStart
        hStarts <- split(hStarts, factor(1:length(rangeHStart)))
        
        if (kernelType %in% c("phylo", "gaussian", "exponential")){
          if (hOptBase == "optimized") {
            if (verbose) {
              print("Now optimizing hyperparameter for estimating haplotype effects...")
            }
            
            maximizeFunc <- function(h) {
              if (kernelType == "phylo") {
                gKernel <- exp(- h * distNodes)
                
                gKernelPart <- gKernel[1:nHaplo, 1:nHaplo]
              } else {
                gKernelPart <- calcGRM(blockInterestUniqueSorted,
                                       methodGRM = kernelType,
                                       kernel.h = h)
              }
              
              ZETANow <- c(ZETA, list(Part = list(Z = ZgKernelPart, K = gKernelPart)))
              EM3Res <- EM3.cpp(y = pheno[, 2], ZETA = ZETANow)
              LL <- EM3Res$LL
              
              return(-LL)
            }
            
            
            if (length(hStarts) >= 2) {
              if (verbose) {
                solnList <- pbmcapply::pbmclapply(X = hStarts,
                                                  FUN = function(h) {
                                                    soln <- nlminb(start = h, objective = maximizeFunc, gradient = NULL, hessian = NULL,
                                                                   lower = 0, upper = 1e06, control = list(iter.max = maxIter))
                                                    
                                                    return(soln)
                                                  }, mc.cores = nCores)
              } else {
                solnList <- pbapply::pblapply(X = hStarts,
                                              FUN = function(h) {
                                                soln <- nlminb(start = h, objective = maximizeFunc, gradient = NULL, hessian = NULL,
                                                               lower = 0, upper = 1e06, control = list(iter.max = maxIter))
                                                
                                                return(soln)
                                              }, mc.cores = nCores)
              }
              solnNo <- which.min(unlist(lapply(solnList, function(x) x$objective)))
              soln <- solnList[[solnNo]]
            } else {
              traceInside <- ifelse(verbose, 1, 0)
              soln <- nlminb(start = hStarts[[1]], objective = maximizeFunc, gradient = NULL, hessian = NULL,
                             lower = 0, upper = 1e06, control = list(trace = traceInside, iter.max = maxIter))
            }
            
            hOpt <- soln$par
          } else if (hOptBase == "tuned") {
            hOpt <- h
          } else if (!is.numeric(hOptBase)) {
            stop("`hOpt` should be either one of 'optimized', 'tuned', or numeric!!")
          }
          
          
          if (kernelType == "phylo") {
            gKernel <- exp(- hOpt * distNodes)
            
            gKernelPart <- gKernel[1:nHaplo, 1:nHaplo]
          } else {
            gKernelPart <- calcGRM(blockInterestUniqueSorted,
                                   methodGRM = kernelType,
                                   kernel.h = hOpt)
          }
        } else {
          hOpt <- NA
          gKernelPart <- calcGRM(blockInterestUniqueSorted,
                                 methodGRM = kernelType)
        }
        
        
        if (verbose) {
          print("Now estimating genotypic values...")
        }
        
        ZETANow <- c(ZETA, list(Part = list(Z = ZgKernelPart, K = gKernelPart)))
        EM3Res <- EM3.cpp(y = pheno[, 2], ZETA = ZETANow)
        LL <- EM3Res$LL
        gvEst <- EM3Res$u[(nLine + 1):(nLine + nHaplo), ]
        EMMRes0 <- EMM.cpp(y = pheno[, 2], ZETA = ZETA)
        LL0 <- EMMRes0$LL
        
        pVal <- pchisq(2 * (LL - LL0), df = 1, lower.tail = FALSE)
        minuslog10p <- - log10(pVal)
        if (EM3Res$weights[2] <= 1e-06) {
          warning("This block seems to have no effect on phenotype...")
          plotNode <- FALSE
          
          gvEstTotal <- c(gvEst, rep(NA, nNode))
          names(gvEstTotal) <- rownames(distNodes)
          
          cexTip <- cexMax / 2
          pchTip <- pchBase[2]
          EMMRes <- NA
          hOpt2 <- NA
          
          cexTipForGG <- rep(cexMaxForGG / sqrt(2), nHaplo)
        } else {
          plotNode <- TRUE
          
          ZgKernel <- diag(nTotal)
          rownames(ZgKernel) <- colnames(ZgKernel) <- rownames(distNodes)
          
          gvEst2 <- matrix(c(gvEst, rep(NA, nNode)))
          rownames(gvEst2) <- rownames(distNodes)
          
          if (hOptBase2 == "optimized") {   
            if (verbose) {
              print("Now optimizing hyperparameter for estimating haplotype effects of nodes...")
            }
            
            maximizeFunc2 <- function(h) {
              gKernel <- exp(- h * distNodes)
              ZETA2 <- list(Part = list(Z = ZgKernel, K = gKernel))
              
              EMMRes <- EMM.cpp(y = gvEst2, ZETA = ZETA2)
              LL <- EMMRes$LL
              
              return(-LL)
            }
            
            
            if (length(hStarts) >= 2) {
              if (verbose) {
                solnList2 <- pbmcapply::pbmclapply(X = hStarts,
                                                   FUN = function(h) {
                                                     soln <- nlminb(start = h, objective = maximizeFunc2, gradient = NULL, hessian = NULL,
                                                                    lower = 0, upper = 1e06, control = list(iter.max = maxIter))
                                                     
                                                     return(soln)
                                                   }, mc.cores = nCores)
              } else {
                solnList2 <- pbapply::pblapply(X = hStarts,
                                               FUN = function(h) {
                                                 soln <- nlminb(start = h, objective = maximizeFunc2, gradient = NULL, hessian = NULL,
                                                                lower = 0, upper = 1e06, control = list(iter.max = maxIter))
                                                 
                                                 return(soln)
                                               }, mc.cores = nCores)
              }
              solnNo2 <- which.min(unlist(lapply(solnList2, function(x) x$objective)))
              soln2 <- solnList2[[solnNo2]]
            } else {
              traceInside <- ifelse(verbose, 1, 0)
              soln2 <- nlminb(start = hStarts[[1]], objective = maximizeFunc2, gradient = NULL, hessian = NULL,
                              lower = 0, upper = 1e06, control = list(trace = traceInside, iter.max = maxIter))
            }
            
            hOpt2 <- soln2$par
          } else if (hOptBase2 == "tuned") {
            hOpt2 <- h
          } else if (!is.numeric(hOptBase2)) {
            stop("`hOpt2` should be either one of 'optimized', 'tuned', or numeric!!")
          }
          
          gKernel2 <- exp(- hOpt2 * distNodes)
          ZETA2 <- list(Part = list(Z = ZgKernel, K = gKernel2))
          
          EMMRes <- EMM.cpp(y = gvEst2, ZETA = ZETA2)
          
          gvNode <- EMMRes$u[(nHaplo + 1):nTotal] + EMMRes$beta
          gvEstTotal <- c(gvEst, gvNode)
          names(gvEstTotal) <- rownames(distNodes)
          gvCentered <- gvEstTotal - mean(gvEstTotal)
          gvScaled <- gvCentered / sd(gvCentered)
          gvScaled4Cex <- gvScaled * (cexMax - cexMin) / max(abs(gvScaled)) 
          
          cexNode <- (abs(gvScaled4Cex) + cexMin)[(nHaplo + 1):nTotal]
          pchNode <- ifelse(gvScaled[(nHaplo + 1):nTotal] > 0, pchBase[1], pchBase[2])
          colNode <- ifelse(gvScaled[(nHaplo + 1):nTotal] > 0, colNodeBase[1], colNodeBase[2])
          
          cexTip <- (abs(gvScaled4Cex) + cexMin)[1:nHaplo]
          pchTip <- ifelse(gvScaled[1:nHaplo] > 0, pchBase[1], pchBase[2])
          
          
          colorForNodeFac <- factor(gvScaled[(nHaplo + 1):nTotal] > 0,
                                    levels = c("TRUE", "FALSE"))
          colorForNodeDF <- data.frame(model.matrix(~ colorForNodeFac - 1))
          rownames(colorForNodeDF) <- names(colorForNodeFac)
          colnames(colorForNodeDF) <- c("gvNode +", "gvNode -")
          
          colorForNodeDF$node <- rownames(colorForNodeDF)
          
          
          gvScaled4CexForGG <- sign(gvScaled) * sqrt(abs(gvScaled)) * 
            (cexMaxForGG - cexMinForGG) / max(sqrt(abs(gvScaled)))
          
          cexNodeForGG <- (abs(gvScaled4CexForGG) + cexMinForGG)[(nHaplo + 1):nTotal]
          cexTipForGG <- (abs(gvScaled4CexForGG) + cexMinForGG)[1:nHaplo]
          
          alphaTip <- ifelse(gvScaled[1:nHaplo] > 0, alphaBase[1], alphaBase[2])
          alphaNode <- alphaBase[2]
        }
      } else {
        plotNode <- FALSE
        
        nNode <- njRes$Nnode
        nTotal <- nNode + nHaplo
        gvEstTotal <- rep(NA, nTotal)
        names(gvEstTotal) <- rownames(distNodes)
        minuslog10p <- NA
        
        cexTip <- cexMax / 2
        pchTip <- pchBase[2]
        EM3Res <- EMMRes <- EMMRes0 <- NA
        hOpt <- hOpt2 <- NA
        
        cexTipForGG <- rep(cexMaxForGG / sqrt(2), nHaplo)
      }
      gvEstTotalForLine <- (gvEstTotal[haploNames])[haploClusterNow]
      names(gvEstTotalForLine) <- lineNames
      
      
      hOpts <- c(hOpt = hOpt,
                 hOpt2 = hOpt2)
      EMMResults <- list(EM3Res = EM3Res,
                         EMMRes = EMMRes,
                         EMMRes0 = EMMRes0)
      
      
      
      gvEstTotals[[kernelType]] <- gvEstTotal 
      gvEstTotalForLines[[kernelType]] <- gvEstTotalForLine 
      minuslog10ps[kernelType] <- minuslog10p 
      hOptsList[[kernelType]] <- hOpts 
      EMMResultsList[[kernelType]] <- EMMResults
      
      if (!is.null(subpopInfo)) {
        if (length(colTipBase) != nGrp) {
          stop("The length of 'colTipBase' should be equal to 'nGrp' or the number of subpopulations!")
        }
        
        colTipNo <- as.numeric(subpopInfo)
        colTip <- colTipBase[colTipNo]
        names(colTipNo) <- names(colTip) <- lineNames
        
        colTip <- tapply(colTip, INDEX = haploCluster, FUN = function(x) {
          as.numeric(names(which.max(table(x) / sum(table(x)))))
        })[haploNames]
        
        clusterNosForHaplotype <- tapply(subpopInfo, INDEX = haploCluster,
                                         FUN = table)[haploNames]
        
        edgeCol <- as.numeric(colTip[njRes$edge[, 2]])
        
        clusterNosDF <- data.frame(do.call(what = rbind, 
                                           args = clusterNosForHaplotype))
        clusterNosRatioDF <- data.frame(t(apply(clusterNosDF, 1, 
                                                function(x) x / sum(x))))
        
        clusterNosRatioDF$node <- 1:nHaplo
      } else {
        colTip <- rep("gray", nHaplo)
        names(colTip) <- haploNames
        colTipBase <- "gray"
        
        edgeCol <- colTip[njRes$edge[, 2]]
        clusterNosForHaplotype <- NA
        
        clusterNosRatioDF <- data.frame(matrix(data = 1,
                                               nrow = nHaplo,
                                               ncol = 1))
        rownames(clusterNosRatioDF) <- haploNames
        colnames(clusterNosRatioDF) <- "Tip"
        clusterNosRatioDF$node <- 1:nHaplo
      }
      
      edgeCol[is.na(edgeCol)] <- "gray"
      
      if (plotTree) {
        if (!is.null(saveName)) {
          savePlotNameBase <- paste0(paste(c(saveName, blockName, kernelType), collapse = "_"),
                                     "_phylogenetic_tree")
          if (saveStyle == "pdf") {
            savePlotName <- paste0(savePlotNameBase, ".pdf")
            pdf(file = savePlotName, width = 12, height = 9)
          } else if (saveStyle == "jpg") {
            savePlotName <- paste0(savePlotNameBase, ".jpg")  
            jpeg(filename = savePlotName, width = 800, height = 600)
          } else if (saveStyle == "tiff") {
            savePlotName <- paste0(savePlotNameBase, ".tiff")
            tiff(filename = savePlotName, width = 800, height = 600)
          } else {
            savePlotName <- paste0(savePlotNameBase, ".png")
            png(filename = savePlotName, width = 800, height = 600)
          }
        }
        
        if (edgeColoring) {
          ape::plot.phylo(njRes, type = "u", show.tip.label = F, edge.color = edgeCol)
        } else {
          ape::plot.phylo(njRes, type = "u", show.tip.label = F, edge.color = "gray30")
        }
        if (plotNode) {
          ape::nodelabels(pch = pchNode, cex = cexNode, col = colNode)
        }
        if (tipLabel) {
          ape::tiplabels(pch = pchTip, cex = cexTip, col = colTip)
        }
        title(main = paste0(paste(c(colnames(pheno)[2], 
                                    blockName, kernelType), collapse = "_"),
                            " (-log10p: ", round(minuslog10p, 2), ")"))
        if (plotNode) {
          if (!is.null(subpopInfo)) {
            legend("topleft", legend = c(paste0(rep(levels(subpopInfo), each = 2),
                                                rep(c(" (gv:+)",  " (gv:-)"), nGrp)),
                                         "Node (gv:+)", "Node (gv:-)"),
                   col = c(rep(colTipBase, each = 2), colNodeBase),
                   pch = rep(pchBase, nGrp + 1))
          } else {
            legend("topleft", legend = c("Tip (gv:+)", "Tip (gv:-)",
                                         "Node (gv:+)", "Node (gv:-)"),
                   col = c(rep(colTipBase, each = 2), colNodeBase),
                   pch = rep(pchBase, 2))
          }
        } else {
          if (!is.null(subpopInfo)) {
            legend("topleft", legend = levels(subpopInfo),
                   col = colTipBase,
                   pch = pchTip)
          } else {
            legend("topleft", legend = "Tip",
                   col = colTipBase,
                   pch = pchTip)
          }   
        }
        if (!is.null(saveName)) {
          dev.off()
        } else {
          Sys.sleep(0.5)
        }
      }
      
      
      
      
      
      if (ggPlotTree) {
        trPhylo41 <- as(njRes, 'phylo4')
        
        if (edgeColoring) {
          dataForTipCol <- data.frame(color = colTip)
        } else {
          dataForTipCol <- data.frame(color = rep("gray30", nHaplo))
        }
        rownames(dataForTipCol) <- njRes$tip.label
        trPhylo42 <- phylobase::phylo4d(trPhylo41, dataForTipCol)
        
        nodeData <- data.frame(randomTrait = gvEstTotal[(nHaplo + 1):nTotal],
                               color = rep("gray", phylobase::nNodes(trPhylo41)),
                               row.names = phylobase::nodeId(trPhylo41, "internal"))
        phylobase::nodeData(trPhylo42) <- nodeData
        
        
        plt <- ggtree::ggtree(tr = trPhylo42,
                              ggtree::aes(col = I(.data$color)),
                              layout = "equal_angle")
        
        
        if (plotNode) {
          piesNode <- ggtree::nodepie(colorForNodeDF, cols = 1:2,
                                      color = colNodeBase,
                                      alpha = alphaNode)
          
          plt <- plt + ggtree::geom_inset(insets = piesNode,
                                          width = cexNodeForGG,
                                          height = cexNodeForGG)
        }
        
        if (tipLabel) {
          if (!is.null(subpopInfo)) {
            colTipBaseForPie <- colTipBase
          } else {
            colTipBaseForPie <- rep("gray", nGrp)
          }
          
          if (!all(is.na(EMMRes))) {
            piesPlus <- ggtree::nodepie(clusterNosRatioDF[alphaTip == alphaBase[1], ],
                                        cols = 1:nGrp,
                                        color = colTipBaseForPie,
                                        alpha = alphaBase[1])
            piesMinus <- ggtree::nodepie(clusterNosRatioDF[alphaTip == alphaBase[2], ],
                                         cols = 1:nGrp,
                                         color = colTipBaseForPie,
                                         alpha = alphaBase[2])
            
            plt <- plt + ggtree::geom_inset(insets = piesPlus,
                                            width = cexTipForGG[alphaTip == alphaBase[1]],
                                            height = cexTipForGG[alphaTip == alphaBase[1]])
            
            plt <- plt + ggtree::geom_inset(insets = piesMinus,
                                            width = cexTipForGG[alphaTip == alphaBase[2]],
                                            height = cexTipForGG[alphaTip == alphaBase[2]])
          } else {
            piesTip <- ggtree::nodepie(clusterNosRatioDF,
                                       cols = 1:nGrp,
                                       color = colTipBaseForPie,
                                       alpha = mean(alphaBase)) 
            
            plt <- plt + ggtree::geom_inset(insets = piesTip,
                                            width = cexTipForGG,
                                            height = cexTipForGG)
          }
        }
        plt <- plt + ggplot2::ggtitle(label = paste0(paste(c(colnames(pheno)[2], 
                                                             blockName, kernelType), collapse = "_"),
                                                     " (-log10p: ", round(minuslog10p, 2), ")"))
        
        
        if (!is.null(saveName)) {
          savePlotNameBase <- paste0(paste(c(saveName, blockName, kernelType), collapse = "_"),
                                     "_phylogenetic_ggtree")
          if (saveStyle == "pdf") {
            savePlotName <- paste0(savePlotNameBase, ".pdf")
            pdf(file = savePlotName, width = 12, height = 9)
          } else if (saveStyle == "jpg") {
            savePlotName <- paste0(savePlotNameBase, ".jpg")  
            jpeg(filename = savePlotName, width = 800, height = 600)
          } else if (saveStyle == "tiff") {
            savePlotName <- paste0(savePlotNameBase, ".tiff")
            tiff(filename = savePlotName, width = 800, height = 600)
          } else {
            savePlotName <- paste0(savePlotNameBase, ".png")
            png(filename = savePlotName, width = 800, height = 600)
          }
        }
        
        print(plt)
        
        if (!is.null(saveName)) {
          dev.off()
        } else {
          Sys.sleep(0.5)
        }
      }
      
    }
    
    
    
    return(list(haplotypeInfo = haplotypeInfo,
                subpopInfo = subpopInfo,
                pValChi2Test = pValChi2Test,
                distMats = list(distMat = distMat,
                                distMatEvol = dist4Nj,
                                distMatNJ = distNodes),
                njRes = njRes,
                gvTotal = gvEstTotals,
                gvTotalForLine = gvEstTotalForLines,
                minuslog10p = minuslog10ps,
                hOpts = hOptsList,
                EMMResults = EMMResultsList,
                clusterNosForHaplotype = clusterNosForHaplotype))
  }
  
  if (is.matrix(blockInterest)) {
    return(estPhyloEachBlock(blockInterest = blockInterest,
                             blockName = blockName))
  } else if (is.list(blockInterest)) {
    nTopRes <- length(blockInterest)
    allResultsList <- rep(list(NULL), nTopRes)
    blockNames <- names(blockInterest)
    
    if (is.null(blockNames)) {
      blockNames <- paste0("block_", 1:nTopRes)
    }
    names(allResultsList) <- blockNames
    
    
    for (topNo in 1:nTopRes) {
      allResultsList[[topNo]] <- estPhyloEachBlock(blockInterest = blockInterest[[topNo]],
                                                   blockName = blockNames[topNo])
    }
    
    return(allResultsList)
  }
}















#' Function to estimate & plot haplotype network
#'
#' @param blockInterest A \eqn{n \times M} matrix representing the marker genotype that belongs to the haplotype block of interest.
#' If this argument is NULL, this argument will automatically be determined by `geno`, 
#' @param gwasRes You can use the results (data.frame) of haplotype-based (SNP-set) GWAS by `RGWAS.multisnp` function. 
#' @param nTopRes Haplotype blocks (or gene sets, SNP-sets) with top `nTopRes` p-values by `gwasRes` will be used.
#' @param gene.set If you have information of gene (or haplotype block), you can use it to perform kernel-based GWAS.
#'            You should assign your gene information to gene.set in the form of a "data.frame" (whose dimension is (the number of gene) x 2).
#'            In the first column, you should assign the gene name. And in the second column, you should assign the names of each marker,
#'            which correspond to the marker names of "geno" argument.
#' @param indexRegion You can specify the haplotype block (or gene set, SNP-set) of interest by the marker index in `geno`.
#' @param chrInterest You can specify the haplotype block (or gene set, SNP-set) of interest by the marker position in `geno`.
#' Please assign the chromosome number to this argument.
#' @param posRegion You can specify the haplotype block (or gene set, SNP-set) of interest by the marker position in `geno`.
#' Please assign the position in the chromosome to this argument.
#' @param blockName You can specify the haplotype block (or gene set, SNP-set) of interest by the name of haplotype block in `geno`.
#' @param pheno Data frame where the first column is the line name (gid). 
#' The remaining columns should be a phenotype to test.
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
#' @param chi2Test If TRUE, chi-square test for the relationship between haplotypes & subpopulations will be performed. 
#' @param thresChi2Test The threshold for the chi-square test.
#' @param plotNetwork If TRUE, the function will return the plot of haplotype network.
#' @param distMat You can assign the distance matrix of the block of interest. 
#' If NULL, the distance matrix will be computed in this function.
#' @param evolutionDist If TRUE, the evolution distance will be used instead of the pure distance.
#' The `distMat` will be converted to the distance matrix by the evolution distance when you use `complementHaplo = "phylo"`.
#' @param distMethod You can choose the method to calculate distance between accessions.
#' This argument corresponds to the `method` argument in the `dist` function.
#' @param complementHaplo how to complement unobserved haplotypes.
#' When `complementHaplo = "all"`, all possible haplotypes will be complemented from the observed haplotypes.
#' When `complementHaplo = "never"`, unobserved haplotypes will not be complemented.
#' When `complementHaplo = "phylo"`, unobserved haplotypes will be complemented as nodes of phylogenetic tree.
#' When `complementHaplo = "TCS"`, unobserved haplotypes will be complemented by TCS methods (Clement et al., 2002).
#' @param subpopInfo The information of subpopulations. This argument should be a vector of factor. 
#' @param groupingMethod If `subpopInfo` argument is NULL, this function estimates subpopulation information from marker genotype.
#' You can choose the grouping method from `kmeans`, `kmedoids`, and `hclust`. 
#' @param nGrp The number of groups (or subpopulations) grouped by `groupingMethod`.
#' If this argument is 0, the subpopulation information will not be estimated.
#' @param nIterClustering If `groupingMethod` = `kmeans`, the clustering will be performed multiple times.
#' This argument specifies the number of classification performed by the function.
#' @param iterRmst The number of iterations for RMST (randomized minimum spanning tree).
#' @param networkMethod Either one of 'mst' (minimum spanning tree),
#'  'msn' (minimum spanning network), and 'rmst' (randomized minimum spanning tree).
#'  'rmst' is recommended.
#' @param autogamous This argument will be valid only when you use `complementHaplo = "all"` or `complementHaplo = "TCS"`.
#' This argument specifies whether the plant is autogamous or not. If autogamous = TRUE, 
#' complemented haplotype will consist of only homozygous sites ({-1, 1}). 
#' If FALSE, complemented haplotype will consist of both homozygous & heterozygous sites ({-1, 0, 1}).
#' @param probParsimony Equal to the argument `prob` in `haplotypes::parsimnet` function: 
#' 
#' A numeric vector of length one in the range [0.01, 0.99] giving the probability of parsimony as defined in Templeton et al. (1992). 
#' In order to set maximum connection steps to Inf (to connect all the haplotypes in a single network), set the probability to NULL.
#' @param nMaxHaplo The maximum number of haplotypes. If the number of total (complemented + original) haplotypes are larger than `nMaxHaplo`, 
#' we will only show the results only for the original haplotypes to reduce the computational time.
#' @param kernelTypes In the function, similarlity matrix between accessions will be computed from marker genotype to estimate genotypic values.
#' This argument specifies the method to compute similarity matrix: 
#' If this argument is `addNOIA` (or one of other options in `methodGRM` in `calcGRM`), 
#' then the `addNOIA` (or corresponding) option in the `calcGRM` function will be used,
#' and if this argument is `diffusion`, the diffusion kernel based on Laplacian matrix will be computed from network.
#' You can assign more than one kernelTypes for this argument; for example, kernelTypes = c("addNOIA", "diffusion").
#' @param nCores The number of cores used for optimization.
#' @param hOpt Optimized hyper parameter for constructing kernel when estimating haplotype effects.
#'  If hOpt = "optimized", hyper parameter will be optimized in the function.
#'  If hOpt is numeric, that value will be directly used in the function.
#' @param hOpt2 Optimized hyper parameter for constructing kernel when estimating complemented haplotype effects.
#'  If hOpt2 = "optimized", hyper parameter will be optimized in the function.
#'  If hOpt2 is numeric, that value will be directly used in the function. 
#' @param maxIter Max number of iterations for optimization algorithm.
#' @param rangeHStart The median of off-diagonal of distance matrix multiplied by rangeHStart will be used 
#' as the initial values for optimization of hyper parameters.
#' @param saveName When drawing any plot, you can save plots in png format. In saveName, you should substitute the name you want to save.
#' When saveName = NULL, the plot is not saved.
#' @param saveStyle This argument specifies how to save the plot of phylogenetic tree.
#' The function offers `png`, `pdf`, `jpg`, and `tiff`.
#' @param plotWhichMDS We will show the MDS (multi-dimensional scaling) plot, 
#' and this argument is a vector of two integers specifying that will define which MDS dimension will be plotted.
#' The first and second integers correspond to the horizontal and vertical axes, respectively.
#' @param colConnection A vector of two integers or characters specifying the colors of connection between nodes for the original and complemented haplotypes, respectively.
#' @param ltyConnection A vector of two characters specifying the line types of connection between nodes for the original and complemented haplotypes, respectively.
#' @param lwdConnection A vector of two integers specifying the line widths of connection between nodes for the original and complemented haplotypes, respectively.
#' @param pchBase A vector of two integers specifying the plot types for the positive and negative genotypic values respectively.
#' @param colCompBase A vector of two integers or characters specifying color of complemented haplotypes for the positive and negative genotypic values respectively.
#' @param colHaploBase A vector of integers or characters specifying color of original haplotypes for the positive and negative genotypic values respectively.
#' The length of the vector should equal to the number of subpopulations.
#' @param cexMax A numeric specifying the maximum point size of the plot.
#' @param cexMin A numeric specifying the minimum point size of the plot.
#' @param ggPlotNetwork If TRUE, the function will return the ggplot version of haplotype network. 
#' It offers the precise information on subgroups for each haplotype.
#' @param cexMaxForGG A numeric specifying the maximum point size of the plot for the ggplot version of haplotype network,
#'  relative to the range of x and y-axes (0 < cexMaxForGG <= 1).
#' @param cexMinForGG A numeric specifying the minimum point size of the plot for the ggplot version of haplotype network,
#'  relative to the range of x and y-axes (0 < cexMaxForGG <= 1).
#' @param alphaBase alpha (parameter that indicates the opacity of a geom) for original haplotype with positive / negative effects.
#' alpha for complemented haplotype will be same as the alpha for original haplotype with negative effects.
#' @param verbose If this argument is TRUE, messages for the current steps will be shown.
#'
#' @return
#' \describe{A list / lists of 
#' \item{$haplotypeInfo}{\describe{A list of haplotype information with 
#' \item{$haploCluster}{A vector indicating each individual belongs to which haplotypes.}
#' \item{$haploBlock}{Marker genotype of haplotype block of interest for the representing haplotypes.}
#' }
#' }
#' \item{$subpopInfo}{The information of subpopulations.}
#' \item{$pValChi2Test}{A p-value of the chi-square test for the dependency between haplotypes & subpopulations.
#' If `chi2Test = FALSE`, `NA` will be returned.}
#' \item{$mstResults}{\describe{A list of estimated results of MST / MSN / RMST: 
#' \item{$mstRes}{Estimated results of MST / MSN / RMST for the data including original haplotypes.}
#' \item{$mstResComp}{Estimated results of MST / MSN / RMST for the data including both original and complemented haplotype.}
#' }
#' }
#' \item{$distMats}{\describe{A list of distance matrix: 
#' \item{$distMat}{Distance matrix between haplotypes.}
#' \item{$distMatComp}{Distance matrix between haplotypes (including unobserved ones).}
#' \item{$laplacianMat}{Laplacian matrix between haplotypes (including unobserved ones).}
#' }
#' }
#' \item{$gvTotal}{Estimated genotypic values by kernel regression for each haplotype.}
#' \item{$gvTotalForLine}{Estimated genotypic values by kernel regression for each individual.}
#' \item{$minuslog10p}{\eqn{-log_{10}(p)} for haplotype block of interest.
#'  p is the p-value for the siginifacance of the haplotype block effect.}
#' \item{$hOpts}{Optimized hyper parameters, hOpt1 & hOpt2.}
#' \item{$EMMResults}{\describe{A list of estimated results of kernel regression: 
#' \item{$EM3Res}{Estimated results of kernel regression for the estimation of haplotype effects. (1st step)}
#' \item{$EMMRes}{Estimated results of kernel regression for the estimation of haplotype effects of nodes. (2nd step)}
#' \item{$EMM0Res}{Estimated results of kernel regression for the null model.}
#' }
#' }
#' \item{$clusterNosForHaplotype}{A list of cluster Nos of individuals that belong to each haplotype.}
#'}
#'
estNetwork <- function(blockInterest = NULL, gwasRes = NULL, nTopRes = 1, gene.set = NULL,
                       indexRegion = 1:10, chrInterest = NULL, posRegion = NULL, blockName = NULL,
                       pheno = NULL, geno = NULL, ZETA = NULL, chi2Test = TRUE, 
                       thresChi2Test = 5e-2,  plotNetwork = TRUE, distMat = NULL,
                       distMethod = "manhattan", evolutionDist = FALSE, complementHaplo = "phylo",
                       subpopInfo = NULL, groupingMethod = "kmedoids", nGrp = 3, 
                       nIterClustering = 100, iterRmst = 100, networkMethod = "rmst",
                       autogamous = FALSE, probParsimony = 0.95, nMaxHaplo = 1000,
                       kernelTypes = "addNOIA", nCores = parallel::detectCores(),
                       hOpt = "optimized", hOpt2 = "optimized", maxIter = 20,
                       rangeHStart = 10 ^ c(-1:1), saveName = NULL, saveStyle = "png",
                       plotWhichMDS = 1:2, colConnection = c("grey40", "grey60"),
                       ltyConnection = c("solid", "dashed"), lwdConnection = c(1.5, 0.8),
                       pchBase = c(1, 16), colCompBase = c(2, 4),
                       colHaploBase = c(3, 5, 6), cexMax = 2, cexMin = 0.7,
                       ggPlotNetwork = FALSE, cexMaxForGG = 0.025, cexMinForGG = 0.008,
                       alphaBase = c(0.9, 0.3), verbose = TRUE) {
  
  if (!is.null(geno)) {
    M <- t(geno[, -c(1:3)])
    map <- geno[, 1:3]
    
    nLine <- nrow(M)
    lineNames <- rownames(M)
    
    if (is.null(ZETA)) {
      K <- calcGRM(M)
      Z <- diag(nLine)
      
      rownames(Z) <- colnames(Z) <- lineNames
      
      ZETA <- list(A = list(Z = Z, K = K))
    }
    
    
    if (is.null(subpopInfo)) {
      if (verbose) {
        print("Now clustering...")
      }
      if (nGrp > 0) {
        if (groupingMethod == "kmeans") {
          bwSSRatios <- rep(NA, nIterClustering)
          kmResList <- NULL
          
          for (iterNo in 1:nIterClustering) {
            kmResNow <- kmeans(x = M, centers = nGrp)
            bwSSRatio <- kmResNow$betweenss / (kmResNow$betweenss + kmResNow$tot.withinss)
            
            bwSSRatios[iterNo] <- bwSSRatio
            kmResList <- c(kmResList, list(kmResNow))
          }
          
          maxNo <- which(bwSSRatios == max(bwSSRatios))
          maxNoNow <- sample(maxNo, 1)
          clusterNos <- match(kmResList[[maxNoNow]]$cluster, unique(kmResList[[maxNoNow]]$cluster))
        } else if (groupingMethod == "kmedoids") {
          kmResNow <- cluster::pam(x = M, k = nGrp)
          clusterNos <- match(kmResNow$clustering, unique(kmResNow$clustering))
        } else if (groupingMethod == "hclust") {
          tre <- hclust(dist(M, method = distMethod))
          cutreeRes <- cutree(tree = tre, k = nGrp)
          clusterNos <- match(cutreeRes, unique(cutreeRes))
        } else {
          stop("We only offer 'kmeans', 'kmedoids', and 'hclust' methods for grouping methods.")
        }
        clusterRes <- factor(paste0("cluster_", clusterNos))
        names(clusterRes) <- lineNames
        
        subpopInfo <- clusterRes
      } else {
        subpopInfo <- NULL
      }
    } else {
      if (!is.factor(subpopInfo)) {
        subpopInfo <- as.factor(subpopInfo)
      }
    }
    
    nGrp <- length(levels(subpopInfo))
    
    
    if (is.null(blockInterest)) {
      if (!is.null(gwasRes)) {
        gwasResOrd <- gwasRes[order(gwasRes[, 4], decreasing = TRUE), ]
        
        if (nTopRes == 1) {
          blockName <- as.character(gwasResOrd[1, 1])
          mrkInBlock <- gene.set[gene.set[, 1] %in% blockName, 2]
          
          blockInterest <- M[, geno[, 1] %in% mrkInBlock]
        } else {
          blockInterest <- NULL
          blockNames <- rep(NA, nTopRes)
          
          for (topNo in 1:nTopRes) {
            blockName <- as.character(gwasResOrd[topNo, 1])
            mrkInBlock <- gene.set[gene.set[, 1] %in% blockName, 2]
            
            blockInterestNow <- M[, geno[, 1] %in% mrkInBlock]
            
            blockNames[topNo] <- blockName
            blockInterest <- c(blockInterest, list(blockInterestNow))
          }
          names(blockInterest) <- blockNames
        }
      } else if (!is.null(posRegion)) {
        if (is.null(chrInterest)) {
          stop("Please input the chromosome number of interest!")
        } else {
          chrCondition <- map[, 2] %in% chrInterest
          posCondition <- (posRegion[1] <= map[, 3]) & (posRegion[2] >= map[, 3])
          
          indexRegion <- which(chrCondition & posCondition)
        }
      }
      blockInterest <- M[, indexRegion]
    } 
  } else {
    blockInterestCheck <- !is.null(blockInterest)
    ZETACheck <- !is.null(ZETA)
    
    if (blockInterestCheck & ZETACheck) {
      lineNames <- rownames(blockInterest)
      nLine <- nrow(blockInterest)
    } else {
      stop("Please input 'geno', or both 'blockInterest' and 'ZETA'.")
    }
  }
  
  estNetworkEachBlock <- function(blockInterest,
                                  blockName = NULL) {
    stringBlock <- apply(blockInterest, 1, function(x) paste0(x, collapse = ""))
    blockInterestUnique <- blockInterest[!duplicated(stringBlock), ]
    nHaplo <- length(unique(stringBlock))
    haploClusterNow <- as.numeric(factor(stringBlock))
    lineNames <- rownames(blockInterest)
    names(haploClusterNow) <- lineNames
    
    
    blockInterestUniqueSorted <- blockInterestUnique[order(unique(stringBlock)), ]
    haploNames <- paste0("h", 1:nHaplo)
    rownames(blockInterestUniqueSorted) <- haploNames
    stringBlockUniqueSorted <- apply(blockInterestUniqueSorted, 1,
                                     function(x) paste0(x, collapse = ""))
    haploCluster <- haploNames[haploClusterNow]
    names(haploCluster) <- lineNames
    
    if (!is.null(subpopInfo)) {
      tableRes <- table(haploCluster = haploCluster,
                        subpop = subpopInfo)[haploNames, ]
      
      if (verbose) {
        cat("Tabular for haplotype x subpopulation: \n")
        print(head(tableRes))
        cat("\nThe number of haplotypes: ", nHaplo, "\n")
      }
      
      if (chi2Test) {
        chi2TestRes <- chisq.test(tableRes)
        pValChi2Test <- chi2TestRes$p.value
        
        if (verbose) {
          cat("\n")
          print(chi2TestRes)
          
          cat("Haplotypes & Subpopulations are: \n")
          if (pValChi2Test <= thresChi2Test) {
            cat("\tSiginificantly dependent\n")
          } else {
            cat("\tindependent\n")
          }
        }
      } else {
        pValChi2Test <- NA
      }
    }
    
    haplotypeInfo <- list(haploCluster = haploCluster,
                          haploBlock = blockInterestUniqueSorted)
    
    
    nMrkInBlock <- ncol(blockInterest)
    
    if (is.null(distMat)) {
      distMat <- as.matrix(dist(blockInterestUniqueSorted, method = distMethod))
    }
    
    if (networkMethod == "rmst") {
      mstRes <- pegas::rmst(d = distMat, B = iterRmst)
      
      mstResPlus <- attr(mstRes, "alter.links")
      mstResAll <- rbind(mstRes,
                         mstResPlus)
    } else if (networkMethod == "mst") {
      mstRes <- pegas::mst(d = distMat)
      mstResPlus <- NULL
      mstResAll <- mstRes
    } else if (networkMethod == "msn") {
      mstRes <- pegas::msn(d = distMat)
      mstResPlus <- attr(mstRes, "alter.links")
      mstResAll <- rbind(mstRes,
                         mstResPlus)
    } else {
      stop("We only offer 'rmst', 'mst', and 'msn' for `networkMethod`!!")
    }
    if (complementHaplo %in% c("all", "never")) {
      if (complementHaplo == "all") {
        blockInterestComp <- do.call(
          what = rbind,
          args = sapply(1:nrow(mstResAll), function(eachComb) {
            
            blocksNow <- blockInterestUniqueSorted[mstResAll[eachComb, 1:2], ]
            diffBlocks <- diff(blocksNow)
            whichDiff <- which(diffBlocks != 0)
            nDiff <- length(whichDiff)
            
            blocksPart <- blocksNow[, whichDiff, drop = FALSE]
            gridList <- lapply(apply(blocksPart, 2, function(eachMrk) {
              gridCand <- min(eachMrk):max(eachMrk)
              
              if (autogamous) {
                gridCand <- gridCand[gridCand != 0]
              } 
              
              return(list(gridCand))
            }), function(x) x[[1]])
            newBlocksPart <- as.matrix(expand.grid(gridList))
            nNewBlocks <- nrow(newBlocksPart)
            
            newBlocks <- matrix(data = rep(blocksNow[1, ], nNewBlocks),
                                nrow = nNewBlocks,
                                ncol = ncol(blocksNow),
                                byrow = TRUE)
            colnames(newBlocks) <- colnames(blocksNow)
            newBlocks[, whichDiff] <- newBlocksPart
            
            return(newBlocks)
          }, simplify = FALSE)
        )
        
        if (nrow(blockInterestComp) > nMaxHaplo) {
          warning("There are too many complemented haplotypes... We will show the results only for the original haplotypes.")
          blockInterestComp <- blockInterestUniqueSorted
        }
      } else if (complementHaplo == "never") {
        blockInterestComp <- blockInterestUniqueSorted
      } 
      
      blockInterestComp <- blockInterestComp[!duplicated(blockInterestComp), ]
      stringBlockComp <- apply(blockInterestComp, 1, 
                               function(x) paste0(x, collapse = ""))
      nHaploComp <- nrow(blockInterestComp)
      
      blockInterestCompSorted <- blockInterestComp[order(stringBlockComp), ]
      stringBlockCompSorted <- apply(blockInterestCompSorted, 1, 
                                     function(x) paste0(x, collapse = ""))
      
      matchString <- match(stringBlockCompSorted, stringBlockUniqueSorted)
      
      namesBlockInterestComp <- rownames(blockInterestUniqueSorted)[matchString]
      nComp <- sum(is.na(namesBlockInterestComp))
      if (nComp >= 1){
        existComp <- TRUE
        compNames <- paste0("c", 1:nComp)
      } else {
        existComp <- FALSE
        compNames <- NULL
      }
      namesBlockInterestComp[is.na(namesBlockInterestComp)] <- compNames
      
      rownames(blockInterestCompSorted) <- namesBlockInterestComp
      
      haplotypeInfo$haploBlockCompSorted <- blockInterestCompSorted
      
      distMatComp <- dist(blockInterestCompSorted, method = distMethod)
    } else if (complementHaplo == "phylo") {
      haplotypeInfo$haploBlockCompSorted <- NULL
      
      if (evolutionDist) {
        inLog <- 1 - 4 * (distMat / (ncol(blockInterest) * 2)) / 3
        
        inLog[inLog <= 0] <- min(inLog[inLog > 0]) / 10
        
        dist4Nj <- - 3 * log(inLog) / 4
      } else {
        dist4Nj <- distMat
      }
      
      njRes <- ape::nj(X = as.dist(dist4Nj))
      
      distMatComp <- ape::dist.nodes(njRes)
      nHaploComp <- nrow(distMatComp)
      nComp <- nHaploComp - nHaplo
      
      if (nComp >= 1){
        existComp <- TRUE
        compNames <- paste0("c", 1:nComp)
      } else {
        existComp <- FALSE
        compNames <- NULL
      }
      namesBlockInterestComp <- c(haploNames, compNames)
      rownames(distMatComp) <- colnames(distMatComp) <- namesBlockInterestComp
    } else if (complementHaplo == "TCS") {
      extractParsimnetRes <- function (parsimnetRes) {
        nHapEachNet <- parsimnetRes@nhap
        distForEachNet <- parsimnetRes@d
        nNet <- length(distForEachNet)
        
        nTotalEachNet <- unlist(lapply(distForEachNet, nrow))
        nCompEachNet <- nTotalEachNet - nHapEachNet
        
        haploNamesByTCS <- unlist(sapply(1:nNet, function(netNo) {
          distNow <- distForEachNet[[netNo]]
          haploNamesNow <- rownames(distNow)
          nHapNow <- nHapEachNet[netNo]
          nTotalNow <- nTotalEachNet[netNo]
          
          if (nCompEachNet[netNo] >= 1) {
            
            for (compNo in (nHapNow + 1):nTotalNow) {
              compNameNow <- haploNamesNow[compNo]
              x <- stringr::str_split(string = compNameNow, pattern = "_")[[1]]
              
              haploNos <- as.numeric(x[1:2])
              if (is.na(x[3])) {
                compNameNew <- paste(haploNamesNow[haploNos],
                                     collapse = "_")
              } else {
                compNameNew <- paste(c(haploNamesNow[haploNos],
                                       x[3]), collapse = "_")
              }
              haploNamesNow[compNo] <- compNameNew 
            }
          }
          
          return(haploNamesNow)
        }, simplify = FALSE))
        
        
        return(list(nCompEachNet = nCompEachNet, 
                    distForEachNet = distForEachNet,
                    haploNamesByTCS = haploNamesByTCS))
      }
      
      haplotypeInfo$haploBlockCompSorted <- NULL
      
      if (autogamous) {
        parsimnetRes <- haplotypes::parsimnet(x = distMat / 2,
                                              seqlength = nMrkInBlock,
                                              prob = probParsimony)
        parsimnetResAll <- haplotypes::parsimnet(x = distMat / 2,
                                                 seqlength = nMrkInBlock,
                                                 prob = NULL)
      } else {
        parsimnetRes <- haplotypes::parsimnet(x = distMat,
                                              seqlength = nMrkInBlock)
        parsimnetResAll <- haplotypes::parsimnet(x = distMat,
                                                 seqlength = nMrkInBlock,
                                                 prob = NULL)
        
      }
      parsimnetResults <- list(parsimnetRes = parsimnetRes,
                               parsimnetResForOneNet = parsimnetResAll)
      haplotypeInfo$parsimnetResults <- parsimnetResults
      
      extractInfoByTCS <- extractParsimnetRes(parsimnetRes = parsimnetRes)
      extractInfoByTCSAll <- extractParsimnetRes(parsimnetRes = parsimnetResAll)
      haploNamesByTCS <- extractInfoByTCS$haploNamesByTCS
      distMatByTCSAll <- extractInfoByTCSAll$distForEachNet$net1
      distMatComp <- distMatByTCSAll[haploNamesByTCS, haploNamesByTCS]
      
      if (autogamous) {
        distMatComp <- distMatComp * 2
      }
      
      nComp <- sum(extractInfoByTCS$nCompEachNet)
      nHaploComp <- nHaplo + nComp
      
      
      if (nComp >= 1){
        existComp <- TRUE
        compNames <- paste0("c", 1:nComp)
      } else {
        existComp <- FALSE
        compNames <- NULL
      }
      namesBlockInterestComp <- haploNamesByTCS
      namesBlockInterestComp[!(namesBlockInterestComp %in% haploNames)] <- compNames
      rownames(distMatComp) <- colnames(distMatComp) <- namesBlockInterestComp
    }
    
    if (networkMethod == "rmst") {
      mstResComp <- pegas::rmst(d = distMatComp, B = iterRmst)
      
      mstResCompPlus <- attr(mstResComp, "alter.links")
      mstResCompAll <- rbind(mstResComp,
                             mstResCompPlus)
    } else if (networkMethod == "mst") {
      mstResComp <- pegas::mst(d = distMatComp)
      mstResCompPlus <- NULL
      mstResCompAll <- mstResComp
    } else if (networkMethod == "msn") {
      mstResComp <- pegas::msn(d = distMatComp)
      mstResCompPlus <- attr(mstResComp, "alter.links")
      mstResCompAll <- rbind(mstResComp,
                             mstResCompPlus)
    } else {
      stop("We only offer 'rmst', 'mst', and 'msn' for `networkMethod`!!")
    }
    
    L <- as.matrix(Matrix::sparseMatrix(i = mstResCompAll[, 1],
                                        j = mstResCompAll[, 2],
                                        x = rep(-1, nrow(mstResCompAll)),
                                        dims = c(nHaploComp, nHaploComp),
                                        dimnames = list(namesBlockInterestComp,
                                                        namesBlockInterestComp)))
    L <- L + t(L)
    diag(L) <- - apply(L, 2, sum)
    
    
    nTotal <- nrow(L)
    plotPlus <- !is.null(mstResCompPlus)
    
    
    minuslog10ps <- c() 
    gvEstTotals <- gvEstTotalForLines <- 
      hOptsList <- EMMResultsList <- list()
    
    hOptBase <- hOpt
    hOptBase2 <- hOpt2
    
    for (kernelType in kernelTypes){
      if (verbose) {
        print(paste0("Now optimizing for kernelType: ", kernelType))
      }
      
      
      if (!is.null(pheno)) {
        ZgKernelPart <- as.matrix(Matrix::sparseMatrix(i = 1:nLine,
                                                       j = haploClusterNow,
                                                       x = rep(1, nLine),
                                                       dims = c(nLine, nHaplo),
                                                       dimnames = list(lineNames, haploNames)))
        
        h <- 1
        hStarts <- h * rangeHStart
        hStarts <- split(hStarts, factor(1:length(rangeHStart)))
        if (kernelType %in% c("diffusion", "gaussian", "exponential")) {
          if (hOptBase == "optimized") {
            if (verbose) {
              print("Now optimizing hyperparameter for estimating haplotype effects...")
            }
            
            
            maximizeFunc <- function(h) {
              if (kernelType == "diffusion") {
                gKernel <- expm::expm(- h * L)
                
                gKernelPart <- gKernel[haploNames, haploNames]
              } else {
                gKernelPart <- calcGRM(blockInterestUniqueSorted,
                                       methodGRM = kernelType,
                                       kernel.h = h)
              }
              
              ZETANow <- c(ZETA, list(Part = list(Z = ZgKernelPart, K = gKernelPart)))
              EM3Res <- EM3.cpp(y = pheno[, 2], ZETA = ZETANow)
              LL <- EM3Res$LL
              
              return(-LL)
            }
            
            
            if (length(hStarts) >= 2) {
              if (verbose) {
                solnList <- pbmcapply::pbmclapply(X = hStarts,
                                                  FUN = function(h) {
                                                    soln <- nlminb(start = h, objective = maximizeFunc, gradient = NULL, hessian = NULL,
                                                                   lower = 0, upper = 1e06, control = list(iter.max = maxIter))
                                                    
                                                    return(soln)
                                                  }, mc.cores = nCores)
              } else {
                solnList <- pbapply::pblapply(X = hStarts,
                                              FUN = function(h) {
                                                soln <- nlminb(start = h, objective = maximizeFunc, gradient = NULL, hessian = NULL,
                                                               lower = 0, upper = 1e06, control = list(iter.max = maxIter))
                                                
                                                return(soln)
                                              }, mc.cores = nCores)
              }
              solnNo <- which.min(unlist(lapply(solnList, function(x) x$objective)))
              soln <- solnList[[solnNo]]
            } else {
              traceInside <- ifelse(verbose, 1, 0)
              soln <- nlminb(start = hStarts[[1]], objective = maximizeFunc, gradient = NULL, hessian = NULL,
                             lower = 0, upper = 1e06, control = list(trace = traceInside, iter.max = maxIter))
            }
            
            hOpt <- soln$par
          } else if (!is.numeric(hOptBase)) {
            stop("`hOpt` should be either one of 'optimized' or numeric!!")
          }
          
          
          if (kernelType == "diffusion") {
            gKernel <- expm::expm(- hOpt * L)
            
            gKernelPart <- gKernel[haploNames, haploNames]
          } else {
            gKernelPart <- calcGRM(blockInterestUniqueSorted,
                                   methodGRM = kernelType,
                                   kernel.h = hOpt)
          }
        } else {
          hOpt <- NA
          gKernelPart <- calcGRM(blockInterestUniqueSorted,
                                 methodGRM = kernelType)
        }
        
        
        if (verbose) {
          print("Now estimating genotypic values...")
        }
        
        ZETANow <- c(ZETA, list(Part = list(Z = ZgKernelPart, K = gKernelPart)))
        EM3Res <- EM3.cpp(y = pheno[, 2], ZETA = ZETANow)
        gvEst <- EM3Res$u[(nLine + 1):(nLine + nHaplo), ]
        LL <- EM3Res$LL
        EMMRes0 <- EMM.cpp(y = pheno[, 2], ZETA = ZETA)
        LL0 <- EMMRes0$LL
        
        pVal <- pchisq(2 * (LL - LL0), df = 1, lower.tail = FALSE)
        minuslog10p <- - log10(pVal)
        
        
        if ((EM3Res$weights[2] <= 1e-06)) {
          warning("This block seems to have no effect on phenotype...")
          gvEstTotal <- rep(NA, nTotal)
          names(gvEstTotal) <- namesBlockInterestComp
          gvEstTotal[haploNames] <- gvEst
          
          cexHaplo <- cexMax / 2
          cexComp <- cexMax / 4
          pchHaplo <- pchComp <- pchBase[2]
          colComp <- colCompBase[2]
          EMMRes <- NA
          hOpt2 <- NA
          
          
          cexHaploForGG <- cexMaxForGG / sqrt(2)
          cexCompForGG <- cexMaxForGG / 2
          cexAll <- rep(NA, nTotal)
          names(cexAll) <- namesBlockInterestComp
          cexAll[haploNames] <- cexHaploForGG
          if (existComp) {
            cexAll[compNames] <- cexCompForGG
          }
        } else {
          if (existComp){
            ZgKernel <- diag(nTotal)
            rownames(ZgKernel) <- colnames(ZgKernel) <- namesBlockInterestComp
            
            gvEst2 <- matrix(rep(NA, nTotal))
            rownames(gvEst2) <- namesBlockInterestComp
            gvEst2[haploNames, ] <- gvEst
            if (hOptBase2 == "optimized") {  
              
              if (verbose) {
                print("Now optimizing hyperparameter for estimating complemented haplotype effects...")
              }
              
              maximizeFunc2 <- function(h) {
                gKernel <- expm::expm(- h * L)
                ZETA2 <- list(Part = list(Z = ZgKernel, K = gKernel))
                
                EMMRes <- EMM.cpp(y = gvEst2, ZETA = ZETA2)
                LL <- EMMRes$LL
                
                return(-LL)
              }
              
              
              if (length(hStarts) >= 2) {
                if (verbose) {
                  solnList2 <- pbmcapply::pbmclapply(X = hStarts,
                                                     FUN = function(h) {
                                                       soln <- nlminb(start = h, objective = maximizeFunc2, gradient = NULL, hessian = NULL,
                                                                      lower = 0, upper = 1e06, control = list(iter.max = maxIter))
                                                       
                                                       return(soln)
                                                     }, mc.cores = nCores)
                } else {
                  solnList2 <- pbapply::pblapply(X = hStarts,
                                                 FUN = function(h) {
                                                   soln <- nlminb(start = h, objective = maximizeFunc2, gradient = NULL, hessian = NULL,
                                                                  lower = 0, upper = 1e06, control = list(iter.max = maxIter))
                                                   
                                                   return(soln)
                                                 }, mc.cores = nCores)
                }
                solnNo2 <- which.min(unlist(lapply(solnList2, function(x) x$objective)))
                soln2 <- solnList2[[solnNo2]]
              } else {
                traceInside <- ifelse(verbose, 1, 0)
                soln2 <- nlminb(start = hStarts[[1]], objective = maximizeFunc2, gradient = NULL, hessian = NULL,
                                lower = 0, upper = 1e06, control = list(trace = traceInside, iter.max = maxIter))
              }
              
              hOpt2 <- soln2$par
            } else if (!is.numeric(hOptBase2)) {
              stop("`hOpt2` should be either one of 'optimized' or numeric!!")
            }
            
            gKernel2 <- expm::expm(- hOpt2 * L)
            ZETA2 <- list(Part = list(Z = ZgKernel, K = gKernel2))
            
            EMMRes <- EMM.cpp(y = gvEst2, ZETA = ZETA2)
            
            u <- EMMRes$u
            names(u) <- namesBlockInterestComp
            gvPlus <- u[compNames] + EMMRes$beta
            gvEstTotal <- gvEst2[, 1]
            gvEstTotal[compNames] <- gvPlus
            gvCentered <- gvEstTotal - mean(gvEstTotal)
            gvScaled <- gvCentered / sd(gvCentered)
            gvScaled4Cex <- gvScaled * (cexMax - cexMin) / max(abs(gvScaled)) 
            
            cexComp <- (abs(gvScaled4Cex) + cexMin)[compNames]
            pchComp <- ifelse(gvScaled[compNames] > 0, pchBase[1], pchBase[2])
            colComp <- ifelse(gvScaled[compNames] > 0, colCompBase[1], colCompBase[2])
            
            cexHaplo <- (abs(gvScaled4Cex) + cexMin)[haploNames]
            pchHaplo <- ifelse(gvScaled[haploNames] > 0, pchBase[1], pchBase[2])
            
            
            
            
            gvScaled4CexForGG <- gvScaled * (cexMaxForGG - cexMinForGG) / max(abs(gvScaled))
            cexAll <- abs(gvScaled4CexForGG) + cexMinForGG
            cexHaploForGG <- cexAll[haploNames]
            cexCompForGG <- cexAll[compNames]
          } else {
            gvEstTotal <- gvEst
            names(gvEstTotal) <- haploNames
            
            gvCentered <- gvEstTotal - mean(gvEstTotal)
            gvScaled <- gvCentered / sd(gvCentered)
            gvScaled4Cex <- gvScaled * (cexMax - cexMin) / max(abs(gvScaled)) 
            
            cexHaplo <- (abs(gvScaled4Cex) + cexMin)[haploNames]
            pchHaplo <- ifelse(gvScaled[haploNames] > 0, pchBase[1], pchBase[2])
            
            cexComp <- cexMax / 4
            colComp <- colCompBase[2]
            pchComp <- pchBase[2]
            EMMRes <- NA
            hOpt2 <- NA
            
            
            gvScaled4CexForGG <- gvScaled * (cexMaxForGG - cexMinForGG) / max(abs(gvScaled))
            cexAll <- abs(gvScaled4CexForGG) + cexMinForGG
            cexHaploForGG <- cexAll[haploNames]
            cexCompForGG <- cexMaxForGG / 2
          }
        }
      } else {
        gvEstTotal <- rep(NA, nTotal)
        names(gvEstTotal) <- namesBlockInterestComp 
        minuslog10p <- NA
        
        cexHaplo <- cexMax / 2
        cexComp <- cexMax / 4
        pchHaplo <- pchComp <- pchBase[2]
        colComp <- colCompBase[2]
        EM3Res <- EMMRes <- EMMRes0 <- NA
        hOpt <- hOpt2 <- NA
        
        
        cexHaploForGG <- cexMaxForGG / sqrt(2)
        cexCompForGG <- cexMaxForGG / 2
        cexAll <- rep(NA, nTotal)
        names(cexAll) <- namesBlockInterestComp
        cexAll[haploNames] <- cexHaploForGG
        if (existComp) {
          cexAll[compNames] <- cexCompForGG
        }
        
      }
      gvEstTotalForLine <- (gvEstTotal[haploNames])[haploClusterNow]
      names(gvEstTotalForLine) <- lineNames
      
      hOpts <- c(hOpt = hOpt,
                 hOpt2 = hOpt2)
      EMMResults <- list(EM3Res = EM3Res,
                         EMMRes = EMMRes,
                         EMMRes0 = EMMRes0)
      
      
      
      gvEstTotals[[kernelType]] <- gvEstTotal 
      gvEstTotalForLines[[kernelType]] <- gvEstTotalForLine 
      minuslog10ps[kernelType] <- minuslog10p 
      hOptsList[[kernelType]] <- hOpts 
      EMMResultsList[[kernelType]] <- EMMResults
      
      
      if (!is.null(subpopInfo)) {
        if (length(colHaploBase) != nGrp) {
          stop("The length of 'colHaploBase' should be equal to 'nGrp' or the number of subpopulations!")
        }
        
        colHaploNo <- as.numeric(subpopInfo)
        colHaplo <- colHaploBase[colHaploNo]
        names(colHaploNo) <- names(colHaplo) <- lineNames
        
        colHaplo <- tapply(colHaplo, INDEX = haploCluster, FUN = function(x) {
          as.numeric(names(which.max(table(x) / sum(table(x)))))
        })[haploNames]
        
        clusterNosForHaplotype <- tapply(subpopInfo, INDEX = haploCluster,
                                         FUN = table)[haploNames]
      } else {
        colHaplo <- rep("gray", nHaplo)
        names(colHaplo) <- haploNames
        
        clusterNosForHaplotype <- NA
      }
      
      
      nDimMDS <- max(plotWhichMDS)
      mdsResComp <- cmdscale(d = distMatComp,
                             k = nDimMDS,
                             eig = TRUE)
      rangeX <- range(mdsResComp$points[, plotWhichMDS[1]])
      rangeY <- range(mdsResComp$points[, plotWhichMDS[2]])
      
      
      if (plotNetwork) {
        if (!is.null(saveName)) {
          savePlotNameBase <- paste0(paste(c(saveName, blockName, kernelType), collapse = "_"),
                                     "_network")
          if (saveStyle == "pdf") {
            savePlotName <- paste0(savePlotNameBase, ".pdf")
            pdf(file = savePlotName, width = 12, height = 9)
          } else if (saveStyle == "jpg") {
            savePlotName <- paste0(savePlotNameBase, ".jpg")  
            jpeg(filename = savePlotName, width = 800, height = 600)
          } else if (saveStyle == "tiff") {
            savePlotName <- paste0(savePlotNameBase, ".tiff")
            tiff(filename = savePlotName, width = 800, height = 600)
          } else {
            savePlotName <- paste0(savePlotNameBase, ".png")
            png(filename = savePlotName, width = 800, height = 600)
          }
        }
        
        
        rangeX[1] <- rangeX[1] - abs(diff(rangeX)) / 4
        rangeY[2] <- rangeY[2] + abs(diff(rangeY)) / 3
        
        plot(x = mdsResComp$points[haploNames, plotWhichMDS[1]],
             y = mdsResComp$points[haploNames, plotWhichMDS[2]],
             col = colHaplo,
             pch = pchHaplo, cex = cexHaplo,
             xlab = paste0("MDS", plotWhichMDS[1]),
             ylab = paste0("MDS", plotWhichMDS[2]),
             xlim = rangeX,
             ylim = rangeY)
        
        if (existComp) {
          points(x = mdsResComp$points[compNames, plotWhichMDS[1]],
                 y = mdsResComp$points[compNames, plotWhichMDS[2]],
                 col = colComp,
                 pch = pchComp, cex = cexComp)
        }
        
        segments(x0 = mdsResComp$points[mstResComp[, 1], plotWhichMDS[1]],
                 y0 = mdsResComp$points[mstResComp[, 1], plotWhichMDS[2]],
                 x1 = mdsResComp$points[mstResComp[, 2], plotWhichMDS[1]],
                 y1 = mdsResComp$points[mstResComp[, 2], plotWhichMDS[2]],
                 col = colConnection[1],
                 lty = ltyConnection[1],
                 lwd = lwdConnection[1])
        if (plotPlus) {
          segments(x0 = mdsResComp$points[mstResCompPlus[, 1], plotWhichMDS[1]],
                   y0 = mdsResComp$points[mstResCompPlus[, 1], plotWhichMDS[2]],
                   x1 = mdsResComp$points[mstResCompPlus[, 2], plotWhichMDS[1]],
                   y1 = mdsResComp$points[mstResCompPlus[, 2], plotWhichMDS[2]],
                   col = colConnection[2],
                   lty = ltyConnection[2],
                   lwd = lwdConnection[2])
        }
        
        
        
        title(main = paste0(paste(c(colnames(pheno)[2], 
                                    blockName, kernelType), collapse = "_"),
                            " (-log10p: ", round(minuslog10p, 2), ")"))
        if (existComp) {
          if (!is.null(subpopInfo)) {
            legend("topleft", legend = c(paste0(rep(levels(subpopInfo), each = 2),
                                                rep(c(" (gv:+)",  " (gv:-)"), nGrp)),
                                         "Complement (gv:+)", "Complement (gv:-)"),
                   col = c(rep(colHaploBase, each = 2), colCompBase),
                   pch = rep(pchBase, nGrp + 1))
          } else {
            legend("topleft", legend = c("Haplotype (gv:+)", "Haplotype (gv:-)",
                                         "Complement (gv:+)", "Complement (gv:-)"),
                   col = c(rep(colHaploBase, each = 2), colCompBase),
                   pch = rep(pchBase, 2))
          }
        } else if (!is.null(pheno)) {
          if (!is.null(subpopInfo)) {
            legend("topleft", legend = paste0(rep(levels(subpopInfo), each = 2),
                                              rep(c(" (gv:+)",  " (gv:-)"), nGrp)),
                   col = rep(colHaploBase, each = 2),
                   pch = rep(pchBase, nGrp))
          } else {
            legend("topleft", legend = c("Haplotype (gv:+)", "Haplotype (gv:-)"),
                   col = rep(colHaploBase, each = 2),
                   pch = pchBase)
          }
        } else {
          if (!is.null(subpopInfo)) {
            legend("topleft", legend = levels(subpopInfo),
                   col = colHaploBase,
                   pch = pchHaplo)
          } else {
            legend("topleft", legend = "Haplotype",
                   col = colHaploBase,
                   pch = pchHaplo)
          }   
        }
        
        
        
        if (!is.null(saveName)) {
          dev.off()
        } else {
          Sys.sleep(0.5)
        }
      }
      
      
      
      if (ggPlotNetwork) {
        mdsPoints <- data.frame(mdsResComp$points[, plotWhichMDS[1:2]])
        colnames(mdsPoints) <- paste0("MDS", plotWhichMDS[1:2])
        mdsPoints$name <- rownames(mdsPoints)
        
        if (!is.null(subpopInfo)) {
          clusterNosDF <- data.frame(do.call(what = rbind, 
                                             args = clusterNosForHaplotype))
          clusterNosRatioDF <- data.frame(t(apply(clusterNosDF, 1, 
                                                  function(x) x / sum(x))))
          colHaploBaseForPie <- colHaploBase
        } else {
          nGrp <- 1
          colHaploBaseForPie <- "gray"
          
          clusterNosRatioDF <- data.frame(matrix(data = 1,
                                                 nrow = nHaplo,
                                                 ncol = 1))
          rownames(clusterNosRatioDF) <- haploNames
          colnames(clusterNosRatioDF) <- "Haplotype"
        }
        
        
        
        if (existComp) {
          if (!all(is.na(EMMRes))) {
            clusterNosRatioDF$`gvComp +` <- 0
            clusterNosRatioDF$`gvComp -` <- 0
            
            
            colorForCompFac <- factor(gvScaled[compNames] > 0,
                                      levels = c("TRUE", "FALSE"))
            colorForCompDF <- data.frame(model.matrix(~ colorForCompFac - 1))
            rownames(colorForCompDF) <- names(colorForCompFac)
            colnames(colorForCompDF) <- c("gvComp +", "gvComp -")
            
            alphaHaplo <- ifelse(gvScaled[haploNames] > 0, alphaBase[1], alphaBase[2])
            alpha1 <- alphaBase[1]
          } else {
            clusterNosRatioDF$`Complement` <- 0
            
            
            colorForCompDF <- data.frame(matrix(data = NA, 
                                                nrow = nComp, ncol = 1))
            rownames(colorForCompDF) <- compNames
            colnames(colorForCompDF) <- "Complement"
            
            alpha1 <- mean(alphaBase)
            alphaHaplo <- rep(alpha1, nHaplo)
            names(alphaHaplo) <- haploNames
            
            colCompBase <- "gray60"
          }
          
          
          colorForCompDFComp <- data.frame(matrix(data = 0,
                                                  nrow = nComp,
                                                  ncol = nGrp,
                                                  dimnames = list(rownames(colorForCompDF),
                                                                  colnames(clusterNosRatioDF)[1:nGrp])))
          colorForCompDF <- cbind(colorForCompDFComp,
                                  colorForCompDF)
          
          colorForAllCompDF <- data.frame(matrix(data = 0,
                                                 nrow = nTotal,
                                                 ncol = ncol(colorForCompDF)))
          rownames(colorForAllCompDF) <- rownames(mdsPoints)
          colnames(colorForAllCompDF) <- colnames(colorForCompDF)
          colorForAllCompDF[haploNames, ] <- clusterNosRatioDF
          colorForAllCompDF[compNames, ] <- colorForCompDF
          
          
          alphaAll <- rep(NA, nTotal)
          names(alphaAll) <- rownames(mdsPoints)
          alphaAll[haploNames] <- alphaHaplo
          
          alphaComp <- rep(alphaBase[2], nComp)
          alphaAll[compNames] <- alphaComp
        } else {
          colorForAllCompDF <- clusterNosRatioDF
          
          if (EM3Res$weights[2] > 1e-06) {
            alphaHaplo <- ifelse(gvScaled[haploNames] > 0, alphaBase[1], alphaBase[2])
            alpha1 <- alphaBase[1]
          } else {
            alpha1 <- mean(alphaBase)
            alphaHaplo <- rep(alpha1, nHaplo)
            names(alphaHaplo) <- haploNames
          }
          
          alphaAll <- alphaHaplo
          
          colCompBase <- NULL
        }
        alpha2 <- alphaBase[2]
        
        
        cexAll <- max(diff(rangeX), diff(rangeY)) * cexAll
        mdsPointsForPlotDF <- cbind(mdsPoints,
                                    colorForAllCompDF)
        colnames(mdsPointsForPlotDF)[1:2] <- paste0("MDS", 1:2)
        mdsPointsForPlotDF$cex <- cexAll
        
        
        plt <- ggplot2::ggplot(data = mdsPointsForPlotDF,
                               ggplot2::aes(x = .data$MDS1,
                                            y = .data$MDS2)) 
        if (any(alphaAll == alpha2)) {
          plt <- plt +
            scatterpie::geom_scatterpie(data = mdsPointsForPlotDF[alphaAll == alpha1, ],
                                        ggplot2::aes(x = .data$MDS1,
                                                     y = .data$MDS2,
                                                     r = .data$cex),
                                        cols = colnames(colorForAllCompDF),
                                        col = NA, alpha = alpha1) +
            scatterpie::geom_scatterpie(data = mdsPointsForPlotDF[alphaAll == alpha2, ],
                                        ggplot2::aes(x = .data$MDS1,
                                                     y = .data$MDS2,
                                                     r = .data$cex),
                                        cols = colnames(colorForAllCompDF),
                                        col = NA, alpha = alpha2) +
            ggplot2::scale_fill_manual(values = c(colHaploBaseForPie,
                                                  colCompBase)[order(colnames(colorForAllCompDF))]) +
            ggplot2::coord_equal() +
            ggplot2::xlab(paste0("MDS", plotWhichMDS[1])) +
            ggplot2::ylab(paste0("MDS", plotWhichMDS[2]))
        } else {
          plt <- plt +
            scatterpie::geom_scatterpie(data = mdsPointsForPlotDF[alphaAll == alpha1, ],
                                        ggplot2::aes(x = .data$MDS1,
                                                     y = .data$MDS2,
                                                     r = .data$cex),
                                        cols = colnames(colorForAllCompDF),
                                        col = NA, alpha = alpha1) +
            ggplot2::scale_fill_manual(values = c(colHaploBaseForPie,
                                                  colCompBase)[order(colnames(colorForAllCompDF))]) +
            ggplot2::coord_equal() +
            ggplot2::xlab(paste0("MDS", plotWhichMDS[1])) +
            ggplot2::ylab(paste0("MDS", plotWhichMDS[2]))
        }
        
        
        
        mdsSegmentsDF <- data.frame(x1 = mdsPoints[mstResComp[, 1], 1],
                                    y1 = mdsPoints[mstResComp[, 1], 2],
                                    x2 = mdsPoints[mstResComp[, 2], 1],
                                    y2 = mdsPoints[mstResComp[, 2], 2])
        
        plt <- plt + ggplot2::geom_segment(ggplot2::aes(x = .data$x1,
                                                        y = .data$y1,
                                                        xend = .data$x2,
                                                        yend = .data$y2),
                                           data = mdsSegmentsDF,
                                           col = colConnection[1],
                                           lty = ltyConnection[1],
                                           lwd = lwdConnection[1])
        if (plotPlus) {
          mdsSegmentsPlusDF <- data.frame(x1 = mdsPoints[mstResCompPlus[, 1], 1],
                                          y1 = mdsPoints[mstResCompPlus[, 1], 2],
                                          x2 = mdsPoints[mstResCompPlus[, 2], 1],
                                          y2 = mdsPoints[mstResCompPlus[, 2], 2])
          
          plt <- plt + ggplot2::geom_segment(ggplot2::aes(x = .data$x1,
                                                          y = .data$y1,
                                                          xend = .data$x2,
                                                          yend = .data$y2),
                                             data = mdsSegmentsPlusDF,
                                             col = colConnection[2],
                                             lty = ltyConnection[2],
                                             lwd = lwdConnection[2])
        }
        
        
        plt <- plt + ggplot2::ggtitle(label = paste0(paste(c(colnames(pheno)[2], 
                                                             blockName, kernelType), collapse = "_"),
                                                     " (-log10p: ", round(minuslog10p, 2), ")"))
        
        
        if (!is.null(saveName)) {
          savePlotNameBase <- paste0(paste(c(saveName, blockName, kernelType), collapse = "_"),
                                     "_network_ggplot")
          if (saveStyle == "pdf") {
            savePlotName <- paste0(savePlotNameBase, ".pdf")
            pdf(file = savePlotName, width = 15.5, height = 9)
          } else if (saveStyle == "jpg") {
            savePlotName <- paste0(savePlotNameBase, ".jpg")  
            jpeg(filename = savePlotName, width = 1300, height = 700)
          } else if (saveStyle == "tiff") {
            savePlotName <- paste0(savePlotNameBase, ".tiff")
            tiff(filename = savePlotName, width = 1300, height = 700)
          } else {
            savePlotName <- paste0(savePlotNameBase, ".png")
            png(filename = savePlotName, width = 1300, height = 700)
          }
        }
        
        print(plt)
        
        if (!is.null(saveName)) {
          dev.off()
        } else {
          Sys.sleep(0.5)
        }
      }
      
      
      
    }
    
    
    return(list(haplotypeInfo = haplotypeInfo,
                subpopInfo = subpopInfo,
                mstResults = list(mstRes = mstRes,
                                  mstResComp = mstResComp),
                distMats = list(distMat = distMat,
                                distMatComp = distMatComp,
                                laplacianMat = L),
                pValChi2Test = pValChi2Test,
                gvTotal = gvEstTotals,
                gvTotalForLine = gvEstTotalForLines,
                minuslog10p = minuslog10ps,
                hOpts = hOptsList,
                EMMResults = EMMResultsList,
                clusterNosForHaplotype = clusterNosForHaplotype))
  }
  
  if (is.matrix(blockInterest)) {
    return(estNetworkEachBlock(blockInterest = blockInterest,
                               blockName = blockName))
  } else if (is.list(blockInterest)) {
    nTopRes <- length(blockInterest)
    allResultsList <- rep(list(NULL), nTopRes)
    blockNames <- names(blockInterest)
    
    if (is.null(blockNames)) {
      blockNames <- paste0("block_", 1:nTopRes)
    }
    names(allResultsList) <- blockNames
    
    
    for (topNo in 1:nTopRes) {
      allResultsList[[topNo]] <- estNetworkEachBlock(blockInterest = blockInterest[[topNo]],
                                                     blockName = blockNames[topNo])
    }
    
    return(allResultsList)
  }
}










#' Function to plot phylogenetic tree from the estimated results
#'
#' @param estPhyloRes The estimated results of phylogenetic analysis by `estPhylo` function for one
#' @param traitName Name of trait of interest. This will be used in the title of the plots.
#' @param blockName You can specify the haplotype block (or gene set, SNP-set) of interest by the name of haplotype block in `geno`.
#' This will be used in the title of the plots.
#' @param plotTree If TRUE, the function will return the plot of phylogenetic tree.
#' @param subpopInfo The information of subpopulations.
#' @param saveName When drawing any plot, you can save plots in png format. In saveName, you should substitute the name you want to save.
#' When saveName = NULL, the plot is not saved.
#' @param saveStyle This argument specifies how to save the plot of phylogenetic tree.
#' The function offers `png`, `pdf`, `jpg`, and `tiff`.
#' @param pchBase A vector of two integers specifying the plot types for the positive and negative genotypic values respectively.
#' @param colNodeBase A vector of two integers or chracters specifying color of nodes for the positive and negative genotypic values respectively.
#' @param colTipBase A vector of integers or chracters specifying color of tips for the positive and negative genotypic values respectively.
#' The length of the vector should equal to the number of subpopulations.
#' @param cexMax A numeric specifying the maximum point size of the plot.
#' @param cexMin A numeric specifying the minimum point size of the plot.
#' @param edgeColoring If TRUE, the edge branch of phylogenetic tree wiil be colored.
#' @param tipLabel If TRUE, lavels for tips will be shown.
#' @param ggPlotTree If TRUE, the function will return the ggplot version of phylogenetic tree. 
#' It offers the precise information on subgroups for each haplotype.
#' @param cexMaxForGG A numeric specifying the maximum point size of the plot for ggtree,
#'  relative to the range of x and y-axes (0 < cexMaxForGG <= 1).
#' @param cexMinForGG A numeric specifying the minimum point size of the plot for ggtree,
#'  relative to the range of x and y-axes (0 < cexMaxForGG <= 1).
#' @param alphaBase alpha (parameter that indicates the opacity of a geom) for tip with positive / negative effects.
#' alpha for node will be same as the alpha for tip with negative effects.
#'
#' @return Draw plots of phylogenetic tree.
#'
plotPhyloTree <- function(estPhyloRes, traitName = NULL, blockName = NULL, plotTree = TRUE,
                          subpopInfo = estPhyloRes$subpopInfo, saveName = NULL, saveStyle = "png",
                          pchBase = c(1, 16), colNodeBase = c(2, 4), colTipBase = c(3, 5, 6),
                          cexMax = 2, cexMin = 0.7, edgeColoring = TRUE, tipLabel = TRUE,
                          ggPlotTree = FALSE, cexMaxForGG = 0.12, cexMinForGG = 0.06, alphaBase = c(0.9, 0.3)) {
  
  haplotypeInfo <- estPhyloRes$haplotypeInfo
  haploCluster <- haplotypeInfo$haploCluster
  blockInterestUniqueSorted <- haplotypeInfo$haploBlock
  
  nHaplo <- nrow(blockInterestUniqueSorted)
  haploNames <- rownames(blockInterestUniqueSorted)
  haploClusterNow <- match(haploCluster, haploNames)
  lineNames <- names(haploCluster)
  names(haploClusterNow) <- lineNames
  
  nGrp <- length(levels(subpopInfo))
  
  nMrkInBlock <- ncol(blockInterestUniqueSorted)
  njRes <- estPhyloRes$njRes
  distNodes <- estPhyloRes$distMats$distMatNJ
  nTotal <- nrow(distNodes)
  nNode <- nTotal - nHaplo
  
  kernelTypes <- names(estPhyloRes$gvTotal)
  
  for (kernelType in kernelTypes){
    EMMResults <- estPhyloRes$EMMResults[[kernelType]]
    EM3Res <- EMMResults$EM3Res
    EMMRes <- EMMResults$EMMRes
    gvEstTotal <- estPhyloRes$gvTotal[[kernelType]]
    minuslog10p <- estPhyloRes$minuslog10p[kernelType]
    nullPheno <- all(is.na(gvEstTotal))
    
    if (!nullPheno) {
      if (EM3Res$weights[2] <= 1e-06) {
        plotNode <- FALSE
        
        cexTip <- cexMax / 2
        pchTip <- pchBase[2]
        hOpt2 <- NA
        
        cexTipForGG <- rep(cexMaxForGG / sqrt(2), nHaplo)
      } else {
        plotNode <- TRUE
        
        gvNode <- gvEstTotal[(nHaplo + 1):nTotal]
        
        gvCentered <- gvEstTotal - mean(gvEstTotal)
        gvScaled <- gvCentered / sd(gvCentered)
        gvScaled4Cex <- gvScaled * (cexMax - cexMin) / max(abs(gvScaled)) 
        
        cexNode <- (abs(gvScaled4Cex) + cexMin)[(nHaplo + 1):nTotal]
        pchNode <- ifelse(gvScaled[(nHaplo + 1):nTotal] > 0, pchBase[1], pchBase[2])
        colNode <- ifelse(gvScaled[(nHaplo + 1):nTotal] > 0, colNodeBase[1], colNodeBase[2])
        
        cexTip <- (abs(gvScaled4Cex) + cexMin)[1:nHaplo]
        pchTip <- ifelse(gvScaled[1:nHaplo] > 0, pchBase[1], pchBase[2])
        
        
        colorForNodeFac <- factor(gvScaled[(nHaplo + 1):nTotal] > 0,
                                  levels = c("TRUE", "FALSE"))
        colorForNodeDF <- data.frame(model.matrix(~ colorForNodeFac - 1))
        rownames(colorForNodeDF) <- names(colorForNodeFac)
        colnames(colorForNodeDF) <- c("gvNode +", "gvNode -")
        
        colorForNodeDF$node <- rownames(colorForNodeDF)
        
        
        gvScaled4CexForGG <- sign(gvScaled) * sqrt(abs(gvScaled)) * 
          (cexMaxForGG - cexMinForGG) / max(sqrt(abs(gvScaled)))
        
        cexNodeForGG <- (abs(gvScaled4CexForGG) + cexMinForGG)[(nHaplo + 1):nTotal]
        cexTipForGG <- (abs(gvScaled4CexForGG) + cexMinForGG)[1:nHaplo]
        
        alphaTip <- ifelse(gvScaled[1:nHaplo] > 0, alphaBase[1], alphaBase[2])
        alphaNode <- alphaBase[2]
      }
    } else {
      plotNode <- FALSE
      
      cexTip <- cexMax / 2
      pchTip <- pchBase[2]
      
      cexTipForGG <- rep(cexMaxForGG / sqrt(2), nHaplo)
    }
    
    
    if (!is.null(subpopInfo)) {
      if (length(colTipBase) != nGrp) {
        stop("The length of 'colTipBase' should be equal to 'nGrp' or the number of subpopulations!")
      }
      
      colTipNo <- as.numeric(subpopInfo)
      colTip <- colTipBase[colTipNo]
      names(colTipNo) <- names(colTip) <- lineNames
      
      colTip <- tapply(colTip, INDEX = haploCluster, FUN = function(x) {
        as.numeric(names(which.max(table(x) / sum(table(x)))))
      })[haploNames]
      
      clusterNosForHaplotype <- tapply(subpopInfo, INDEX = haploCluster,
                                       FUN = table)[haploNames]
      
      edgeCol <- as.numeric(colTip[njRes$edge[, 2]])
      
      clusterNosDF <- data.frame(do.call(what = rbind, 
                                         args = clusterNosForHaplotype))
      clusterNosRatioDF <- data.frame(t(apply(clusterNosDF, 1, 
                                              function(x) x / sum(x))))
      
      clusterNosRatioDF$node <- 1:nHaplo
    } else {
      colTip <- rep("gray", nHaplo)
      names(colTip) <- haploNames
      colTipBase <- "gray"
      
      edgeCol <- colTip[njRes$edge[, 2]]
      clusterNosForHaplotype <- NA
      
      clusterNosRatioDF <- data.frame(matrix(data = 1,
                                             nrow = nHaplo,
                                             ncol = 1))
      rownames(clusterNosRatioDF) <- haploNames
      colnames(clusterNosRatioDF) <- "Tip"
      clusterNosRatioDF$node <- 1:nHaplo
    }
    
    edgeCol[is.na(edgeCol)] <- "gray"
    
    if (plotTree) {
      if (!is.null(saveName)) {
        savePlotNameBase <- paste0(paste(c(saveName, blockName, kernelType), collapse = "_"),
                                   "_phylogenetic_tree")
        if (saveStyle == "pdf") {
          savePlotName <- paste0(savePlotNameBase, ".pdf")
          pdf(file = savePlotName, width = 12, height = 9)
        } else if (saveStyle == "jpg") {
          savePlotName <- paste0(savePlotNameBase, ".jpg")  
          jpeg(filename = savePlotName, width = 800, height = 600)
        } else if (saveStyle == "tiff") {
          savePlotName <- paste0(savePlotNameBase, ".tiff")
          tiff(filename = savePlotName, width = 800, height = 600)
        } else {
          savePlotName <- paste0(savePlotNameBase, ".png")
          png(filename = savePlotName, width = 800, height = 600)
        }
      }
      
      if (edgeColoring) {
        ape::plot.phylo(njRes, type = "u", show.tip.label = F, edge.color = edgeCol)
      } else {
        ape::plot.phylo(njRes, type = "u", show.tip.label = F, edge.color = "gray30")
      }
      if (plotNode) {
        ape::nodelabels(pch = pchNode, cex = cexNode, col = colNode)
      }
      if (tipLabel) {
        ape::tiplabels(pch = pchTip, cex = cexTip, col = colTip)
      }
      title(main = paste0(paste(c(traitName[2], 
                                  blockName, kernelType), collapse = "_"),
                          " (-log10p: ", round(minuslog10p, 2), ")"))
      if (plotNode) {
        if (!is.null(subpopInfo)) {
          legend("topleft", legend = c(paste0(rep(levels(subpopInfo), each = 2),
                                              rep(c(" (gv:+)",  " (gv:-)"), nGrp)),
                                       "Node (gv:+)", "Node (gv:-)"),
                 col = c(rep(colTipBase, each = 2), colNodeBase),
                 pch = rep(pchBase, nGrp + 1))
        } else {
          legend("topleft", legend = c("Tip (gv:+)", "Tip (gv:-)",
                                       "Node (gv:+)", "Node (gv:-)"),
                 col = c(rep(colTipBase, each = 2), colNodeBase),
                 pch = rep(pchBase, 2))
        }
      } else {
        if (!is.null(subpopInfo)) {
          legend("topleft", legend = levels(subpopInfo),
                 col = colTipBase,
                 pch = pchTip)
        } else {
          legend("topleft", legend = "Tip",
                 col = colTipBase,
                 pch = pchTip)
        }   
      }
      if (!is.null(saveName)) {
        dev.off()
      } else {
        Sys.sleep(0.5)
      }
    }
    
    
    
    
    
    if (ggPlotTree) {
      trPhylo41 <- as(njRes, 'phylo4')
      
      if (edgeColoring) {
        dataForTipCol <- data.frame(color = colTip)
      } else {
        dataForTipCol <- data.frame(color = rep("gray30", nHaplo))
      }
      rownames(dataForTipCol) <- njRes$tip.label
      trPhylo42 <- phylobase::phylo4d(trPhylo41, dataForTipCol)
      
      nodeData <- data.frame(randomTrait = gvEstTotal[(nHaplo + 1):nTotal],
                             color = rep("gray", phylobase::nNodes(trPhylo41)),
                             row.names = phylobase::nodeId(trPhylo41, "internal"))
      phylobase::nodeData(trPhylo42) <- nodeData
      
      
      plt <- ggtree::ggtree(tr = trPhylo42,
                            ggtree::aes(col = I(.data$color)),
                            layout = "equal_angle")
      
      
      if (plotNode) {
        piesNode <- ggtree::nodepie(colorForNodeDF, cols = 1:2,
                                    color = colNodeBase,
                                    alpha = alphaNode)
        
        plt <- plt + ggtree::geom_inset(insets = piesNode,
                                        width = cexNodeForGG,
                                        height = cexNodeForGG)
      }
      
      if (tipLabel) {
        if (!is.null(subpopInfo)) {
          colTipBaseForPie <- colTipBase
        } else {
          colTipBaseForPie <- rep("gray", nGrp)
        }
        
        if (!all(is.na(EMMRes))) {
          piesPlus <- ggtree::nodepie(clusterNosRatioDF[alphaTip == alphaBase[1], ],
                                      cols = 1:nGrp,
                                      color = colTipBaseForPie,
                                      alpha = alphaBase[1])
          piesMinus <- ggtree::nodepie(clusterNosRatioDF[alphaTip == alphaBase[2], ],
                                       cols = 1:nGrp,
                                       color = colTipBaseForPie,
                                       alpha = alphaBase[2])
          
          plt <- plt + ggtree::geom_inset(insets = piesPlus,
                                          width = cexTipForGG[alphaTip == alphaBase[1]],
                                          height = cexTipForGG[alphaTip == alphaBase[1]])
          
          plt <- plt + ggtree::geom_inset(insets = piesMinus,
                                          width = cexTipForGG[alphaTip == alphaBase[2]],
                                          height = cexTipForGG[alphaTip == alphaBase[2]])
        } else {
          piesTip <- ggtree::nodepie(clusterNosRatioDF,
                                     cols = 1:nGrp,
                                     color = colTipBaseForPie,
                                     alpha = mean(alphaBase)) 
          
          plt <- plt + ggtree::geom_inset(insets = piesTip,
                                          width = cexTipForGG,
                                          height = cexTipForGG)
        }
      }
      plt <- plt + ggplot2::ggtitle(label = paste0(paste(c(traitName[2], 
                                                           blockName, kernelType), collapse = "_"),
                                                   " (-log10p: ", round(minuslog10p, 2), ")"))
      
      
      if (!is.null(saveName)) {
        savePlotNameBase <- paste0(paste(c(saveName, blockName, kernelType), collapse = "_"),
                                   "_phylogenetic_ggtree")
        if (saveStyle == "pdf") {
          savePlotName <- paste0(savePlotNameBase, ".pdf")
          pdf(file = savePlotName, width = 12, height = 9)
        } else if (saveStyle == "jpg") {
          savePlotName <- paste0(savePlotNameBase, ".jpg")  
          jpeg(filename = savePlotName, width = 800, height = 600)
        } else if (saveStyle == "tiff") {
          savePlotName <- paste0(savePlotNameBase, ".tiff")
          tiff(filename = savePlotName, width = 800, height = 600)
        } else {
          savePlotName <- paste0(savePlotNameBase, ".png")
          png(filename = savePlotName, width = 800, height = 600)
        }
      }
      
      print(plt)
      
      if (!is.null(saveName)) {
        dev.off()
      } else {
        Sys.sleep(0.5)
      }
    }
  }
  
  
}






#' Function to plot haplotype network from the estimated results
#'
#' @param estNetworkRes The estimated results of haplotype network by `estNetwork` function for one
#' @param traitName Name of trait of interest. This will be used in the title of the plots.
#' @param blockName You can specify the haplotype block (or gene set, SNP-set) of interest by the name of haplotype block in `geno`.
#' This will be used in the title of the plots.
#' @param plotNetwork If TRUE, the function will return the plot of haplotype network.
#' @param subpopInfo The information of subpopulations.
#' @param saveName When drawing any plot, you can save plots in png format. In saveName, you should substitute the name you want to save.
#' When saveName = NULL, the plot is not saved.
#' @param saveStyle This argument specifies how to save the plot of phylogenetic tree.
#' The function offers `png`, `pdf`, `jpg`, and `tiff`.
#' @param plotWhichMDS We will show the MDS (multi-dimensional scaling) plot, 
#' and this argument is a vector of two integers specifying that will define which MDS dimension will be plotted.
#' The first and second integers correspond to the horizontal and vertical axes, respectively.
#' @param colConnection A vector of two integers or characters specifying the colors of connection between nodes for the original and complemented haplotypes, respectively.
#' @param ltyConnection A vector of two characters specifying the line types of connection between nodes for the original and complemented haplotypes, respectively.
#' @param lwdConnection A vector of two integers specifying the line widths of connection between nodes for the original and complemented haplotypes, respectively.
#' @param pchBase A vector of two integers specifying the plot types for the positive and negative genotypic values respectively.
#' @param colCompBase A vector of two integers or characters specifying color of complemented haplotypes for the positive and negative genotypic values respectively.
#' @param colHaploBase A vector of integers or characters specifying color of original haplotypes for the positive and negative genotypic values respectively.
#' The length of the vector should equal to the number of subpopulations.
#' @param cexMax A numeric specifying the maximum point size of the plot.
#' @param cexMin A numeric specifying the minimum point size of the plot.
#' @param ggPlotNetwork If TRUE, the function will return the ggplot version of haplotype network. 
#' It offers the precise information on subgroups for each haplotype.
#' @param cexMaxForGG A numeric specifying the maximum point size of the plot for the ggplot version of haplotype network,
#'  relative to the range of x and y-axes (0 < cexMaxForGG <= 1).
#' @param cexMinForGG A numeric specifying the minimum point size of the plot for the ggplot version of haplotype network,
#'  relative to the range of x and y-axes (0 < cexMaxForGG <= 1).
#' @param alphaBase alpha (parameter that indicates the opacity of a geom) for original haplotype with positive / negative effects.
#' alpha for complemented haplotype will be same as the alpha for original haplotype with negative effects.
#'
#' @return Draw plot of haplotype network.
#'
plotHaploNetwork <- function(estNetworkRes, traitName = NULL, blockName = NULL, 
                             plotNetwork = TRUE, subpopInfo = estNetworkRes$subpopInfo, 
                             saveName = NULL, saveStyle = "png",
                             plotWhichMDS = 1:2, colConnection = c("grey40", "grey60"),
                             ltyConnection = c("solid", "dashed"), lwdConnection = c(1.5, 0.8),
                             pchBase = c(1, 16), colCompBase = c(2, 4),
                             colHaploBase = c(3, 5, 6), cexMax = 2, cexMin = 0.7,
                             ggPlotNetwork = FALSE, cexMaxForGG = 0.025, 
                             cexMinForGG = 0.008, alphaBase = c(0.9, 0.3)) {
  
  haplotypeInfo <- estNetworkRes$haplotypeInfo
  haploCluster <- haplotypeInfo$haploCluster
  blockInterestUniqueSorted <- haplotypeInfo$haploBlock
  blockInterestCompSorted <- haplotypeInfo$haploBlockCompSorted
  
  nHaplo <- nrow(blockInterestUniqueSorted)
  haploNames <- rownames(blockInterestUniqueSorted)
  haploClusterNow <- match(haploCluster, haploNames)
  lineNames <- names(haploCluster)
  names(haploClusterNow) <- lineNames
  distMat <- estNetworkRes$distMats$distMat
  distMatComp <- estNetworkRes$distMats$distMatComp
  
  nGrp <- length(levels(subpopInfo))
  
  nMrkInBlock <- ncol(blockInterestUniqueSorted)
  mstResults <- estNetworkRes$mstResults
  mstResComp <- mstResults$mstResComp
  
  nTotal <- nrow(distMatComp)
  nComp <- nTotal - nHaplo
  
  kernelTypes <- names(estNetworkRes$gvTotal)
  namesBlockInterestComp <- rownames(distMatComp)
  
  if (nComp >= 1){
    existComp <- TRUE
    compNames <- paste0("c", 1:nComp)
  } else {
    existComp <- FALSE
    compNames <- NULL
  }
  
  mstResCompPlus <- attr(mstResComp, "alter.links")
  mstResCompAll <- rbind(mstResComp,
                         mstResCompPlus)
  
  L <- estNetworkRes$distMats$laplacianMat
  plotPlus <- !is.null(mstResCompPlus)
  
  
  
  for (kernelType in kernelTypes){
    EMMResults <- estNetworkRes$EMMResults[[kernelType]]
    EM3Res <- EMMResults$EM3Res
    EMMRes <- EMMResults$EMMRes
    gvEstTotal <- estNetworkRes$gvTotal[[kernelType]]
    minuslog10p <- estNetworkRes$minuslog10p[kernelType]
    nullPheno <- all(is.na(gvEstTotal))    
    
    if (!nullPheno) {
      if (EM3Res$weights[2] <= 1e-06) {
        cexHaplo <- cexMax / 2
        cexComp <- cexMax / 4
        pchHaplo <- pchComp <- pchBase[2]
        colComp <- colCompBase[2]
        
        cexHaploForGG <- cexMaxForGG / sqrt(2)
        cexCompForGG <- cexMaxForGG / 2
        cexAll <- rep(NA, nTotal)
        names(cexAll) <- namesBlockInterestComp
        cexAll[haploNames] <- cexHaploForGG
        if (existComp) {
          cexAll[compNames] <- cexCompForGG
        }
      } else {
        if (existComp){
          gvCentered <- gvEstTotal - mean(gvEstTotal)
          gvScaled <- gvCentered / sd(gvCentered)
          gvScaled4Cex <- gvScaled * (cexMax - cexMin) / max(abs(gvScaled)) 
          
          cexComp <- (abs(gvScaled4Cex) + cexMin)[compNames]
          pchComp <- ifelse(gvScaled[compNames] > 0, pchBase[1], pchBase[2])
          colComp <- ifelse(gvScaled[compNames] > 0, colCompBase[1], colCompBase[2])
          
          cexHaplo <- (abs(gvScaled4Cex) + cexMin)[haploNames]
          pchHaplo <- ifelse(gvScaled[haploNames] > 0, pchBase[1], pchBase[2])
          
          
          gvScaled4CexForGG <- gvScaled * (cexMaxForGG - cexMinForGG) / max(abs(gvScaled))
          cexAll <- abs(gvScaled4CexForGG) + cexMinForGG
          cexHaploForGG <- cexAll[haploNames]
          cexCompForGG <- cexAll[compNames]
        } else {
          gvCentered <- gvEstTotal - mean(gvEstTotal)
          gvScaled <- gvCentered / sd(gvCentered)
          gvScaled4Cex <- gvScaled * (cexMax - cexMin) / max(abs(gvScaled)) 
          
          cexHaplo <- (abs(gvScaled4Cex) + cexMin)[haploNames]
          pchHaplo <- ifelse(gvScaled[haploNames] > 0, pchBase[1], pchBase[2])
          
          cexComp <- cexMax / 4
          colComp <- colCompBase[2]
          pchComp <- pchBase[2]
          
          
          gvScaled4CexForGG <- gvScaled * (cexMaxForGG - cexMinForGG) / max(abs(gvScaled))
          cexAll <- abs(gvScaled4CexForGG) + cexMinForGG
          cexHaploForGG <- cexAll[haploNames]
          cexCompForGG <- cexMaxForGG / 2
        }
      }
    } else {
      cexHaplo <- cexMax / 2
      cexComp <- cexMax / 4
      pchHaplo <- pchComp <- pchBase[2]
      colComp <- colCompBase[2]
      
      cexHaploForGG <- cexMaxForGG / sqrt(2)
      cexCompForGG <- cexMaxForGG / 2
      cexAll <- rep(NA, nTotal)
      names(cexAll) <- namesBlockInterestComp
      cexAll[haploNames] <- cexHaploForGG
      if (existComp) {
        cexAll[compNames] <- cexCompForGG
      }
      
    }
    
    
    if (!is.null(subpopInfo)) {
      if (length(colHaploBase) != nGrp) {
        stop("The length of 'colHaploBase' should be equal to 'nGrp' or the number of subpopulations!")
      }
      
      colHaploNo <- as.numeric(subpopInfo)
      colHaplo <- colHaploBase[colHaploNo]
      names(colHaploNo) <- names(colHaplo) <- lineNames
      
      colHaplo <- tapply(colHaplo, INDEX = haploCluster, FUN = function(x) {
        as.numeric(names(which.max(table(x) / sum(table(x)))))
      })[haploNames]
      
      clusterNosForHaplotype <- tapply(subpopInfo, INDEX = haploCluster,
                                       FUN = table)[haploNames]
    } else {
      colHaplo <- rep("gray", nHaplo)
      names(colHaplo) <- haploNames
      
      clusterNosForHaplotype <- NA
    }
    
    
    nDimMDS <- max(plotWhichMDS)
    mdsResComp <- cmdscale(d = distMatComp,
                           k = nDimMDS,
                           eig = TRUE)
    
    rangeX <- range(mdsResComp$points[, plotWhichMDS[1]])
    rangeY <- range(mdsResComp$points[, plotWhichMDS[2]])
    
    
    if (plotNetwork) {
      if (!is.null(saveName)) {
        savePlotNameBase <- paste0(paste(c(saveName, blockName, kernelType), collapse = "_"),
                                   "_network")
        if (saveStyle == "pdf") {
          savePlotName <- paste0(savePlotNameBase, ".pdf")
          pdf(file = savePlotName, width = 12, height = 9)
        } else if (saveStyle == "jpg") {
          savePlotName <- paste0(savePlotNameBase, ".jpg")  
          jpeg(filename = savePlotName, width = 800, height = 600)
        } else if (saveStyle == "tiff") {
          savePlotName <- paste0(savePlotNameBase, ".tiff")
          tiff(filename = savePlotName, width = 800, height = 600)
        } else {
          savePlotName <- paste0(savePlotNameBase, ".png")
          png(filename = savePlotName, width = 800, height = 600)
        }
      }
      
      rangeX[1] <- rangeX[1] - abs(diff(rangeX)) / 4
      rangeY[2] <- rangeY[2] + abs(diff(rangeY)) / 3
      
      plot(x = mdsResComp$points[haploNames, plotWhichMDS[1]],
           y = mdsResComp$points[haploNames, plotWhichMDS[2]],
           col = colHaplo,
           pch = pchHaplo, cex = cexHaplo,
           xlab = paste0("MDS", plotWhichMDS[1]),
           ylab = paste0("MDS", plotWhichMDS[2]),
           xlim = rangeX,
           ylim = rangeY)
      
      if (existComp) {
        points(x = mdsResComp$points[compNames, plotWhichMDS[1]],
               y = mdsResComp$points[compNames, plotWhichMDS[2]],
               col = colComp,
               pch = pchComp, cex = cexComp)
      }
      
      segments(x0 = mdsResComp$points[mstResComp[, 1], plotWhichMDS[1]],
               y0 = mdsResComp$points[mstResComp[, 1], plotWhichMDS[2]],
               x1 = mdsResComp$points[mstResComp[, 2], plotWhichMDS[1]],
               y1 = mdsResComp$points[mstResComp[, 2], plotWhichMDS[2]],
               col = colConnection[1],
               lty = ltyConnection[1],
               lwd = lwdConnection[1])
      if (plotPlus) {
        segments(x0 = mdsResComp$points[mstResCompPlus[, 1], plotWhichMDS[1]],
                 y0 = mdsResComp$points[mstResCompPlus[, 1], plotWhichMDS[2]],
                 x1 = mdsResComp$points[mstResCompPlus[, 2], plotWhichMDS[1]],
                 y1 = mdsResComp$points[mstResCompPlus[, 2], plotWhichMDS[2]],
                 col = colConnection[2],
                 lty = ltyConnection[2],
                 lwd = lwdConnection[2])
      }
      
      
      
      title(main = paste0(paste(c(traitName, 
                                  blockName, kernelType), collapse = "_"),
                          " (-log10p: ", round(minuslog10p, 2), ")"))
      if (existComp) {
        if (!is.null(subpopInfo)) {
          legend("topleft", legend = c(paste0(rep(levels(subpopInfo), each = 2),
                                              rep(c(" (gv:+)",  " (gv:-)"), nGrp)),
                                       "Complement (gv:+)", "Complement (gv:-)"),
                 col = c(rep(colHaploBase, each = 2), colCompBase),
                 pch = rep(pchBase, nGrp + 1))
        } else {
          legend("topleft", legend = c("Haplotype (gv:+)", "Haplotype (gv:-)",
                                       "Complement (gv:+)", "Complement (gv:-)"),
                 col = c(rep(colHaploBase, each = 2), colCompBase),
                 pch = rep(pchBase, 2))
        }
      } else if (!nullPheno) {
        if (!is.null(subpopInfo)) {
          legend("topleft", legend = paste0(rep(levels(subpopInfo), each = 2),
                                            rep(c(" (gv:+)",  " (gv:-)"), nGrp)),
                 col = rep(colHaploBase, each = 2),
                 pch = rep(pchBase, nGrp))
        } else {
          legend("topleft", legend = c("Haplotype (gv:+)", "Haplotype (gv:-)"),
                 col = rep(colHaploBase, each = 2),
                 pch = pchBase)
        }
      } else {
        if (!is.null(subpopInfo)) {
          legend("topleft", legend = levels(subpopInfo),
                 col = colHaploBase,
                 pch = pchHaplo)
        } else {
          legend("topleft", legend = "Haplotype",
                 col = colHaploBase,
                 pch = pchHaplo)
        }   
      }
      
      
      
      if (!is.null(saveName)) {
        dev.off()
      } else {
        Sys.sleep(0.5)
      }
    }
    
    
    
    if (ggPlotNetwork) {
      mdsPoints <- data.frame(mdsResComp$points[, plotWhichMDS[1:2]])
      colnames(mdsPoints) <- paste0("MDS", plotWhichMDS[1:2])
      mdsPoints$name <- rownames(mdsPoints)
      
      if (!is.null(subpopInfo)) {
        clusterNosDF <- data.frame(do.call(what = rbind, 
                                           args = clusterNosForHaplotype))
        clusterNosRatioDF <- data.frame(t(apply(clusterNosDF, 1, 
                                                function(x) x / sum(x))))
        colHaploBaseForPie <- colHaploBase
      } else {
        nGrp <- 1
        colHaploBaseForPie <- "gray"
        
        clusterNosRatioDF <- data.frame(matrix(data = 1,
                                               nrow = nHaplo,
                                               ncol = 1))
        rownames(clusterNosRatioDF) <- haploNames
        colnames(clusterNosRatioDF) <- "Haplotype"
      }
      
      
      
      if (existComp) {
        if (!all(is.na(EMMRes))) {
          clusterNosRatioDF$`gvComp +` <- 0
          clusterNosRatioDF$`gvComp -` <- 0
          
          
          colorForCompFac <- factor(gvScaled[compNames] > 0,
                                    levels = c("TRUE", "FALSE"))
          colorForCompDF <- data.frame(model.matrix(~ colorForCompFac - 1))
          rownames(colorForCompDF) <- names(colorForCompFac)
          colnames(colorForCompDF) <- c("gvComp +", "gvComp -")
          
          alphaHaplo <- ifelse(gvScaled[haploNames] > 0, alphaBase[1], alphaBase[2])
          alpha1 <- alphaBase[1]
        } else {
          clusterNosRatioDF$`Complement` <- 0
          
          
          colorForCompDF <- data.frame(matrix(data = NA, 
                                              nrow = nComp, ncol = 1))
          rownames(colorForCompDF) <- compNames
          colnames(colorForCompDF) <- "Complement"
          
          alpha1 <- mean(alphaBase)
          alphaHaplo <- rep(alpha1, nHaplo)
          names(alphaHaplo) <- haploNames
          
          colCompBase <- "gray60"
        }
        
        
        colorForCompDFComp <- data.frame(matrix(data = 0,
                                                nrow = nComp,
                                                ncol = nGrp,
                                                dimnames = list(rownames(colorForCompDF),
                                                                colnames(clusterNosRatioDF)[1:nGrp])))
        colorForCompDF <- cbind(colorForCompDFComp,
                                colorForCompDF)
        
        colorForAllCompDF <- data.frame(matrix(data = 0,
                                               nrow = nTotal,
                                               ncol = ncol(colorForCompDF)))
        rownames(colorForAllCompDF) <- rownames(mdsPoints)
        colnames(colorForAllCompDF) <- colnames(colorForCompDF)
        colorForAllCompDF[haploNames, ] <- clusterNosRatioDF
        colorForAllCompDF[compNames, ] <- colorForCompDF
        
        
        alphaAll <- rep(NA, nTotal)
        names(alphaAll) <- rownames(mdsPoints)
        alphaAll[haploNames] <- alphaHaplo
        
        alphaComp <- rep(alphaBase[2], nComp)
        alphaAll[compNames] <- alphaComp
      } else {
        colorForAllCompDF <- clusterNosRatioDF
        
        if (EM3Res$weights[2] > 1e-06) {
          alphaHaplo <- ifelse(gvScaled[haploNames] > 0, alphaBase[1], alphaBase[2])
          alpha1 <- alphaBase[1]
        } else {
          alpha1 <- mean(alphaBase)
          alphaHaplo <- rep(alpha1, nHaplo)
          names(alphaHaplo) <- haploNames
        }
        
        alphaAll <- alphaHaplo
        
        colCompBase <- NULL
      }
      alpha2 <- alphaBase[2]
      
      
      cexAll <- max(diff(rangeX), diff(rangeY)) * cexAll
      mdsPointsForPlotDF <- cbind(mdsPoints,
                                  colorForAllCompDF)
      colnames(mdsPointsForPlotDF)[1:2] <- paste0("MDS", 1:2)
      mdsPointsForPlotDF$cex <- cexAll
      
      
      plt <- ggplot2::ggplot(data = mdsPointsForPlotDF,
                             ggplot2::aes(x = .data$MDS1,
                                          y = .data$MDS2)) 
      if (any(alphaAll == alpha2)) {
        plt <- plt +
          scatterpie::geom_scatterpie(data = mdsPointsForPlotDF[alphaAll == alpha1, ],
                                      ggplot2::aes(x = .data$MDS1,
                                                   y = .data$MDS2,
                                                   r = .data$cex),
                                      cols = colnames(colorForAllCompDF),
                                      col = NA, alpha = alpha1) +
          scatterpie::geom_scatterpie(data = mdsPointsForPlotDF[alphaAll == alpha2, ],
                                      ggplot2::aes(x = .data$MDS1,
                                                   y = .data$MDS2,
                                                   r = .data$cex),
                                      cols = colnames(colorForAllCompDF),
                                      col = NA, alpha = alpha2) +
          ggplot2::scale_fill_manual(values = c(colHaploBaseForPie,
                                                colCompBase)[order(colnames(colorForAllCompDF))]) +
          ggplot2::coord_equal() +
          ggplot2::xlab(paste0("MDS", plotWhichMDS[1])) +
          ggplot2::ylab(paste0("MDS", plotWhichMDS[2]))
      } else {
        plt <- plt +
          scatterpie::geom_scatterpie(data = mdsPointsForPlotDF[alphaAll == alpha1, ],
                                      ggplot2::aes(x = .data$MDS1,
                                                   y = .data$MDS2,
                                                   r = .data$cex),
                                      cols = colnames(colorForAllCompDF),
                                      col = NA, alpha = alpha1) +
          ggplot2::scale_fill_manual(values = c(colHaploBaseForPie,
                                                colCompBase)[order(colnames(colorForAllCompDF))]) +
          ggplot2::coord_equal() +
          ggplot2::xlab(paste0("MDS", plotWhichMDS[1])) +
          ggplot2::ylab(paste0("MDS", plotWhichMDS[2]))
      }
      
      
      
      mdsSegmentsDF <- data.frame(x1 = mdsPoints[mstResComp[, 1], 1],
                                  y1 = mdsPoints[mstResComp[, 1], 2],
                                  x2 = mdsPoints[mstResComp[, 2], 1],
                                  y2 = mdsPoints[mstResComp[, 2], 2])
      
      plt <- plt + ggplot2::geom_segment(ggplot2::aes(x = .data$x1,
                                                      y = .data$y1,
                                                      xend = .data$x2,
                                                      yend = .data$y2),
                                         data = mdsSegmentsDF,
                                         col = colConnection[1],
                                         lty = ltyConnection[1],
                                         lwd = lwdConnection[1])
      if (plotPlus) {
        mdsSegmentsPlusDF <- data.frame(x1 = mdsPoints[mstResCompPlus[, 1], 1],
                                        y1 = mdsPoints[mstResCompPlus[, 1], 2],
                                        x2 = mdsPoints[mstResCompPlus[, 2], 1],
                                        y2 = mdsPoints[mstResCompPlus[, 2], 2])
        
        plt <- plt + ggplot2::geom_segment(ggplot2::aes(x = .data$x1,
                                                        y = .data$y1,
                                                        xend = .data$x2,
                                                        yend = .data$y2),
                                           data = mdsSegmentsPlusDF,
                                           col = colConnection[2],
                                           lty = ltyConnection[2],
                                           lwd = lwdConnection[2])
      }
      
      
      plt <- plt + ggplot2::ggtitle(label = paste0(paste(c(traitName, 
                                                           blockName, kernelType), collapse = "_"),
                                                   " (-log10p: ", round(minuslog10p, 2), ")"))
      
      
      if (!is.null(saveName)) {
        savePlotNameBase <- paste0(paste(c(saveName, blockName, kernelType), collapse = "_"),
                                   "_network_ggplot")
        if (saveStyle == "pdf") {
          savePlotName <- paste0(savePlotNameBase, ".pdf")
          pdf(file = savePlotName, width = 15.5, height = 9)
        } else if (saveStyle == "jpg") {
          savePlotName <- paste0(savePlotNameBase, ".jpg")  
          jpeg(filename = savePlotName, width = 1300, height = 700)
        } else if (saveStyle == "tiff") {
          savePlotName <- paste0(savePlotNameBase, ".tiff")
          tiff(filename = savePlotName, width = 1300, height = 700)
        } else {
          savePlotName <- paste0(savePlotNameBase, ".png")
          png(filename = savePlotName, width = 1300, height = 700)
        }
      }
      
      print(plt)
      
      if (!is.null(saveName)) {
        dev.off()
      } else {
        Sys.sleep(0.5)
      }
    }
  }
}
