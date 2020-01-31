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
          end.other <- narray
          
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
#' @param blockInterest
#' @param gwasRes
#' @param nTopRes
#' @param gene.set
#' @param indexRegion
#' @param chrInterest
#' @param posRegion
#' @param pheno
#' @param geno
#' @param ZETA
#' @param plotTree
#' @param distMat
#' @param distMethod
#' @param evolutionDist
#' @param subpopInfo
#' @param groupingMethod
#' @param nGrp
#' @param nIterClustering
#' @param saveName
#' @param saveStyle
#' @param pchBase
#' @param colNodeBase
#' @param colTipCand
#' @param cexMax
#' @param edgeColoring
#' @param tipLabel
#'
#' @return
#' \describe{
#' \item{$gvTotal}{a}
#' \item{$minuslog10p}{a}
#'}
#'
estPhylo <- function(blockInterest = NULL, gwasRes = NULL, nTopRes = 1, gene.set = NULL,
                     indexRegion = 1:10, chrInterest = NULL, posRegion = NULL,
                     pheno, geno, ZETA = NULL, plotTree = TRUE,
                     distMat = NULL, distMethod = "manhattan", evolutionDist = FALSE,
                     subpopInfo = NULL, groupingMethod = "kmeans",
                     nGrp = 4, nIterClustering = 100,
                     saveName = NULL, saveStyle = "png",
                     pchBase = c(1, 16), colNodeBase = c(2, 4),
                     colTipCand = 2:8, cexMax = 2,
                     edgeColoring = TRUE, tipLabel = TRUE) {
  M <- t(geno[, -c(1:3)])
  map <- geno[, 1:3]
  
  nLine <- nrow(M)
  lineNames <- rownames(M)
  
  if (is.null(ZETA)) {
    K <- A.mat(M)
    Z <- diag(nLine)
    
    rownames(Z) <- colnames(Z) <- lineNames
    
    ZETA <- list(A = list(Z = Z, K = K))
  }
  
  
  if (is.null(subpopInfo)) {
    print("Now clustering...")
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
        kmResNow <- pam(x = M, k = nGrp)
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
  }
  
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
    blockName <- NULL
  } else {
    blockName <- NULL
  }
  
  
  if (is.matrix(blockInterest)) {
    nMrkInBlock <- ncol(blockInterest)
    
    if (is.null(distMat)) {
      distMat <- as.matrix(dist(blockInterest, method = distMethod))
    }
    
    if (evolutionDist) {
      inLog <- 1 - 4 * (distMat / (ncol(blockInterest) * 2)) / 3
      
      inLog[inLog <= 0] <- min(inLog[inLog > 0]) / 10
      
      dist4Nj <- - 3 * log(inLog) / 4
    } else {
      dist4Nj <- distMat
    }
    
    njRes <- nj(X = as.dist(dist4Nj))
    
    
    distNodes <- dist.nodes(njRes) / nMrkInBlock
    hInv <- median(distNodes[upper.tri(distNodes)])
    h <- 1 / hInv
    
    gKernel <- exp(- h * distNodes)
    rownames(gKernel)[1:nLine] <- colnames(gKernel)[1:nLine] <- lineNames
    nTotal <- nrow(gKernel)
    nNode <- nTotal - nLine
    
    gKernelPart <- gKernel[1:nLine, 1:nLine]
    
    ZgKernelPart <- diag(nLine)
    rownames(ZgKernelPart) <- colnames(ZgKernelPart) <- lineNames
    
    ZETANow <- c(ZETA, list(Part = list(Z = ZgKernelPart, K = gKernelPart)))
    
    print("Now estimating genotypic values...")
    EM3Res <- EM3.cpp(y = pheno[, 2], ZETA = ZETANow)
    gvEst <- EM3Res$u[(nLine + 1):(2 * nLine), ]
    LL <- EM3Res$LL
    EMMRes0 <- EMM.cpp(y = pheno[, 2], ZETA = ZETA)
    LL0 <- EMMRes0$LL
    
    pVal <- pchisq(2 * (LL - LL0), df = 1, lower.tail = FALSE)
    minuslog10p <- - log10(pVal)
    
    ZgKernel <- diag(nTotal)
    rownames(ZgKernel) <- colnames(ZgKernel) <- rownames(gKernel)
    
    ZETA2 <- list(Part = list(Z = ZgKernel, K = gKernel))
    gvEst2 <- matrix(c(gvEst, rep(NA, nNode)))
    rownames(gvEst2) <- rownames(gKernel)
    
    EMMRes <- EMM.cpp(y = gvEst2, ZETA = ZETA2)
    gvNode <- EMMRes$u[(nLine + 1):nTotal]
    gvEstTotal <- c(gvEst, gvNode)
    names(gvEstTotal) <- rownames(gKernel)
    gvCentered <- gvEstTotal - mean(gvEstTotal)
    gvScaled <- gvCentered / sd(gvCentered)
    gvScaled4Cex <- gvScaled * cexMax / max(abs(gvScaled)) 
    
    
    cexNode <- abs(gvScaled4Cex)[(nLine + 1):nTotal]
    pchNode <- ifelse(gvNode > 0, pchBase[1], pchBase[2])
    colNode <- ifelse(gvNode > 0, colNodeBase[1], colNodeBase[2])
    
    cexTip <- abs(gvScaled4Cex)[1:nLine]
    pchTip <- ifelse(gvEst > 0, pchBase[1], pchBase[2])
    
    if (!is.null(subpopInfo)) {
      colTip <- as.numeric(subpopInfo) + 1
      names(colTip) <- lineNames
      colUsed <- unique(c(colTip, colNodeBase))
      colTipCand2 <- colTipCand[!(colTipCand %in% colUsed)]
      
      colTip[colTip == colNodeBase[1]] <- colTipCand2[1]
      colTip[colTip == colNodeBase[2]] <- colTipCand2[2]
      edgeCol <- as.numeric(colTip[njRes$edge[, 2]])
    } else {
      colTip <- rep("gray", nLine)
      edgeCol <- colTip[njRes$edge[, 2]]
    }
    
    edgeCol[is.na(edgeCol)] <- "gray"
    
    
    if (!is.null(saveName)) {
      if (saveStyle == "pdf") {
        savePlotName <- paste0(saveName, "_phylogenetic_tree.pdf")
        pdf(file = savePlotName, width = 12, height = 9)
      } else if (saveStyle == "jpg") {
        savePlotName <- paste0(saveName, "_phylogenetic_tree.jpg")  
        jpeg(file = savePlotName, width = 800, height = 600)
      } else if (saveStyle == "tiff") {
        savePlotName <- paste0(saveName, "_phylogenetic_tree.tiff")
        tiff(file = savePlotName, width = 800, height = 600)
      } else {
        savePlotName <- paste0(saveName, "_phylogenetic_tree.png")
        png(file = savePlotName, width = 800, height = 600)
      }
    }
    
    if (edgeColoring) {
      plot.phylo(njRes, type = "u", show.tip.label = F, edge.color = edgeCol)
    } else {
      plot.phylo(njRes, type = "u", show.tip.label = F, edge.color = "gray30")
    }
    nodelabels(pch = pchNode, cex = cexNode, col = colNode)
    if (tipLabel) {
      tiplabels(pch = pchTip, cex = cexTip, col = colTip)
    }
    title(main = paste0(paste(c(colnames(pheno)[2], 
                                blockName), collapse = "_"),
                        " (-log10p: ", round(minuslog10p, 2), ")"))
    if (!is.null(subpopInfo)) {
      legend("topleft", legend = c(paste0(rep(unique(subpopInfo), each = 2),
                                          rep(c(" (gv:+)",  " (gv:-)"), nGrp)),
                                   "Node (gv:+)", "Node (gv:-)"),
             col = c(rep(unique(colTip), each = 2), colNodeBase),
             pch = rep(pchBase, nGrp + 1))
    } else {
      legend("topleft", legend = c("Tip (gv:+)", "Tip (gv:-)",
                                   "Node (gv:+)", "Node (gv:-)"),
             col = c(rep(unique(colTip), each = 2), colNodeBase),
             pch = rep(pchBase, 2))
    }
    if (!is.null(saveName)) {
      dev.off()
    }
    
    
    return(list(gvTotal = gvEstTotal,
                minuslog10p = minuslog10p))
  } else {
    gvEstTotals <- matrix(NA, nrow = nLine, ncol = nTopRes)
    minuslog10ps <- rep(NA, nTopRes)
    colnames(gvEstTotals) <- names(minuslog10ps) <- blockNames
    rownames(gvEstTotals) <- lineNames
    
    for (topNo in 1:nTopRes) {
      blockName <- blockNames[topNo]
      print(paste0("block No: ", topNo, "; ", blockName))
      
      blockInterestNow <- blockInterest[[topNo]]
      
      nMrkInBlock <- ncol(blockInterestNow)
      
      if (is.null(distMat)) {
        distMat <- as.matrix(dist(blockInterestNow, method = distMethod))
      }
      
      if (evolutionDist) {
        inLog <- 1 - 4 * (distMat / (ncol(blockInterestNow) * 2)) / 3
        
        inLog[inLog <= 0] <- min(inLog[inLog > 0]) / 10
        
        dist4Nj <- - 3 * log(inLog) / 4
      } else {
        dist4Nj <- distMat
      }
      
      njRes <- nj(X = as.dist(dist4Nj))
      
      
      distNodes <- dist.nodes(njRes) / nMrkInBlock
      hInv <- median(distNodes[upper.tri(distNodes)])
      h <- 1 / hInv
      
      gKernel <- exp(- h * distNodes)
      rownames(gKernel)[1:nLine] <- colnames(gKernel)[1:nLine] <- lineNames
      nTotal <- nrow(gKernel)
      nNode <- nTotal - nLine
      
      gKernelPart <- gKernel[1:nLine, 1:nLine]
      
      ZgKernelPart <- diag(nLine)
      rownames(ZgKernelPart) <- colnames(ZgKernelPart) <- lineNames
      
      ZETANow <- c(ZETA, list(Part = list(Z = ZgKernelPart, K = gKernelPart)))
      
      print("Now estimating genotypic values...")
      EM3Res <- EM3.cpp(y = pheno[, 2], ZETA = ZETANow)
      gvEst <- EM3Res$u[(nLine + 1):(2 * nLine), ]
      LL <- EM3Res$LL
      EMMRes0 <- EMM.cpp(y = pheno[, 2], ZETA = ZETA)
      LL0 <- EMMRes0$LL
      
      pVal <- pchisq(2 * (LL - LL0), df = 1, lower.tail = FALSE)
      minuslog10p <- - log10(pVal)
      
      ZgKernel <- diag(nTotal)
      rownames(ZgKernel) <- colnames(ZgKernel) <- rownames(gKernel)
      
      ZETA2 <- list(Part = list(Z = ZgKernel, K = gKernel))
      gvEst2 <- matrix(c(gvEst, rep(NA, nNode)))
      rownames(gvEst2) <- rownames(gKernel)
      
      EMMRes <- EMM.cpp(y = gvEst2, ZETA = ZETA2)
      gvNode <- EMMRes$u[(nLine + 1):nTotal]
      gvEstTotal <- c(gvEst, gvNode)
      names(gvEstTotal) <- rownames(gKernel)
      gvCentered <- gvEstTotal - mean(gvEstTotal)
      gvScaled <- gvCentered / sd(gvCentered)
      gvScaled4Cex <- gvScaled * cexMax / max(abs(gvScaled)) 
      
      
      cexNode <- abs(gvScaled4Cex)[(nLine + 1):nTotal]
      pchNode <- ifelse(gvNode > 0, pchBase[1], pchBase[2])
      colNode <- ifelse(gvNode > 0, colNodeBase[1], colNodeBase[2])
      
      cexTip <- abs(gvScaled4Cex)[1:nLine]
      pchTip <- ifelse(gvEst > 0, pchBase[1], pchBase[2])
      
      if (!is.null(subpopInfo)) {
        colTip <- as.numeric(subpopInfo) + 1
        names(colTip) <- lineNames
        colUsed <- unique(c(colTip, colNodeBase))
        colTipCand2 <- colTipCand[!(colTipCand %in% colUsed)]
        
        colTip[colTip == colNodeBase[1]] <- colTipCand2[1]
        colTip[colTip == colNodeBase[2]] <- colTipCand2[2]
        edgeCol <- as.numeric(colTip[njRes$edge[, 2]])
      } else {
        colTip <- rep("gray", nLine)
        edgeCol <- colTip[njRes$edge[, 2]]
      }
      
      edgeCol[is.na(edgeCol)] <- "gray"
      
      if (!is.null(saveName)) {
        if (saveStyle == "pdf") {
          savePlotName <- paste0(saveName, "_phylogenetic_tree.pdf")
          pdf(file = savePlotName, width = 12, height = 9)
        } else if (saveStyle == "jpg") {
          savePlotName <- paste0(saveName, "_phylogenetic_tree.jpg")  
          jpeg(file = savePlotName, width = 800, height = 600)
        } else if (saveStyle == "tiff") {
          savePlotName <- paste0(saveName, "_phylogenetic_tree.tiff")
          tiff(file = savePlotName, width = 800, height = 600)
        } else {
          savePlotName <- paste0(saveName, "_phylogenetic_tree.png")
          png(file = savePlotName, width = 800, height = 600)
        }
      }
      
      if (edgeColoring) {
        plot.phylo(njRes, type = "u", show.tip.label = F, edge.color = edgeCol)
      } else {
        plot.phylo(njRes, type = "u", show.tip.label = F, edge.color = "gray30")
      }
      nodelabels(pch = pchNode, cex = cexNode, col = colNode)
      if (tipLabel) {
        tiplabels(pch = pchTip, cex = cexTip, col = colTip)
      }
      title(main = paste0(paste(c(colnames(pheno)[2], 
                                  blockName), collapse = "_"),
                          " (-log10p: ", round(minuslog10p, 2), ")"))
      if (!is.null(subpopInfo)) {
        legend("topleft", legend = c(paste0(rep(unique(subpopInfo), each = 2),
                                            rep(c(" (gv:+)",  " (gv:-)"), nGrp)),
                                     "Node (gv:+)", "Node (gv:-)"),
               col = c(rep(unique(colTip), each = 2), colNodeBase),
               pch = rep(pchBase, nGrp + 1))
      } else {
        legend("topleft", legend = c("Tip (gv:+)", "Tip (gv:-)",
                                     "Node (gv:+)", "Node (gv:-)"),
               col = c(rep(unique(colTip), each = 2), colNodeBase),
               pch = rep(pchBase, 2))
      }
      
      if (!is.null(saveName)) {
        dev.off()
      }
      
      gvEstTotals[, topNo] <- gvEstTotal
      minuslog10ps[topNo] <- minuslog10p
    }
    
    return(list(gvTotal = gvEstTotals,
                minuslog10p = minuslog10ps))
  }
}
