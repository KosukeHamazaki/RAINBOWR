#' Calculate some summary statistics of GWAS (genome-wide association studies) for simulation study
#'
#' @param res Data frame of GWAS results where the first column is the marker names,
#' the second and third column is the chromosome amd map position, and the forth column is -log10(p) for each marker.
#' @param x A N (lines) x M (markers) marker genotype data (matrix), coded as {-1, 0, 1} = {aa, Aa, AA}.
#' @param map.x Data frame with the marker names in the first column. The second and third columns contain the chromosome and map position.
#' @param qtn.candidate A vector of causal markers. You should assign where those causal markers are positioned in our marker genotype,
#' rather than physical position of those causal markers.
#' @param gene.set If you have information of gene (or haplotype block), and if you used it to perform kernel-based GWAS,
#'            you should assign your gene information to gene.set in the form of a "data.frame" (whose dimension is (the number of gene) x 2).
#'            In the first column, you should assign the gene name. And in the second column, you should assign the names of each marker,
#'            which correspond to the marker names of "x" argument.
#' @param n.top.false.block We will calculate the mean of -log10(p) values of top 'n.top.false.block' blocks
#' to evaluate the inflation level of results. The default is 10.
#' @param sig.level Significance level for the threshold. The default is 0.05.
#' @param method.thres Method for detemining threshold of significance. "BH" and "Bonferroni are offered.
#' @param LD_length SNPs within the extent of LD are regareded as one set. This LD_length determines the size of LD block,
#' and 2 x LD_length (b.p.) will be the size of LD block.
#' @param cor.thres SNPs within the extent of LD are regareded as one set. This cor.thres also determines the size of LD block,
#' and the region with square of correlation coefficients >= cor.thres is regareded as one LD block. More precisely, the regions
#' which satisfies both LD_length and cor.thres condition is rearded as one LD block.
#' @param inflator.plus If `the -log10(p) value for each marker` exceeds (`the inflation level` + `inflator.plus`),
#' that marker is regarded as significant.
#' @param window.size If you peform SNP-set analysis with slinding window, we can consider the effect of window size by this argument.
#' @param saveName When drawing any plot, you can save plots in png format. In saveName, you should substitute the name you want to save.
#' When saveName = NULL, the plot is not saved.
#' @param plot.ROC If this argunent is TRUE, ROC (Reciever Operating Characteristic) curve will be drawn with AUC (Area Under the Curve).
#' @return
#' \describe{
#' \item{$log.p}{-log10(p)) values of the causals.}
#' \item{$qtn.logp.order}{The rank of -log10(p) of causals.}
#' \item{$thres}{A vector which contains the information of threshold.}
#' \item{$overthres}{The number of markers which exceed the threshold.}
#' \item{$AUC}{Area under the curve.}
#' \item{$AUC.relax}{Area under the curve calculated with LD block units.}
#' \item{$FDR}{False discovery rate. 1 - Precision.}
#' \item{$FPR}{False positive rate.}
#' \item{$FNR}{False negative rate. 1 - Recall.}
#' \item{$Recall}{The proportion of the number of causals dected by GWAS to the number of causals you set.}
#' \item{$Precision}{The proportion of the number of causals dected by GWAS to the number of markers detected by GWAS.}
#' \item{$Accuracy}{The accuracy of GWAS results.}
#' \item{$Hm}{Harmonic mean of Recall and Precision.}
#' \item{$haplo.name}{The haplotype block name which correspond to causals.}
#' \item{$mean.false}{The mean of -log10(p) values of top 'n.top.false.block' blocks.}
#' \item{$max.trues}{Maximum of the -log10(p) values of the region near causals.}
#' }
#'
#'
#'
SS_gwas <- function(res, x, map.x, qtn.candidate, gene.set = NULL, n.top.false.block = 10,
                    sig.level = c(0.05, 0.01), method.thres = "BH", inflator.plus = 2,
                    LD_length = 150000, cor.thres = 0.35, window.size = 0, saveName = NULL, plot.ROC = TRUE){
  trait.name <- colnames(res)[4]


  if(is.null(gene.set)){
    ###### 1. For normal Results ######
    ##### 1.1. Chromosome and position ######
    res <- res[order(res[, 2], res[, 3]), ]
    marker <- as.character(res[, 1])
    chr <- res[, 2]
    pos <- res[, 3]

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

    ##### 1.2. Calculate theshold #####
    thress <- CalcThreshold(res, sig.level = sig.level, method = method.thres)
    overthress <- sapply(thress, function(x) ifelse(is.na(x), 0, sum(res[, 4] >= x)))


    ##### 1.3. The rank of -log10(p) of QTNs #####
    jun <- order(res[, 4], decreasing = T)
    qtn.logp.order.seg <- match(qtn.candidate, jun)
    True.num <- length(qtn.candidate)
    position.qtn <- cum.pos[qtn.candidate]

    True.position.list <- NULL
    over.qtn.list <- NULL
    for(i in 1:True.num){
      True.position0 <- which(cum.pos >= (position.qtn[i] - LD_length) &
                                cum.pos <= (position.qtn[i] + LD_length))

      cor.position <- cor(x[,qtn.candidate[i]], x[, True.position0]) ^ 2
      True.position <- (True.position0[min(which(cor.position >= cor.thres))] - window.size):(True.position0[max(which(cor.position >= cor.thres))] + window.size)
      True.position <- True.position[True.position >= 1 & True.position <= ncol(x)]

      True.position.list <- c(True.position.list, list(True.position))

      over.qtn <- True.position[which(res[True.position, 4] >= res[qtn.candidate[i], 4])]
      over.qtn.list <- c(over.qtn.list, list(over.qtn))
    }



    #### 1.4. Calculation of AUC and draw ROC curve ####
    over.nums <- rep(NA, True.num)
    overs <- NULL
    for(l in 1:True.num){
      overs <- c(overs, over.qtn.list[[order(qtn.logp.order.seg)[l]]])
      over.nums[l] <- length(unique(overs))
    }


    AUC.TP <- c(0, 1:True.num, True.num) / True.num
    AUC.FP <- (c(0, sort(qtn.logp.order.seg), ncol(x)) -
                 c(0, over.nums, over.nums[True.num])) /
      (ncol(x) - over.nums[True.num])
    AUC.FP[AUC.FP <= 0] <- 0


    AUC <- AUC.TP[2] * AUC.FP[2] / 2

    for(m in 2:(length(AUC.FP) - 1)){
      AUC.segment <- (AUC.TP[m] + AUC.TP[m + 1]) * (AUC.FP[m + 1] - AUC.FP[m]) / 2
      AUC <- AUC + AUC.segment
    }

    if(plot.ROC){
      if(is.null(saveName)){
        plot(AUC.FP, AUC.TP, type="l")
        text(0.83, 0.15, paste("AUC = ", round(AUC, 3), sep = ""), cex = 1.5)
      }else{
        png(paste(saveName, trait.name, "_ROC_curve.png", sep = ""), width = 700)
        plot(AUC.FP, AUC.TP, type="l")
        text(0.83, 0.15, paste("AUC = ", round(AUC, 3), sep = ""), cex = 1.5)
        dev.off()
      }
    }


    #### 1.5. Calculate FPR, FNR, FDR and Harmonic mean of precision and recall ####
    res.overthres <- sapply(thress, function(x) res[which(res[, 4] >= x), 4], simplify = FALSE)
    cand.overthres <- sapply(thress, function(x) which(res[, 4] >= x), simplify = FALSE)

    jun <- order(res[, 4], decreasing = TRUE)
    cand.jun <- sapply(overthress, function(x) jun[0:x], simplify = FALSE)

    count0 <- count1 <- count2 <- count3 <- count4 <- 0
    count.no <- rep(NA, True.num)
    qtn.candidate2s <- NULL
    false.block.log10ps <- NULL
    while(length(jun) >= 1){
      cand.now <- jun[1]
      cand.around0 <- which(cum.pos >= (cum.pos[cand.now] - LD_length) &
                              cum.pos <= (cum.pos[cand.now] + LD_length))


      cor.now <- cor(x[, cand.now], x[, cand.around0]) ^ 2

      cand.around <- (cand.around0[min(which(cor.now >= cor.thres))] - window.size):(cand.around0[max(which(cor.now >= cor.thres))] + window.size)

      delete.here <- jun[!is.na(match(jun, cand.around))]
      jun <- jun[is.na(match(jun, cand.around))]

      cand.delete <- lapply(cand.jun, function(x) x[!is.na(match(x, cand.around))])
      cand.jun <- lapply(cand.jun, function(x) x[is.na(match(x, cand.around))])


      count0 <- count0 + 1

      check.true <- any(!is.na(match(qtn.candidate, cand.around)))
      if(check.true){
        qtn.candidate2 <- qtn.candidate[!is.na(match(qtn.candidate, cand.around))]

        if(any(is.na(match(qtn.candidate2, qtn.candidate2s)))){
          count1 <- count1 + sapply(cand.jun, function(x) ifelse(length(x) >= 1, 1, 0))
          count2 <- count2 + 1
          count.no[which((qtn.candidate %in% qtn.candidate2) &
                           (!(qtn.candidate %in% qtn.candidate2s)))] <- count0
          qtn.candidate2s <- unique(c(qtn.candidate2s, qtn.candidate2))
        }

      }else{
        count3 <- count3 + sapply(cand.jun, function(x) ifelse(length(x) >= 1, 1, 0))
        count4 <- count4 + 1
        false.block.log10ps <- c(false.block.log10ps, res[cand.now, 4])
      }
    }

    TP <- count1
    FP <- count3
    FN <- count2 - TP
    TN <- count0 - TP - FP - FN


    max.trues <- which.max.trues <- rep(NA, True.num)
    for(i in 1:True.num){
      max.trues[i] <- max(res[True.position.list[[i]], 4])
      which.max.trues[i] <- which.max(res[True.position.list[[i]], 4])

    }
    names(max.trues) <- names(which.max.trues) <- paste0("qtn_", 1:True.num)
    mean.false <- sapply(n.top.false.block, function(x) mean(false.block.log10ps[1:x], na.rm = TRUE))
    names(mean.false) <- n.top.false.block

    haplo.name <- NULL


    FDR <- FP / (FP + TP)
    FPR <- FP / (FP + TN)
    FNR <- FN / (TP + FN)
    Recall <- TP / (TP + FN)
    Precision <- TP / (TP + FP)
    Accuracy <- (TP + TN) / (TP + FP + TN + FN)
    Hm <- 2 * Recall * Precision / (Recall + Precision)

    if(any(is.nan(Precision))){
      FDR[is.nan(Precision)] <- NA
      Hm[is.nan(Precision)] <- NA
      Precision[is.nan(Precision)] <- NA
    }

    if(any(Recall == 0 | Precision == 0)){
      Hm[Recall == 0 | Precision == 0] <- 0
    }

    log.p <- res[qtn.candidate, 4]


    TP.other <- sapply(mean.false, function(x) sum(log.p[!duplicated(max.trues)] >= (x + inflator.plus)))
    FP.other <- sapply(mean.false, function(x) sum(false.block.log10ps >= (x + inflator.plus)))
    FN.other <- sapply(TP.other, function(x) length(unique(which.max.trues)) - x)
    TN.other <- rep(count0, length(TP.other)) - TP.other - FP.other - FN.other


    FDR.other <- FP.other/ (FP.other + TP.other)
    FPR.other <- FP.other / (FP.other + TN.other)
    FNR.other <- FN.other / (TP.other + FN.other)
    Recall.other <- TP.other / (TP.other + FN.other)
    Precision.other <- TP.other / (TP.other + FP.other)
    Accuracy.other <- (TP.other + TN.other) / (TP.other + FP.other + TN.other + FN.other)
    Hm.other <- 2 * Recall.other * Precision.other / (Recall.other + Precision.other)

    if(any(is.nan(Precision.other))){
      FDR.other[is.nan(Precision.other)] <- NA
      Hm.other[is.nan(Precision.other)] <- NA
      Precision.other[is.nan(Precision.other)] <- NA
    }

    if(any(Recall.other == 0 | Precision.other == 0)){
      Hm.other[Recall.other == 0 | Precision.other == 0] <- 0
    }





    AUC.relax.TP <- c(0, 1:count2, count2) / count2
    AUC.relax.FP <- c(0, sort(unique(count.no)), count4) / count4
    AUC.relax.FP[AUC.relax.FP <= 0] <- 0


    AUC.relax <- AUC.relax.TP[2] * AUC.relax.FP[2] / 2

    for(m in 2:(length(AUC.relax.FP) - 1)){
      AUC.relax.segment <- (AUC.relax.TP[m] + AUC.relax.TP[m + 1]) * (AUC.relax.FP[m + 1] - AUC.relax.FP[m]) / 2
      AUC.relax <- AUC.relax + AUC.relax.segment
    }

    if(plot.ROC){
      if(is.null(saveName)){
        plot(AUC.relax.FP, AUC.relax.TP, type = "l")
        text(0.83, 0.15, paste("AUC.relax = ", round(AUC.relax, 3), sep = ""), cex = 1.5)
      }else{
        png(paste(saveName, trait.name, "_ROC_curve_relax.png", sep = ""), width = 700)
        plot(AUC.relax.FP, AUC.relax.TP, type="l")
        text(0.83, 0.15, paste("AUC.relax = ", round(AUC.relax, 3), sep = ""), cex = 1.5)
        dev.off()
      }
    }

  }else{
    ###### 2. For SNP-set GWAS results ######
    ##### 2.1. Chromosome and position ####
    marker <- as.character(map.x[, 1])
    chr <- map.x[, 2]
    pos <- map.x[, 3]

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

    gene.name <- unique(gene.set[, 1])
    gene.names <- gene.set[, 1]
    mark.id <- gene.set[, 2]



    ##### 2.2. Calculate theshold ####
    thress <- CalcThreshold(res, sig.level = sig.level, method = method.thres)
    overthress <- sapply(thress, function(x) ifelse(is.na(x), 0, sum(res[, 4] >= x)))



    ##### 2.3. The rank of -log10(p) of QTNs #####
    jun <- order(res[, 4], decreasing = T)
    haplo.name <- gene.set[match(marker[qtn.candidate], mark.id), 1]

    ### If there is no corresponding haplotype block for causals,
    for(j in 1:length(qtn.candidate)){
      if(is.na(haplo.name)[j]){
        cum.pos.geneset <- cum.pos[match(mark.id, marker)]
        nearest <- which.min(abs(cum.pos[qtn.candidate[j]] - cum.pos.geneset))

        haplo.name[j] <- gene.set[nearest, 1]
      }
    }


    haplo.no <- match(haplo.name, gene.name)
    gene.candidate <- unique(haplo.no)
    qtn.logp.order.seg <- match(gene.candidate, jun)

    True.num <- length(gene.candidate)
    position.qtn <- cum.pos[qtn.candidate]

    over.haplo.list <- True.haplo.list <- NULL
    for(i in 1:True.num){
      gene.candidate.now <- gene.candidate[i]
      gene.name.now <- gene.name[gene.candidate.now]
      haplo.position.no <- match(gene.candidate.now, haplo.no)

      True.position0 <- which(cum.pos >= (position.qtn[haplo.position.no] - LD_length) &
                                cum.pos <= (position.qtn[haplo.position.no] + LD_length))

      cor.position <- cor(x[, qtn.candidate[haplo.position.no]], x[, True.position0]) ^ 2
      True.position <- (True.position0[min(which(cor.position >= cor.thres))] - window.size):(True.position0[max(which(cor.position >= cor.thres))] + window.size)
      True.position <- True.position[True.position >= 1 & True.position <= ncol(x)]

      cum.pos.geneset <- cum.pos[match(mark.id, marker)]
      cum.pos.True.min <- min(cum.pos[True.position])
      cum.pos.True.max <- max(cum.pos[True.position])

      haplo.name.min <- gene.set[min(which(cum.pos.geneset - cum.pos.True.min >= 0)), 1]
      haplo.name.max <- gene.set[max(which(cum.pos.True.max - cum.pos.geneset >= 0)), 1]

      True.haplo <- match(haplo.name.min, gene.name):match(haplo.name.max, gene.name)
      True.haplo.list <- c(True.haplo.list, list(True.haplo))

      over.haplo <- True.haplo[which(res[True.haplo, 4] >= res[gene.candidate[i], 4])]
      over.haplo.list <- c(over.haplo.list, list(over.haplo))
    }



    ##### 2.4. Calculation of AUC and draw ROC curve #####
    over.nums <- rep(NA, True.num)
    overs <- NULL
    for(l in 1:True.num){
      overs <- c(overs, over.haplo.list[[order(qtn.logp.order.seg)[l]]])
      over.nums[l] <- length(unique(overs))
    }


    AUC.TP <- c(0, 1:True.num, True.num) / True.num
    AUC.FP <- (c(0, sort(qtn.logp.order.seg), ncol(x)) -
                 c(0, over.nums, over.nums[True.num])) /
      (ncol(x) - over.nums[True.num])
    AUC.FP[AUC.FP <= 0] <- 0


    AUC <- AUC.TP[2] * AUC.FP[2] / 2

    for(m in 2:(length(AUC.FP) - 1)){
      AUC.segment <- (AUC.TP[m] + AUC.TP[m + 1]) * (AUC.FP[m + 1] - AUC.FP[m]) / 2
      AUC <- AUC + AUC.segment
    }

    if(plot.ROC){
      if(is.null(saveName)){
        plot(AUC.FP, AUC.TP, type="l")
        text(0.83, 0.15, paste("AUC = ", round(AUC, 3), sep = ""), cex = 1.5)
      }else{
        png(paste(saveName, trait.name, "_ROC_curve.png", sep = ""), width = 700)
        plot(AUC.FP, AUC.TP, type="l")
        text(0.83, 0.15, paste("AUC = ", round(AUC, 3), sep = ""), cex = 1.5)
        dev.off()
      }
    }


    #### 2.5. Calculate FPR, FNR, FDR and Harmonic mean of precision and recall ####
    res.overthres <- sapply(thress, function(x) res[which(res[, 4] >= x), 4], simplify = FALSE)
    cand.overthres <- sapply(thress, function(x) which(res[, 4] >= x), simplify = FALSE)

    jun <- order(res[, 4], decreasing = TRUE)
    cand.jun <- sapply(overthress, function(x) jun[0:x], simplify = FALSE)

    count0 <- count1 <- count2 <- count3 <- count4 <- 0
    count.no <- rep(NA, True.num)
    gene.candidate2s <- NULL
    false.block.log10ps <- NULL
    while(length(jun) >= 1){
      cand.haplo.now <- jun[1]
      cand.nows <- gene.set[which(gene.names == gene.name[cand.haplo.now]), 2]

      cand.arounds <- NULL
      for(cand.no in 1:length(cand.nows)){
        cand.now <- cand.nows[cand.no]
        cand.around0 <- which(cum.pos >= (cum.pos[cand.now] - LD_length) &
                                cum.pos <= (cum.pos[cand.now] + LD_length))

        cor.now <- cor(x[, cand.now], x[, cand.around0]) ^ 2

        cand.around.now <- (cand.around0[min(which(cor.now >= cor.thres))] - window.size):(cand.around0[max(which(cor.now >= cor.thres))] + window.size)
        cand.arounds <- c(cand.arounds, cand.around.now)
      }
      cand.around <- min(cand.arounds):max(cand.arounds)

      cum.pos.geneset <- cum.pos[match(mark.id, marker)]
      cum.pos.cand.min <- min(cum.pos[cand.around])
      cum.pos.cand.max <- max(cum.pos[cand.around])

      haplo.name.min <- gene.set[min(which(cum.pos.geneset - cum.pos.cand.min >= 0)), 1]
      haplo.name.max <- gene.set[max(which(cum.pos.cand.max - cum.pos.geneset >= 0)), 1]

      cand.haplo.around <- match(haplo.name.min, gene.name):match(haplo.name.max, gene.name)

      delete.here <- jun[!is.na(match(jun, cand.haplo.around))]
      jun <- jun[is.na(match(jun, cand.haplo.around))]

      cand.delete <- lapply(cand.jun, function(x) x[!is.na(match(x, cand.haplo.around))])
      cand.jun <- lapply(cand.jun, function(x) x[is.na(match(x, cand.haplo.around))])

      count0 <- count0 + 1

      check.true <- any(!is.na(match(gene.candidate, cand.haplo.around)))
      if(check.true){
        gene.candidate2 <- gene.candidate[!is.na(match(gene.candidate, cand.haplo.around))]

        if(any(is.na(match(gene.candidate2, gene.candidate2s)))){
          count1 <- count1 + sapply(cand.jun, function(x) ifelse(length(x) >= 1, 1, 0))
          count2 <- count2 + 1
          count.no[which((gene.candidate %in% gene.candidate2) &
                           (!(gene.candidate %in% gene.candidate2s)))] <- count0
          gene.candidate2s <- unique(c(gene.candidate2s, gene.candidate2))
        }
      }else{
        count3 <- count3 + sapply(cand.jun, function(x) ifelse(length(x) >= 1, 1, 0))
        count4 <- count4 + 1
        false.block.log10ps <- c(false.block.log10ps, res[cand.haplo.now, 4])
      }
    }

    TP <- count1
    FP <- count3
    FN <- count2 - TP
    TN <- count0 - TP - FP - FN



    max.trues <- which.max.trues <- rep(NA, True.num)
    for(i in 1:True.num){
      max.trues[i] <- max(res[True.haplo.list[[i]], 4])
      which.max.trues[i] <- which.max(res[True.haplo.list[[i]], 4])
    }
    names(max.trues) <- names(which.max.trues) <- paste0("haplo_", 1:True.num)
    mean.false <- sapply(n.top.false.block, function(x) mean(false.block.log10ps[1:x], na.rm = TRUE))
    names(mean.false) <- n.top.false.block

    FDR <- FP / (FP + TP)
    FPR <- FP / (FP + TN)
    FNR <- FN / (TP + FN)
    Recall <- TP / (TP + FN)
    Precision <- TP / (TP + FP)
    Accuracy <- (TP + TN) / (TP + FP + TN + FN)
    Hm <- 2 * Recall * Precision / (Recall + Precision)

    if(any(is.nan(Precision))){
      FDR[is.nan(Precision)] <- NA
      Hm[is.nan(Precision)] <- NA
      Precision[is.nan(Precision)] <- NA
    }

    if(any(Recall == 0 | Precision == 0)){
      Hm[Recall == 0 | Precision == 0] <- 0
    }
    log.p <- res[gene.candidate, 4]


    TP.other <- sapply(mean.false, function(x) sum(log.p[!duplicated(max.trues)] >= (x + inflator.plus)))
    FP.other <- sapply(mean.false, function(x) sum(false.block.log10ps >= (x + inflator.plus)))
    FN.other <- sapply(TP.other, function(x) length(unique(which.max.trues)) - x)
    TN.other <- rep(count0, length(TP.other)) - TP.other - FP.other - FN.other


    FDR.other <- FP.other/ (FP.other + TP.other)
    FPR.other <- FP.other / (FP.other + TN.other)
    FNR.other <- FN.other / (TP.other + FN.other)
    Recall.other <- TP.other / (TP.other + FN.other)
    Precision.other <- TP.other / (TP.other + FP.other)
    Accuracy.other <- (TP.other + TN.other) / (TP.other + FP.other + TN.other + FN.other)
    Hm.other <- 2 * Recall.other * Precision.other / (Recall.other + Precision.other)

    if(any(is.nan(Precision.other))){
      FDR.other[is.nan(Precision.other)] <- NA
      Hm.other[is.nan(Precision.other)] <- NA
      Precision.other[is.nan(Precision.other)] <- NA
    }

    if(any(Recall.other == 0 | Precision.other == 0)){
      Hm.other[Recall.other == 0 | Precision.other == 0] <- 0
    }






    AUC.relax.TP <- c(0, 1:count2, count2) / count2
    AUC.relax.FP <- c(0, sort(unique(count.no)), count4) / count4
    AUC.relax.FP[AUC.relax.FP <= 0] <- 0


    AUC.relax <- AUC.relax.TP[2] * AUC.relax.FP[2] / 2

    for(m in 2:(length(AUC.relax.FP) - 1)){
      AUC.relax.segment <- (AUC.relax.TP[m] + AUC.relax.TP[m + 1]) * (AUC.relax.FP[m + 1] - AUC.relax.FP[m]) / 2
      AUC.relax <- AUC.relax + AUC.relax.segment
    }


    if(plot.ROC){
      if(is.null(saveName)){
        plot(AUC.relax.FP, AUC.relax.TP, type = "l")
        text(0.83, 0.15, paste("AUC.relax = ", round(AUC.relax, 3), sep = ""), cex = 1.5)
      }else{
        png(paste(saveName, trait.name, "_ROC_curve_relax.png", sep = ""), width = 700)
        plot(AUC.relax.FP, AUC.relax.TP, type="l")
        text(0.83, 0.15, paste("AUC.relax = ", round(AUC.relax, 3), sep = ""), cex = 1.5)
        dev.off()
      }
    }

  }


  return(list(log.p = log.p, qtn.logp.order = qtn.logp.order.seg, thres = thress, overthres = overthress,
              AUC = AUC, AUC.relax = AUC.relax, FDR = FDR, FPR = FPR, FNR = FNR, Recall = Recall, Precision = Precision,
              Accuracy = Accuracy, Hm = Hm, FDR.other = FDR.other, FPR.other = FPR.other, FNR.other = FNR.other,
              Recall.other = Recall.other, Precision.other = Precision.other, Accuracy.other = Accuracy.other,
              Hm.other = Hm.other, haplo.name = haplo.name, mean.false = mean.false, max.trues = max.trues,
              false.block.log10ps = false.block.log10ps))
}
