#' Generate pseudo phenotypic values
#'
#' @description This function generates pseudo phenotypic values according to the following formula.
#'
#' \deqn{y = X \beta + Z u + e}
#'
#' where effects of major genes are regarded as fixed effects \eqn{\beta} and
#' polygenetic effects are regarded as random effects \eqn{u}.
#' The variances of \eqn{u} and \eqn{e} are automatically determined by the heritability.
#'
#'
#' @param x A n.sample x n.mark genotype matrix where n.sample is sample size and n.mark is the number of markers.
#' @param sample.sets A n.sample x n.mark genotype matrix. Markers with fixed effects (QTNs) are chosen from sample.sets.
#' If sample.sets = NULL, sample.sets = x.
#' @param candidate If you want to fix QTN postitions, please set the number where SNPs to be fixed are located in your data (so not position).
#' If candidate = NULL, QTNs were randomly sampled from sample.sets or x.
#' @param pos A n.mark x 1 vector. Cumulative position (over chromosomes) of each marker.
#' @param x.par If you don't want to match the sampling population and the genotype data to QTN effects, then use this argument as the latter.
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
#' @param x2 A genotype matrix to calculate additive relationship matrix when Z.ETA = NULL.
#' If Z.ETA = NULL & x2 = NULL, calcGRM(x) will be calculated as kernel matrix.
#' @param num.qtn The number of QTNs
#' @param weight The weights for each QTN by their standard deviations. Negative value is also allowed.
#' @param qtn.effect Additive of dominance for each marker effect. This argument should be the same length as num.qtn.
#' @param prop The proportion of effects of QTNs to polygenetic effects.
#' @param polygene.weight If there are multiple kernels, this argument determines the weights of each kernel effect.
#' @param polygene If polygene = FALSE, pseudo phenotypes with only QTN effects will be generated.
#' @param h2 The wide-sense heritability for generating phenotypes. 0 <= h2 < 1
#' @param h.correction If TRUE, this function will generate phenotypes to match the genomic heritability and "h2".
#' @param seed If seed is not NULL, some fixed phenotypic values will be generated according to set.seed(seed)
#' @param plot If TRUE, boxplot for generated phenotypic values will be drawn.
#' @param saveAt When drawing any plot, you can save plots in png format. In saveAt, you should substitute the name you want to save.
#' When saveAt = NULL, the plot is not saved.
#' @param subpop If there is subpopulation structure, you can draw boxpots divide by subpopulations.
#' n.sample x n.subpop matrix. Please indicate the subpopulation information by (0, 1) for each element.
#' (0 means that line doen't belong to that subpopulation, and 1 means that line belongs to that subpopulation)
#' @param return.all If FALSE, only returns generated phenotypic values.
#' If TRUE, this function will return other information such as positions of candidate QTNs.
#' @param seed.env If TRUE, this function will generate different environment effects every time.
#'
#' @return
#' \describe{
#' \item{trait}{Generated phenotypic values}
#' \item{u}{Generated genotyope values}
#' \item{e}{Generated environmental effects}
#' \item{candidate}{The numbers where QTNs are located in your data (so not position).}
#' \item{qtn.position}{QTN positions}
#' \item{heritability}{Genomic heritability for generated phenotypic values.}
#' }
genetrait <- function(x, sample.sets = NULL, candidate = NULL, pos = NULL, x.par = NULL, ZETA = NULL, x2 = NULL,
                      num.qtn = 3, weight = c(2, 1, 1), qtn.effect = rep("A", num.qtn), prop = 1, polygene.weight = 1, polygene = TRUE, h2 = 0.6, h.correction = FALSE,
                      seed = NULL, plot = TRUE, saveAt = NULL, subpop = NULL, return.all = FALSE, seed.env = TRUE){
  test.effects <- c("A", "D")

  #### The constraint of the heritability ####
  if (h2 < 0 | h2 > 1) {
    stop("The heritability must be in [0, 1]!!!")
  }

  if (h2 == 0) {
    num.qtn <- 0
    polygene <- FALSE
    h.correction <- FALSE
  }

  if (is.null(ZETA)) {
    if (is.null(x2)) {
      K <- calcGRM(genoMat = x)
    } else {
      K <- calcGRM(genoMat = x2)
    }
    
    Z <- diag(nrow(K))
    
    ZETA <- list(A = list(Z = Z, K = K))
  }
  lz <- length(ZETA)
  
  Z <- ZETA[[lz]]$Z
  n <- nrow(Z)
  
  if (sum(num.qtn) != 0) {
    #### If there are some QTNs ####
    if (is.null(candidate)) {
      if (!is.null(seed)) {
        set.seed(seed)
      }
      if (!is.null(pos)) {
        #### If there is map information ####
        if (is.null(sample.sets)) {
          qtn.candidates <-
            sample(1:length(pos), num.qtn)    ### Draw QTNs from genotype data randomly.
        } else {
          qtn.candidates <- NULL
          for (i in 1:ncol(sample.sets)) {
            qtn.candidate <- sample(sample.sets[, i], num.qtn[i])
            while (sum(is.na(qtn.candidate)) != 0) {
              qtn.candidate <- sample(sample.sets[, i], num.qtn[i])
            }
            qtn.candidates <- c(qtn.candidates, qtn.candidate)
          }
        }
      } else {
        #### If there is no map information ####
        if (is.null(sample.sets)) {
          qtn.candidates <-
            sample(1:ncol(x), num.qtn)    ### Draw QTNs from genotype data randomly.
        } else {
          qtn.candidates <- NULL
          for (i in 1:ncol(sample.sets)) {
            qtn.candidate <- sample(sample.sets[, i], num.qtn[i])
            while (sum(is.na(qtn.candidate)) != 0) {
              qtn.candidate <- sample(sample.sets[, i], num.qtn[i])
            }
            qtn.candidates <- c(qtn.candidates, qtn.candidate)
          }
        }
      }
    } else {
      qtn.candidates <- candidate
    }
    if (!is.null(pos)) {
      pos.qtns <- pos[qtn.candidates]
    }

    #### Calculate QTN effects ####
    if(sum(num.qtn) != 1){
      if(is.null(x.par)){
        x.cand <- x[, qtn.candidates]
      }else {
        x.cand <- x.par[, qtn.candidates]
      }
      match.effect <- match(qtn.effect, test.effects)
      x.cand.eff <- x.cand
      x.cand.eff[, match.effect == 2] <- apply(x.cand[, match.effect == 2, drop = F], 2, abs)
      effect.qtns.0 <-  1 / apply(x.cand.eff, 2, sd)
    }else {
      x.cand <- x[, qtn.candidates, drop = F]
      x.cand.eff <- x.cand
      x.cand.eff[, match.effect == 2] <- apply(x.cand[, match.effect == 2, drop = F], 2, abs)
      effect.qtns.0 <-  1 / sd(c(x.cand.eff))
    }
    if(length(weight) != sum(num.qtn)){
      stop("The length of weight should be equal to the number of qtns!!")
    }
    effect.qtns <- effect.qtns.0 * weight

    #### Calculate QTN effects of each lineage ####
    if(sum(num.qtn) != 1){
      g1  <- apply(Z %*% x.cand.eff, 1, function(x) sum(x * effect.qtns))
    }else {
      g1 <- c(Z %*% x.cand.eff) * effect.qtns
    }


    #### Calculate polygenetic effects ####
    if(!is.null(seed)){set.seed(seed)}
    if(polygene){
      g2s <- matrix(NA, nrow = n, ncol = lz)
      for(i in 1:lz){
        Z.now <- ZETA[[i]]$Z
        K.now <- ZETA[[i]]$K
        g2.now  <- MASS::mvrnorm(1, rep(0, nrow(x)), K.now)
        g2s[, i] <- Z.now %*% as.matrix(g2.now)
      }

      g2 <- c(g2s %*% polygene.weight)
      g2_2 <- g2 / (sd(g2) * sqrt(prop) / (sd(g1)))   ## match sd(g1) and sd(g2_2)

      g.all <- g1 + g2_2   ### genotypic value!
    }else {
      g.all <- g1
    }
  }else {
    qtn.candidates <- pos.qtns <- NULL

    if (!is.null(seed)) {
      set.seed(seed)
    }
    if (polygene) {
      g2s <- matrix(NA, nrow = n, ncol = lz)
      for (i in 1:lz) {
        Z.now <- ZETA[[i]]$Z
        K.now <- ZETA[[i]]$K
        g2.now  <- MASS::mvrnorm(1, rep(0, nrow(x)), K.now)
        g2s[, i] <- Z.now %*% as.matrix(g2.now)
      }
      
      g2 <- c(g2s %*% polygene.weight)
      g.all <- g2
    } else {
      g.all <- 0
    }
  }
  
  if (sum(num.qtn) != 0 | polygene) {
    g.all2 <- g.all / sd(g.all)
  } else {
    g.all2 <- 0
  }
  
  
  #### Calculate environmental effects using heritability ####
  if (h2 > 0) {
    wariai <- (1 - h2) / h2
  }
  if (!is.null(seed)) {
    if (!seed.env) {
      seed <- NULL
      set.seed(seed)
    }
  }
  
  if (sum(num.qtn) != 0 | polygene) {
    e <- rnorm(n, 0, sqrt(wariai))
  } else {
    h.correction <- FALSE
    e <- rnorm(n, 0, 1)
  }
  
  
  #### Calculate phenotypic value and the real heritability ####
  trait <- g.all2 + e
  true_h <- var(g.all2) / var(trait)
  
  
  #### Correct the environmental effects to match the real and expected heritability ####
  if (h.correction) {
    e.sca <- e / sqrt(sum(e ^ 2))
    g.sca <- g.all2 / sqrt(sum(g.all2 ^ 2))
    
    e2.dir <-
      e.sca - (sum(e.sca * g.sca) * g.sca)   ### Use Schmidt orthgonalization method
    e2 <- e2.dir * sqrt(wariai) / sd(e2.dir)
    
    trait <- g.all2 + e2
    true_h <- var(g.all2) / var(trait)
  } else {
    e2 <- e
  }
  
  #### Show boxplot of phenotypic values ####
  if (plot) {
    if (!is.null(saveAt)) {
      png(paste(saveAt, "_boxplot.png", sep = ""))
      boxplot(trait)
      dev.off()
      if (!is.null(subpop)) {
        subpop.level <- unique(subpop)
        subpop.max <- max(table(subpop))
        
        trait.subpops <-
          matrix(rep(NA, length(subpop.level) * subpop.max), nrow = subpop.max)
        for (i in 1:length(subpop.level)) {
          trait.subpop <- trait[which(subpop == subpop.level[i])]
          if (length(trait.subpop) < subpop.max) {
            trait.subpops[, i] <-
              c(trait.subpop, rep(NA, subpop.max - length(trait.subpop)))
          } else {
            trait.subpops[, i] <- trait.subpop
          }
        }
        
        png(paste(saveAt, "_boxplot_by_subpop.png", sep = ""))
        boxplot(trait.subpops, names = subpop.level)
        dev.off()
      }
    } else {
      boxplot(trait)
      
      if (!is.null(subpop)) {
        subpop.level <- unique(subpop)
        subpop.max <- max((table(subpop)))
        
        trait.subpops <-
          matrix(rep(NA, length(subpop.level) * subpop.max), nrow = subpop.max)
        for (i in 1:length(subpop.level)) {
          trait.subpop <- trait[which(subpop == subpop.level[i])]
          if (length(trait.subpop) < subpop.max) {
            trait.subpops[, i] <-
              c(trait.subpop, rep(NA, subpop.max - length(trait.subpop)))
          } else {
            trait.subpops[, i] <- trait.subpop
          }
        }
        
        boxplot(trait.subpops, names = subpop.level)
        
      }
      
      
    }
  }
  if (return.all) {
    if (!is.null(pos)) {
      return(
        list(
          trait = trait,
          u = g.all2,
          e = e2,
          candidate = qtn.candidates,
          qtn.position = pos.qtns,
          heritability = true_h
        )
      )
    } else {
      return(
        list(
          trait = trait,
          u = g.all2,
          e = e2,
          candidate = qtn.candidates,
          heritability = true_h
        )
      )
    }
  } else {
    return(trait)
  }
}
