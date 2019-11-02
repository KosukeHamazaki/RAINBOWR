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
