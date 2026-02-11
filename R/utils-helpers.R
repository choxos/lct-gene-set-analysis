#' Convert Gene Set DataFrame to List Format
#'
#' @description
#' Internal helper function to convert gene sets from data frame format
#' (binary matrix with genes as rows and sets as columns) to list format
#' (named list of character vectors). If input is already a list, returns
#' it unchanged.
#'
#' @param GS Gene set as dataframe (genes x sets with binary indicators)
#'   or list of character vectors
#'
#' @return List of gene sets where each element is a character vector of
#'   gene identifiers belonging to that set
#'
#' @keywords internal
#' @noRd
GS.format.dataframe.to.list <- function(GS) {
  if (is.data.frame(GS)) {
    genes <- rownames(GS)
    L <- NULL
    for (ags in names(GS)) {
      w <- which(GS[, ags] == 1)
      if (length(w) > 0) {
        L <- c(L, list(genes[w]))
        names(L)[length(L)] <- ags
      }
    }
    return(L)
  } else {
    return(GS)
  }
}

#' Calculate T²-like Test Statistic
#'
#' @description
#' Internal function to compute the sum of squared mean differences between
#' two phenotype groups. This is the T²-like test statistic used in the
#' Linear Combination Test after eigenvalue transformation.
#'
#' @param DATA Expression data matrix (genes x samples)
#' @param cl Factor or vector defining phenotype groups (must have exactly
#'   two levels)
#'
#' @return Numeric value: sum of squared mean differences across genes
#'
#' @details
#' For a gene set with p genes, computes:
#' \deqn{\sum_{i=1}^{p} (\bar{x}_{i1} - \bar{x}_{i2})^2}
#' where \eqn{\bar{x}_{i1}} and \eqn{\bar{x}_{i2}} are the mean expression
#' levels of gene i in groups 1 and 2, respectively.
#'
#' @keywords internal
#' @noRd
T2.like.SAMGS <- function(DATA, cl) {
  if (!is.null(cl)) {
    cl <- as.factor(cl)
    C1 <- which(cl == levels(cl)[1])
    C2 <- which(cl == levels(cl)[2])
  }

  if (is.null(C1) || is.null(C2)) {
    stop("Error - T2.like.SAMGS: classes 1 and 2 are undefined. ",
         "Phenotype must have exactly two levels.")
  }

  GS.size <- dim(DATA)[1]

  if (GS.size >= 2) {
    means.C1 <- rowMeans(DATA[, C1, drop = FALSE])
    means.C2 <- rowMeans(DATA[, C2, drop = FALSE])
    diffmean.C1C2 <- means.C1 - means.C2
  } else {
    means.C1 <- mean(DATA[C1])
    means.C2 <- mean(DATA[C2])
    diffmean.C1C2 <- means.C1 - means.C2
  }

  return(sum(diffmean.C1C2^2))
}
