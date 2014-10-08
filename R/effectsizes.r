#'
#' Cliff's Delta
#'
#' This function computes Cliff's Delta for two samples. It measures the effect size in a 
#' non parametric way. It is usually employed when a non-parametric test is employed, and a 
#' parametric effect size estimation, like Cohen's d, would make no sense. This naive implementation 
#' is quite computationally expensive, so it should be used with caution, as it could lead to memory
#' shortage problems.
#'
#' @param x First sample
#' @param y Second sample
#'
cliffDelta <- function(x, y) {
  .Deprecated(msg=paste0('This implementation is very inefficient. As of March, 2014, a package ',
                         'named effsize has been added to CRAN, including a Cliff Delta efficient ',
                         'implementation.'))

  if (!is.numeric(x) || !is.numeric(y)) {
    stop('Input samples should be numeric.')
  }

  return(mean(rowMeans(sign(outer(x, y, FUN="-")))))
}

#'
#' Camblor's AC
#'
#' This function computes Camblor's AC, a metric for evaluating the effect size in a non-parametric
#' way between two samples. It was implemented as a efficient replacement for the Cliff's delta in 
#' this package. Main difference is that AC is contained in the interval [0, 1] and it amounts to 
#' the area of the two densities (one for each of the input samples) intersection. So, a value of 1
#' would indicate that the two distributions are the same, hence there would be no effect size. A
#' value of 0 would indicate no intersection between distributions, and amount to a greater effect
#' size. It is an unsigned metric, so no information about the direction of the effect is provided.
#' 
#' @param x First sample.
#' @param y Second sample.
#' @param n Number of points to use for the densities estimation.
#'
camblorAC <- function(x, y, n=512) {
  
  if (n < 5) {
    stop('Number of points for density estimation (n) must be higher or equal to 5')
  }

  minValue <- min(c(x, y), na.rm=TRUE)
  maxValue <- max(c(x, y), na.rm=TRUE)
  densX <- density(x, n=n, from=minValue, to=maxValue, na.rm=TRUE)
  densY <- density(y, n=n, from=minValue, to=maxValue, na.rm=TRUE)
  densityMaximum <- pmax(densX$y, densY$y)
  densityMinimum <- pmin(densX$y, densY$y)
  delta <- diff(densX$x)[1]
  areaMax <- sum(delta * densityMaximum)
  areaMin <- sum(delta * densityMinimum)
  AC <- areaMin / areaMax 

  if (AC > 1.0) {
    return(1.0)
  } else {
    return(AC)
  }
}

#'
#' Odds Ratio
#'
#' This function computes the odds ratio for a 2x2 matrix. It is assumed that the true positive 
#' element count is located in the main diagonal. If that is not the case, the additional parameter
#' changeRows can be used to reflect that situation.
#'
#' @param x Input 2x2 matrix.
#' @param changeRows boolean indicating if rows have to be exchanged prior to OR computation.
#' @export
#'
oddsRatio <- function(x, changeRows=FALSE) {

  if (!is.matrix(x) || nrow(x) != 2 || ncol(x) != 2) {
    stop('Input must be a 2x2 matrix.')
  }

  if (!is.logical(changeRows) || length(changeRows) != 1) {
    stop('changeRows parameter must be of length 1.')
  }

  if (changeRows) {
    x <- rbind(x[2, ], x[1, ])
  } 
  x <- as.numeric(x)
  return((x[1] * x[4]) / (x[2] * x[3])) 
}

#'
#' Odds Ratio for a given level
#'
#' It is not uncommon to work with nx2 contingency matrices, and it is quite often interesting to 
#' be able to compute the odds ratios at the different levels the discrete variable can take. This
#' function does exactly that. It just collapses a bigger matrix onto a 2x2 matrix, like a binary
#' variable representing the pertenence of the data to the given level. The it computes the odds
#' ratio in the traditional way.
#'
#' @param bigTable The nx2 input matrix.
#' @param aLevel The level to be studied. It must be a valid rowname for the matrix.
#'
oddsRatioLevel <- function(bigTable, aLevel) {

  if (!is.matrix(bigTable) || ncol(bigTable) != 2 || nrow(bigTable) <= 2) {
    stop('Input must be a nx2 matrix, with n > 2.')
  }

  indexLevel <- which(rownames(bigTable) == aLevel)
  collapsedTable <- rbind(bigTable[indexLevel, ], colSums(bigTable[-indexLevel, ]))
  return((collapsedTable[1, 1] * collapsedTable[2, 2]) / 
         (collapsedTable[1, 2] * collapsedTable[2, 1])) 
}

