#'
#' averagePerBin function from Hervé Pagès
#'
#' This function is an implementation of the code that Hervé Pagès generously provided at
#' http://permalink.gmane.org/gmane.science.biology.informatics.conductor/45929. It works by
#' averaging the specified metadata columns in 'mcolnames' using a given binsize over an input
#' ranges object.
#'
#' @param x A GenomicRanges object with non-NA seqlengths.
#' @param binsize A single positive integer.
#' @param mcolnames Names of numeric metadata columns in 'x' to "average" i.e. to propagate to the 
#' result after averaging them on each bin.
#' @return A GRanges object with: (a) the same seqinfo as 'x', (b) ranges of width 'binsize' 
#' covering all the sequences in 'seqinfo(x)' and (c) the "averaged" metadata columns specified
#' in 'mcolnames'.
#'
#' @importFrom GenomicRanges seqlengths
#' @export
#'
averagePerBin <- function(x, binsize, mcolnames=NULL)
{
  if (!is(x, "GenomicRanges")) {
    stop("'x' must be a GenomicRanges object")
  }
  if (any(is.na(seqlengths(x)))) {
    stop("'seqlengths(x)' contains NAs")
  }
  bins <- IRangesList(lapply(seqlengths(x), 
                             function(seqlen) IRanges(breakInChunks(seqlen, binsize))))
  ans <- as(bins, "GRanges")
  seqinfo(ans) <- seqinfo(x)
  if (is.null(mcolnames)) {
    return(ans)
  }
  mcols(ans) <- DataFrame(lapply(mcols(x)[mcolnames], .averageMCol, x=x, bins=bins))
  ans
}

# 
# Internal helper function for averagePerBin
#
.averageMCol <- function(colname, x, bins)
{
  cvg <- coverage(x, weight=colname)
  views_list <- RleViewsList(lapply(names(cvg), 
                                    function(seqname) Views(cvg[[seqname]], bins[[seqname]])))
  unlist(viewMeans(views_list), use.names=FALSE)
}

