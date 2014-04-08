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
#' @encoding ISO-8859-2
#'
averagePerBin <- function(x, binsize, mcolnames)
{
  if (!is(x, "GenomicRanges")) {
    stop('x must be a GenomicRanges object')
  }
  if (length(x) == 0) {
    stop('x must be a non-empty object.')
  }
  if (any(is.na(seqlengths(x)))) {
    stop('seqlengths(x) contains NAs.')
  }
  if (!is.character(mcolnames)) {
    stop('mcolnames must be a character vector.')
  }
  bins <- IRangesList(lapply(seqlengths(x), 
                             function(seqlen) IRanges(breakInChunks(seqlen, binsize))))
  ans <- as(bins, 'GRanges')
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

#'
#' Generate CIRCOS suitable data from a GRanges (or a subset)
#'
#' This function is just a wrapper in order to generate data in the format that CIRCOS accepts. It 
#' allows to generate data from the whole GRanges object, or from a subset of it. It also accepts
#' a numerical vector of values to use as scores of the generated track information.
#'
#' @param ranges The input GRanges object.
#' @param ids Character vector indicating the names of the elements to be generated.
#' @param values Numerical vector providing the corresponding elements' scores.
#'
#' @export
#'
generateCircosFromRanges <- function(ranges, ids=names(ranges), values=NULL) {

  if (length(ranges) == 0) {
    stop('Input ranges must be a non-empty object.')
  }
  if (!is.character(ids)) {
    stop('Element ids must be of character type')
  }
  if (!is.null(values) && !is.numeric(values)) {
    stop('Element values must be of numerical type')
  }
  selectedRanges <- ranges[ids]
  selectedDf <- as(selectedRanges, 'data.frame')[, 1:3]
  selectedDf$seqnames <- gsub('chr', 'hs', selectedDf$seqnames)
  
  if (!is.null(values)) {
    selectedDf$values <- values
  }

  return(selectedDf)
}

#'
#' Fill sequence length information
#'
#' This function fills the seqlength field of a given GRanges object with information extracted 
#' from several sources, mainly a BSgenome annotation package. This is just a convenience function
#' for some of our workflows, and it currently supports only hg19 annotation package.
#'
#' @param ranges A GRanges input object.
#' @param type Character value indicating the annotation to use.
#' @return The same GRanges with the seqlengths field updated for the given annotation.
#'
#' @export
#'
updateSeqLengthsFromBSGenome <- function(ranges, type=c('hg19')) {
  type <- match.arg(type)

  if (type == 'hg19') {
    sLengths <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)
  } else {
    stop('BSgenome identifier not supported.')
  }

  chrs <- names(seqlengths(ranges))
  seqlengths(ranges) <- sLengths[chrs]

  return(ranges)
}

#'
#' Generate CIRCOS density data from GRanges
#'
#' This function generates a data.frame in CIRCOS compatible format from an input GRanges object.
#' In this case, an averaging step is performed in order to generate a windowed, averaged version
#' of a given input variable associated with the elements of the input genomic regions. It is also
#' possible to execute the function over a subset of element identifiers, and to specify one or
#' more metadata columns for averaging. It is useful for the display of genome-wide scale Circos
#' diagrams, when we want to see the information with a fixed zoom level.
#'
#' @param ranges An input GRanges object. It must have sequence length information.
#' @param ids Character vector indicating the names of the elements to average over.
#' @param wsize Numeric value indicating the size of the window used for averaging.
#' @param mcolname Character vector indicating the name of the column for averaging.
#' @return A data.frame with the averaged information in Circos friendly format.
#'
#' @export
#'
generateCircosDensityFromRanges <- function(ranges, ids=names(ranges), wsize=1e6, mcolname=NULL) {
  if (is.null(mcolname)) {
    mcols(ranges)$dummy <- 1
    mcolname <- 'dummy'
  }
  if (is.null(ids)) {
    selectedRanges <- ranges
  } else {
    selectedRanges <- ranges[ids]
  }
  result <- averagePerBin(selectedRanges, binsize=wsize, mcolnames=mcolname)

  resultDf <- as(result, 'data.frame')[, c(1:3, 6)]

  levels(resultDf$seqnames) <- gsub('chr', 'hs', levels(resultDf$seqnames))
  
  return(resultDf)
}

#'
#' Write CIRCOS data
#'
#' This is just a wrapper around write.table to ensure that CIRCOS data gets written to disk in good
#' shape. Just a way of grouping the formatting options and isolating the side effects for testing.
#'
#' @param circosData A data.frame containing information in CIRCOS format.
#' @param filename Name of the file to write the information to.
#'
#' @export
#'
writeCircos <- function(circosData, filename) {
  write.table(circosData, file=filename, quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
}

