#' Defunct functions in epilab package
#'
#' The functions or variables listed here are no longer part of R as they are no
#'  longer needed.
#'
#' @name epilab-defunct
#' @rdname epilab-defunct
#' @aliases {epilab-defunct,GenomicRegionCommand,genomicRegionCommand,.executeGenomicRegionCommand,.executeNearestGeneCommand}
#' @docType package
#' @section details:
#' \tabular{rl}{
#'   \code{GenomicRegionCommand} \tab have been replaced by \code{\link{GenRegCommand}}\cr
#'   \code{genomicRegionCommand} \tab have been replaced by \code{\link{genRegCommand}}\cr
#'   \code{NearestGeneCommand} \tab have been replaced by \code{\link{NearestGenCommand}}\cr
#'   \code{nearestGeneCommand} \tab have been replaced by \code{\link{nearestGenCommand}, \link{nearestTSSGeneCommand} and \link{nearestTXGeneCommand}}\cr
#' }
#'




#'
#'
#' @export
#'
setClass('GenomicRegionCommand', contains='AnnotationCommand')

#'
#'
.executeGenomicRegionCommand <- function(command, object) {
  .Defunct('.executeGenRegCommand', package = 'epilab')
  object <- callNextMethod()
  futr <-
    reduce(unlist(fiveUTRsByTranscript(TxDb.Hsapiens.UCSC.hg19.knownGene)))
  exons <- unlist(exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by='tx'))
  exons1 <- reduce(exons[exons$exon_rank == 1])
  exons.no1 <- reduce(exons[exons$exon_rank != 1])
  tss2000 <- flank(exons1, 2000)
  introns <-
    reduce(unlist(intronsByTranscript(TxDb.Hsapiens.UCSC.hg19.knownGene)))
  mcols(object)[[paste0(command@colName, 'Prom')]] <-
    (countOverlaps(object, futr) > 0 |
       countOverlaps(object, exons1) > 0 |
       countOverlaps(object, tss2000) > 0)
  mcols(object)[[paste0(command@colName, 'Intra')]] <-
    (countOverlaps(object, exons.no1) > 0 |
       countOverlaps(object, introns) > 0)
  mcols(object)[[paste0(command@colName, 'Inter')]] <-
    !(mcols(object)[[paste0(command@colName, 'Intra')]] |
        mcols(object)[[paste0(command@colName, 'Prom')]])
  return(object)
}

#'
#'
#' @importFrom GenomicRanges reduce
#' @import GenomicFeatures
#' @importFrom IRanges unlist
#'
setMethod('execute', c('GenomicRegionCommand', 'GRanges'), .executeGenomicRegionCommand)

#'
#'
#' @export
#'
genomicRegionCommand <- function(colName) {
  .Defunct('genRegCommand', package = 'epilab')
  return(new('GenomicRegionCommand', colName=colName))
}


#'
#'
#' @export
#'
setClass('NearestGeneCommand', contains='AnnotationCommand')

#'
#'
#' @export
#'
nearestGeneCommand <- function(colName) {
  .Defunct('nearestGenCommand', package = 'epilab')
  return(new('NearestGeneCommand', colName=colName))
}

#'
#'
#'
.executeNearestGeneCommand <- function(command, object) {
  .Defunct('.executeNearestTSSGeneCommand', package = 'epilab')
  object <- callNextMethod()

  rawTranscriptList <- reduce(transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, 'gene'))
  transcriptsWithoutNames <- unlist(rawTranscriptList)
  transcriptsWithNames <- rawTranscriptList
  names(transcriptsWithNames) <- mget(names(rawTranscriptList), org.Hs.egSYMBOL, ifnotfound=NA)
  transcriptsWithNames <- unlist(transcriptsWithNames)

  overlapDistances <- distanceToNearest(object, resize(transcriptsWithoutNames, 1))
  absoluteDistanceToTSS <- elementMetadata(overlapDistances)$distance
  prec <- precede(object, resize(transcriptsWithoutNames, 1))
  foll <- follow(object, resize(transcriptsWithoutNames, 1))

  mcols(object)[[paste0(command@colName, 'DTSS')]] <-
    ifelse(subjectHits(overlapDistances) == foll, absoluteDistanceToTSS, -absoluteDistanceToTSS)
  mcols(object)[[paste0(command@colName, 'GeneSymbol')]] <-
    names(transcriptsWithNames)[subjectHits(overlapDistances)]
  mcols(object)[[paste0(command@colName, 'GeneId')]] <-
    as.numeric(names(transcriptsWithoutNames)[subjectHits(overlapDistances)])

  return(object)
}

#'
#'
#' @importFrom AnnotationDbi mget
#' @importFrom TxDb.Hsapiens.UCSC.hg19.knownGene TxDb.Hsapiens.UCSC.hg19.knownGene
#' @importFrom org.Hs.eg.db org.Hs.egSYMBOL
#'
setMethod('execute', c('NearestGeneCommand', 'GRanges'), .executeNearestGeneCommand)
