#' Defunct functions in epilab

#'
#' Genomic Region AnnotationCommand
#'
#' This AnnotationCommand adds a column containing the genomic region for a given probe location.
#'
#' @export
#' @rdname GenomicRegionCommand-defunct
#'
setClass('GenomicRegionCommand', contains='AnnotationCommand')

#'
#' Genomic Region AnnotationCommand
#'
#' This class is not valid from now
#'
#' @export
#' @rdname GenomicRegionCommand-defunct
setValidity('GenomicRegionCommand', function(object) {
  # .Defunct('GenRegCommand', package = 'epilab')
  TRUE
})

#
# Internal GenomicRegionCommand implementation of execute
#
#' @rdname .executeGenomicRegionCommand-defunct
#'
.executeGenomicRegionCommand <- function(command, object) {
  # .Defunct(package = 'epilab')
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
#' GenomicRegionCommand implementation of execute
#'
#' GenomicRegionCommand tries to label the input genomic regions according to their relative
#' position with respect to the TSS. We define the Promoter region as the union of a 2kbp region
#' upstream the TSS, the first exon and the 5'UTR. Intragenic region is defined as the union of
#' the remaining exons and introns, those regions inside a gene which are not assigned to the
#' previous definition of Promoter region. Finally, Intergenic region is assigned when neither of
#' the former labels can be applied. This process is always executed at a transcript level, so
#' it is possible for an input region to be at Promoter and Intragenic regions at the same time.
#' Again, this command is better suited for very small input regions.
#'
#' @importFrom GenomicRanges reduce
#' @import GenomicFeatures
#' @importFrom IRanges unlist
#' @param command A GenomicRegionCommand.
#' @param object A GRanges object containing the genomic regions to annotate.
#' @rdname .executeGenomicRegionCommand-defunct
#'
setMethod('execute', c('GenomicRegionCommand', 'GRanges'), .executeGenomicRegionCommand)

#'
#' GenomicRegionCommand constructor
#'
#' This function builds a GenomicRegionCommand with a given column name.
#'
#' @export
#' @param colName Prefix used in order to generate the column name for the annotation.
#' @rdname .genomicRegionCommand-defunct
#'
genomicRegionCommand <- function(colName) {
  # .Defunct('genRegCommand', package = 'epilab')
  return(new('GenomicRegionCommand', colName=colName))
}


#'
#' Nearest gene AnnotationCommand
#'
#' This AnnotationCommand adds columns with information regarding the nearest TSS and gene.
#'
#' @export
#' @rdname NearestGeneCommand-defunct
#'
setClass('NearestGeneCommand', contains='AnnotationCommand')

#'
#' NearestGeneCommand constructor
#'
#' This function builds a NearestGeneCommand with a given column name.
#'
#' @export
#' @param colName Prefix used in order to generate the column name for the annotation.
#' @rdname nearestGeneCommand-defunct
#'
nearestGeneCommand <- function(colName) {
  # .Defunct('nearestGenCommand', package = 'epilab')
  return(new('NearestGeneCommand', colName=colName))
}

#
# Internal NearestGeneCommand implementation of execute
#
#' @rdname .executeNearestGeneCommand-defunct
.executeNearestGeneCommand <- function(command, object) {
  # .Defunct('.executeNearestGenCommand', package = 'epilab')
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
#' NearestGeneCommand implementation of execute
#'
#' NearestGeneCommand labels each input genomic region with information regarding to the nearest
#' gene. Gene information is obtained from the TxDb.Hsapiens.UCSC.hg19.knownGene transcripts
#' database. A gene region is defined as the union of all its transcript regions. The gene symbol
#' is obtained from the org.Hs.eg.db package.
#'
#' @importFrom AnnotationDbi mget
#' @importFrom TxDb.Hsapiens.UCSC.hg19.knownGene TxDb.Hsapiens.UCSC.hg19.knownGene
#' @importFrom org.Hs.eg.db org.Hs.egSYMBOL
#' @param command A NearestGeneCommand.
#' @param object A GRanges object containing the genomic regions to annotate.
#' @rdname NearestGeneCommand-defunct
#'
setMethod('execute', c('NearestGeneCommand', 'GRanges'), .executeNearestGeneCommand)
