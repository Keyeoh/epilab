#'
#' Base AnnotationCommand
#'
#' Abstract class that implements the execute interface for Command Design 
#' Pattern
#'
setClass('AnnotationCommand',
         representation(colName='character')
         )

#' 
#' Execute generic
#'
#' Generic definition of the execute interface for AnnotationCommand objects
#'
#' @docType methods
#' 
setGeneric('execute', function(object, ranges) standardGeneric('execute'))

#'
#' AnnotationCommand base implementation of execute
#'
#' Implementation of the logic that is common to all AnnotationCommands. For 
#' example, error management and input control.
#'
setMethod('execute', c('AnnotationCommand', 'GRanges'),
          function(object, ranges) {
            if(length(ranges) == 0) {
              stop('Cannot execute annotation command on empty GRanges object.')
            } else {
              return(ranges)
            }
          })

#' 
#' Density of CpG AnnotationCommand
#' 
#' This AnnotationCommand adds a column for the density of CpG around a given
#' genomic region
#' 
setClass('DensCpGCommand',
         representation(windowSize='numeric'),
         prototype(windowSize=2000),
         contains='AnnotationCommand')

#'
#' DensCpGCommand validity check
#'
#' This function checks for the validity of a DensCpGCommand object.
#'
setValidity('DensCpGCommand',
            function(object) {
              return(object@windowSize >= 2)
            })
            
#'
#' DensCpG constructor
#'
#' This function builds a DensCpGCommand with a given windowSize
#'
#' @param windowSize Size of the window centered on the input genomic 
#'  region
#'
densCpGCommand <- function(colName, windowSize=2000) {
  return(new('DensCpGCommand', colName=colName, windowSize=windowSize))
}

#'
#' DensCpG implementation of execute
#f
#' @importFrom Biostrings getSeq vcountPattern
#' @importFrom BSgenome.Hsapiens.UCSC.hg19 BSgenome.Hsapiens.UCSC.hg19
#'
setMethod('execute', c('DensCpGCommand', 'GRanges'),
          function(object, ranges) {
            ranges <- callNextMethod()
            rois <- resize(ranges, object@windowSize, fix='center')
            seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, rois)
            counts <- vcountPattern('CG', seqs)
            mcols(ranges)[[object@colName]] <- counts / (object@windowSize / 2)
            return(ranges)
          })

#' 
#' CPGI Status AnnotationCommand
#' 
#' This AnnotationCommand adds a column for the CpG Island status
#' 
setClass('CPGICommand',
         representation(discardDirection='logical'),
         prototype(discardDirection=FALSE),
         contains='AnnotationCommand')

#'
#' CPGICommand constructor
#'
#' This function builds a CPGICommand
#'
#' @param discardDirection Boolean switch indicating if we should distinguish 
#' between North and South Shelves and Shores
#'
cpgiCommand <- function(colName, discardDirection=FALSE) {
  return(new('CPGICommand', colName=colName, discardDirection=discardDirection))
}

#'
#' CPGICommand implementation of execute
#'
#' @importFrom functional Curry
#'
.executeCPGICommand <- function(object, ranges) {
  ranges <- callNextMethod()
  data(hg19.islands, envir=environment(), package='FDb.InfiniumMethylation.hg19')
  hg19.islands <- get('hg19.islands', envir=environment())

  # Region definition
  # NOTE: The ordering of the regions is important, as there are overlapping
  # regions between the different types.
  cpgi.regions <- list('CGI-N-Shelf'=flank(shift(hg19.islands, -2000), 2000),
                       'CGI-S-Shelf'=flank(shift(hg19.islands, 2000), 2000, 
                                           start=FALSE),
                       'CGI-N-Shore'=flank(hg19.islands, 2000),
                       'CGI-S-Shore'=flank(hg19.islands, 2000, start=FALSE),
                       'CGI'=hg19.islands)
  cpgi.region.kv <- mapply(function(xx, yy) list(list(key=xx, value=yy)),
                           names(cpgi.regions), cpgi.regions)

  # Function that adds information about a region
  local.add.region.indices <- function(cpgi, region.kv, ranged.annot) {
    cpgi[countOverlaps(ranged.annot, region.kv$value) > 0] <- region.kv$key
    return(cpgi) 
  }
  add.region.indices <- Curry(local.add.region.indices,
                              ranged.annot=ranges)

  cpgiStatus <- Reduce(add.region.indices, cpgi.region.kv, NA)
  cpgiStatus[is.na(cpgiStatus)] <- 'Non-CGI'

  if (object@discardDirection) {
    cpgiStatus[grep('Shore', cpgiStatus)] <- 'CGI-Shore'  
    cpgiStatus[grep('Shelf', cpgiStatus)] <- 'CGI-Shelf'  
  } 
  
  mcols(ranges)[[object@colName]] <- cpgiStatus

  return(ranges)
}

setMethod('execute', c('CPGICommand', 'GRanges'), .executeCPGICommand)

#'
#' Gap AnnotationCommand
#'
#' This AnnotationCommand adds a column for the distances to centromeres and
#' telomeres
#'
setClass('GapCommand',
         contains='AnnotationCommand')

#'
#' GapCommand constructor
#' 
#' This function builds a GapCommand with a given column name
#'
gapCommand <- function(colName) {
  return(new('GapCommand', colName=colName))
} 

#' 
#' GapCommand implementation of execute
#' 
#' @importFrom rtracklayer browserSession ucscTableQuery getTable
#'
.executeGapCommand <- function(object, ranges) {
    .distance <- function(x) {
      as.data.frame(x)[[3]]
    }

    ranges <- callNextMethod()

    # Get gap table
    s <- browserSession()
    track.name <- 'gap'
    table.name <- 'gap'

    q <- ucscTableQuery(s, track=track.name, table=table.name)

    # Build GRanges object
    raw.table <- getTable(q)

    # Be careful with the 0-based data!
    gr <- with(raw.table, GRanges(chrom, IRanges(chromStart + 1, chromEnd),
                                  type=type))
    telomeres <- gr[gr$type == 'telomere']
    centromeres <- gr[gr$type == 'centromere']

    d.telomeres <- distanceToNearest(ranges, telomeres)
    d.centromeres <- distanceToNearest(ranges, centromeres)

    mcols(ranges)[[paste0(object@colName, 'Cent')]] <- .distance(d.centromeres)
    mcols(ranges)[[paste0(object@colName, 'Telo')]] <- .distance(d.telomeres)

    return(ranges)
}

setMethod('execute', c('GapCommand', 'GRanges'), .executeGapCommand)

#'
#' Genomic Region AnnotationCommand
#' 
#' This AnnotationCommand adds a column containing the genomic region for 
#' a given probe location
#'
setClass('GenomicRegionCommand',
         contains='AnnotationCommand')

#'
#' GenomicRegionCommand constructor
#'
#' This function builds a GenomicRegionCommand with a given column name
#'
genomicRegionCommand <- function(colName) {
    return(new('GenomicRegionCommand', colName=colName))
}

#' 
#' GenomicRegionCommand implementation of execute
#' 
#' @importFrom GenomicRanges reduce
#' @importFrom GenomicFeatures fiveUTRsByTranscript
#' @importFrom IRanges unlist
#'
.executeGenomicRegionCommand <- function(object, ranges) {
  ranges <- callNextMethod()
  futr <- 
    reduce(unlist(fiveUTRsByTranscript(TxDb.Hsapiens.UCSC.hg19.knownGene)))
  exons <- unlist(exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by='tx'))
  exons1 <- reduce(exons[exons$exon_rank == 1])
  exons.no1 <- reduce(exons[exons$exon_rank != 1])
  tss2000 <- flank(exons1, 2000)
  introns <- 
    reduce(unlist(intronsByTranscript(TxDb.Hsapiens.UCSC.hg19.knownGene)))
  mcols(ranges)[[paste0(object@colName, 'Prom')]] <- 
    (countOverlaps(ranges, futr) > 0 | 
     countOverlaps(ranges, exons1) > 0 | 
     countOverlaps(ranges, tss2000) > 0)
  mcols(ranges)[[paste0(object@colName, 'Intra')]] <- 
    (countOverlaps(ranges, exons.no1) > 0 |
     countOverlaps(ranges, introns) > 0)
  mcols(ranges)[[paste0(object@colName, 'Inter')]] <- 
    !(mcols(ranges)[[paste0(object@colName, 'Intra')]] | 
      mcols(ranges)[[paste0(object@colName, 'Prom')]])
  return(ranges) 
}

setMethod('execute', c('GenomicRegionCommand', 'GRanges'), 
          .executeGenomicRegionCommand)

#'
#' Nearest gene AnnotationCommand
#'
#' This AnnotationCommand adds three columns with information regarding 
#' the nearest TSS and gene
#'
setClass('NearestGeneCommand',
         contains='AnnotationCommand')

#'
#' NearestGeneCommand constructor
#'
#' This function builds a NearestGeneCommand with a given column name
#'
nearestGeneCommand <- function(colName) {
  return(new('NearestGeneCommand', colName=colName))
}

#'
#' NearestGeneCommand implementation of execute
#'
#' @importFrom AnnotationDbi mget
#' @importFrom TxDb.Hsapiens.UCSC.hg19.knownGene TxDb.Hsapiens.UCSC.hg19.knownGene
#' @importFrom org.Hs.eg.db org.Hs.egSYMBOL
#'
.executeNearestGeneCommand <- function(object, ranges) {
  ranges <- callNextMethod()
  txs.i <- reduce(transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, 'gene'))
  txs <- unlist(txs.i)
  txs.wn <- txs.i
  names(txs.wn) <- mget(names(txs.i), org.Hs.egSYMBOL, ifnotfound=NA)
  txs.wn <- unlist(txs.wn)
  annot.gene <- distanceToNearest(ranges, resize(txs, 1), keep.sign=TRUE)
  mcols(ranges)[[paste0(object@colName, 'GeneSymbol')]] <- 
    names(txs.wn)[subjectHits(annot.gene)]
  mcols(ranges)[[paste0(object@colName, 'GeneId')]] <- 
    as.numeric(names(txs)[subjectHits(annot.gene)])
  return(ranges)
}

setMethod('execute', c('NearestGeneCommand', 'GRanges'), 
          .executeNearestGeneCommand)

#' 
#' AnnotationCommandList
#' 
#' This AnnotationCommand contains a list of several AnnotationCommands, and 
#' executes them in order
#' 
setClass('AnnotationCommandList',
         representation(commandList='list'),
         prototype(commandList=list()),
         contains='AnnotationCommand',
         validity=function(object) {
           if (!is.list(object@commandList)) {
             stop('commandList slot must be a list.')
           } else if (!all(sapply(object@commandList, function(xx) is(xx, 'AnnotationCommand')))) {
             stop('All of commandList members must be AnnotationCommand objects.')
           } else {
             return(TRUE)
           }
         }
         )

#'
#' AnnotationCommandList constructor
#'
#' This function builds an AnnotationCommandList from a list of 
#' AnnotationCommand objects
#'
annotationCommandList <- function(...) {
  commandList <- list(...)
  return(new('AnnotationCommandList', colName='', commandList=commandList))
}

#'
#' AnnotationCommandList implementation of execute
#'
.executeAnnotationCommandList <- function(object, ranges) {
  ranges <- callNextMethod()
  newRanges <- Reduce(function(xx, yy) execute(yy, xx), object@commandList, 
                      ranges)
  return(newRanges)
}

setMethod('execute', c('AnnotationCommandList', 'GRanges'), 
          .executeAnnotationCommandList)


