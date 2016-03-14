#'
#' Base AnnotationCommand
#'
#' Abstract class that implements the execute interface for Command Design Pattern.
#'
#' @slot colName Prefix used in order to generate the column name for the annotation.
#'
setClass('AnnotationCommand', representation(colName='character'))

#'
#' getColName generic
#'
#' Generic definition of the getColName getter methods.
#'
#' @param object An object that probably contains a colName slot.
#'
setGeneric('getColName', function(object) standardGeneric('getColName'))

#'
#' AnnotationCommand base implementation of getColName
#'
#' Getter method for colName slot for all AnnotationCommand objects.
#'
#' @param object An AnnotationCommand object.
#'
setMethod('getColName', 'AnnotationCommand',
          function(object) {
            return(object@colName)
          })

#'
#' setColName generic
#'
#' Generic definition of the setColName setter methods.
#'
#' @param object An object that probably contains a colName slot.
#' @param value The new column name.
#'
setGeneric('setColName<-', function(object, value) standardGeneric('setColName<-'))

#'
#' AnnotationCommand base implementation of setColName
#'
#' Getter method for colName slot for all AnnotationCommand objects.
#'
#' @param object An AnnotationCommand object.
#' @param value The new column name.
#'
setReplaceMethod('setColName', 'AnnotationCommand',
          function(object, value) {
            object@colName <- value
            validObject(object)
            return(object)
          })

#'
#' Execute generic
#'
#' Generic definition of the execute interface for a command on an object. An abstraction for
#' implementing hierarchies of commands based on the Command Design Pattern
#'
#' @export
#' @param command A command acting on an object and returning the same kind of object.
#' @param object The object the command acts on.
#'
setGeneric('execute', function(command, object) standardGeneric('execute'))

#'
#' AnnotationCommand base implementation of execute
#'
#' Implementation of the logic that is common to all AnnotationCommands. For
#' example, error management and input control. For now, it only prevents the command from being
#' executed on an empty object.
#'
#' @param command An AnnotationCommand object.
#' @param object A GRanges object containing the genomic regions to annotate.
#'
setMethod('execute', c('AnnotationCommand', 'GRanges'),
          function(command, object) {
            if(length(object) == 0) {
              stop('Cannot execute annotation command on empty GRanges command.')
            } else {
              print(class(command))
              return(object)
            }
          })

#'
#' Density of CpG AnnotationCommand
#'
#' This AnnotationCommand adds a column for the density of CpG around a given genomic region. The
#' size of the region is variable and controlled by the user.
#'
#' @export
#' @slot windowSize Size of the window (centered on the probe) for which the density is computed.
#'
setClass('DensCpGCommand',
         representation(windowSize='numeric'),
         prototype(windowSize=2000),
         contains='AnnotationCommand',
         validity=function(object) {
           return(object@windowSize >= 2)
         }
         )

#'
#' DensCpG constructor
#'
#' This function builds a DensCpGCommand with a given windowSize.
#'
#' @param colName Prefix used in order to generate the column name for the annotation.
#' @param windowSize Size of the window centered on the input genomic region.
#' @export
#'
densCpGCommand <- function(colName, windowSize=2000) {
  return(new('DensCpGCommand', colName=colName, windowSize=windowSize))
}

#'
#' getWindowSize generic
#'
#' Generic definition of the getWindowSize getter methods.
#'
#' @param object An object that probably contains a windowSize slot.
#'
setGeneric('getWindowSize', function(object) standardGeneric('getWindowSize'))

#'
#' DensCpGCommand base implementation of getWindowSize
#'
#' Getter method for windowSize slot for all AnnotationCommand objects.
#'
#' @param object An DensCpGCommand object.
#'
setMethod('getWindowSize', 'DensCpGCommand',
          function(object) {
            return(object@windowSize)
          })

#'
#' DensCpG implementation of execute
#'
#' DensCpGCommand computes the density of CpG in a neighbourhood of the genomic input regions.
#' First, it resizes the input regions to the value of the windowSize slot. Be careful, however,
#' that in the case that the windowSize is smaller than the size of the input region the latter is
#' going to be narrowed by this action. For now, this is enough, because our input regions are
#' mostly GRanges objects containing methylation probes from Illumina microarrays. But it should be
#' fixed in order to use it for annotating NGS peaks, for example.
#'
#' @importFrom Biostrings getSeq vcountPattern
#' @importFrom BSgenome.Hsapiens.UCSC.hg19 BSgenome.Hsapiens.UCSC.hg19
#' @param command A DensCpGCommand.
#' @param object A GRanges object containing the genomic regions to annotate.
#'
setMethod('execute', c('DensCpGCommand', 'GRanges'),
          function(command, object) {
            object <- callNextMethod()
            rois <- resize(object, command@windowSize, fix='center')
            seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, rois)
            counts <- vcountPattern('CG', seqs)
            mcols(object)[[command@colName]] <- counts / (command@windowSize / 2)
            return(object)
          })

#'
#' CPGI Status AnnotationCommand
#'
#' This AnnotationCommand adds a column for the CpG Island status.
#'
#' @export
#' @slot discardDirection Logical value controlling whether the upstream(north)-downstream(south)
#' status of a CpG island-related region is important or not.
#'
setClass('CPGICommand',
         representation(discardDirection='logical'),
         prototype(discardDirection=FALSE),
         contains='AnnotationCommand')

#'
#' CPGICommand constructor
#'
#' This function builds a CPGICommand.
#'
#' @param colName Prefix used in order to generate the column name for the annotation.
#' @param discardDirection Logical value controlling whether the upstream(north)-downstream(south)
#' status of a CpG island-related region is important or not.
#' @export
#'
cpgiCommand <- function(colName, discardDirection=FALSE) {
  return(new('CPGICommand', colName=colName, discardDirection=discardDirection))
}

#'
#' getDiscardDirection generic
#'
#' Generic definition of the getDiscardDirection getter methods.
#'
#' @param object An object that probably contains a discardDirection slot.
#'
setGeneric('getDiscardDirection', function(object) standardGeneric('getDiscardDirection'))

#'
#' CPGICommand base implementation of getDiscardDirection
#'
#' Getter method for discardDirection slot for all CPGICommand objects.
#'
#' @param object An CPGICommand object.
#'
setMethod('getDiscardDirection', 'CPGICommand',
          function(object) {
            return(object@discardDirection)
          })

#
# Internal CPGICommand implementation of execute
#
.executeCPGICommand <- function(command, object) {
  object <- callNextMethod()
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
                              ranged.annot=object)

  cpgiStatus <- Reduce(add.region.indices, cpgi.region.kv, NA)
  cpgiStatus[is.na(cpgiStatus)] <- 'Non-CGI'

  if (command@discardDirection) {
    cpgiStatus[grep('Shore', cpgiStatus)] <- 'CGI-Shore'
    cpgiStatus[grep('Shelf', cpgiStatus)] <- 'CGI-Shelf'
  }

  mcols(object)[[command@colName]] <- cpgiStatus

  return(object)
}

#'
#' CPGICommand implementation of execute
#'
#' CPGICommand uses the CpG island definitions of Irizarry et al. in order to label the input
#' genomic regions according to their island status. Shores are defined as both 2kbp regions
#' flanking the island, and Shelves as the two external 2kbp regions flanking the Shores. As in the
#' DensCpGCommand, this command is best suited for the annotation of very short (methylation probes)
#' input regions.
#'
#' @importFrom functional Curry
#' @import GenomicRanges
#' @param command A CPGICommand.
#' @param object A GRanges object containing the genomic regions to annotate.
#'
setMethod('execute', c('CPGICommand', 'GRanges'), .executeCPGICommand)

#'
#' Gap AnnotationCommand
#'
#' This AnnotationCommand adds a column for the distances to centromeres and telomeres.
#'
#' @export
#'
setClass('GapCommand', contains='AnnotationCommand')

#'
#' GapCommand constructor
#'
#' This function builds a GapCommand with a given column name.
#'
#' @export
#' @param colName Prefix used in order to generate the column name for the annotation.
#'
gapCommand <- function(colName) {
  return(new('GapCommand', colName=colName))
}

#
# Internal GapCommand implementation of execute.
#
.executeGapCommand <- function(command, object) {
    .distance <- function(x) {
      if(length(x) > 0) {
        as.data.frame(x)[[3]] }
    }

    object <- callNextMethod()

    # Get gap table
    s <- browserSession()
    track.name <- 'gap'
    table.name <- 'gap'

    q <- ucscTableQuery(s, track=track.name, table=table.name)

    # Build GRanges command
    raw.table <- getTable(q)

    # Be careful with the 0-based data!
    gr <- GRanges(raw.table$chrom, IRanges(raw.table$chromStart + 1, raw.table$chromEnd),
                                  type=raw.table$type)
    telomeres <- gr[gr$type == 'telomere']
    centromeres <- gr[gr$type == 'centromere']

    d.telomeres <- distanceToNearest(object, telomeres)
    d.centromeres <- distanceToNearest(object, centromeres)

    mcols(object)[[paste0(command@colName, 'Cent')]] <- .distance(d.centromeres)
    mcols(object)[[paste0(command@colName, 'Telo')]] <- .distance(d.telomeres)

    return(object)
}

#'
#' GapCommand implementation of execute
#'
#' GapCommand uses rtracklayer for connecting to the UCSC database and get the information
#' associated with the Gap track, which contains information about centromeres and telomeres.
#' Distance to the nearest telomere and centromere is then computed for the input regions.
#'
#' @importFrom rtracklayer browserSession ucscTableQuery getTable
#' @param command A GapCommand.
#' @param object A GRanges object containing the genomic regions to annotate.
#'
setMethod('execute', c('GapCommand', 'GRanges'), .executeGapCommand)

#'
#' Genomic Region AnnotationCommand
#'
#' This AnnotationCommand adds a column containing the genomic region for a given probe location.
#'
#' @export
#'
setClass('GenRegCommand', contains='AnnotationCommand')

#
# Internal GenRegCommand implementation of execute
#
.executeGenRegCommand <- function(command, object) {
  object <- callNextMethod()
  futr <- unlist(fiveUTRsByTranscript(TxDb.Hsapiens.UCSC.hg19.knownGene))
  exons <- unlist(exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by='tx'))
  exons1 <- exons[exons$exon_rank == 1]
  exons.no1 <- exons[exons$exon_rank != 1]
  tss2000 <- flank(exons1, 2000)
  introns <- unlist(intronsByTranscript(TxDb.Hsapiens.UCSC.hg19.knownGene))
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
#' GenRegCommand constructor
#'
#' This function builds a GenRegCommand with a given column name.
#'
#' @export
#' @param colName Prefix used in order to generate the column name for the annotation.
#'
genRegCommand <- function(colName) {
  return(new('GenRegCommand', colName = colName))
}

#'
#' GenRegCommand implementation of execute
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
#'
setMethod('execute', c('GenRegCommand', 'GRanges'), .executeGenRegCommand)

#'
#' Disjoint Genomic Region AnnotationCommand
#'
#' This AnnotationCommand adds a column containing the genomic region for a given probe location. It
#' is different from the GenomicRegionCommand in that it generates disjoint categories and is more
#' specific.
#'
#' @export
#'
setClass('DGenomicRegionCommand', contains='AnnotationCommand')

#'
#' DGenomicRegionCommand constructor
#'
#' This function builds a DGenomicRegionCommand with a given column name.
#'
#' @export
#' @param colName Prefix used in order to generate the column name for the annotation.
#'
dGenomicRegionCommand <- function(colName) {
    return(new('DGenomicRegionCommand', colName=colName))
}

#
# Internal DGenomicRegionCommand implementation of execute
#
.executeDGenomicRegionCommand <- function(command, object) {
  object <- callNextMethod()

  futr <- unlist(fiveUTRsByTranscript(TxDb.Hsapiens.UCSC.hg19.knownGene))
  tutr <- unlist(threeUTRsByTranscript(TxDb.Hsapiens.UCSC.hg19.knownGene))
  exons <- unlist(exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by='tx'))
  exons1 <- exons[exons$exon_rank == 1]
  exons.no1 <- exons[exons$exon_rank != 1]
  tss2000 <- flank(exons1, 2000)
  introns <- unlist(intronsByTranscript(TxDb.Hsapiens.UCSC.hg19.knownGene))

  mcols(object)[[command@colName]] <- 'Intergenic'
  mcols(object)[[command@colName]][countOverlaps(object, introns) > 0] <- 'Intron'
  mcols(object)[[command@colName]][countOverlaps(object, exons.no1) > 0] <- 'Exon'
  mcols(object)[[command@colName]][countOverlaps(object, exons1) > 0] <- 'FirstExon'
  mcols(object)[[command@colName]][countOverlaps(object, tutr) > 0] <- '3UTR'
  mcols(object)[[command@colName]][countOverlaps(object, futr) > 0] <- '5UTR'
  mcols(object)[[command@colName]][countOverlaps(object, tss2000) > 0] <- 'Promoter'

  return(object)
}

#'
#' DGenomicRegionCommand implementation of execute
#'
#' DGenomicRegionCommand tries to label the input genomic regions according to their relative
#' position with respect to the TSS. Regions defined include the Promoter region (2kbp upstream the
#' TSS), 5'UTR, First Exon, a different Exon, an Intron and the 3'UTR. Finally, Intergenic region is
#' assigned when neither of the former labels can be applied. This process is always executed at a
#' transcript level.
#'
#' @importFrom GenomicRanges reduce
#' @import GenomicFeatures
#' @importFrom IRanges unlist
#' @param command A GenomicRegionCommand.
#' @param object A GRanges object containing the genomic regions to annotate.
#'
setMethod('execute', c('DGenomicRegionCommand', 'GRanges'), .executeDGenomicRegionCommand)


#'
#' Nearest gene AnnotationCommand
#'
#' This AnnotationCommand adds columns with information regarding the nearest TSS, TX and gene.
#'
#' @export
#'
setClass('NearestGenCommand', representation(what = 'character'), contains='AnnotationCommand')

#'
#' nearestTSSGeneCommand constructor
#'
#' This function builds a NearestGenCommand with a given column name and the type of nearest (TSS)
#'
#' @export
#' @param colName Prefix used in order to generate the column name for the annotation.
#'
nearestTSSGeneCommand <- function(colName) {
  return(new('NearestGenCommand', colName=colName, what = 'tss'))
}

#'
#' nearestTXGeneCommand constructor
#'
#' This function builds a NearestGenCommand with a given column name and the type of nearest (transcript)
#'
#' @export
#' @param colName Prefix used in order to generate the column name for the annotation.
#'
nearestTXGeneCommand <- function(colName) {
  return(new('NearestGenCommand', colName=colName, what = 'transcript'))
}

#'
#' nearestGenCommand constructor
#'
#' This function builds a NearestGenCommand with a given column name and the type of nearest (gene)
#'
#' @export
#' @param colName Prefix used in order to generate the column name for the annotation.
#'
nearestGenCommand <- function(colName) {
  return(new('NearestGenCommand', colName=colName, what = 'gene'))
}


#
# Internal NearestGenCommand implementation of execute
#
.executeNearestGenCommand <- function(command, object) {
  object <- callNextMethod()

  if (!is(object, "GenomicRanges")) {
    stop("This function expects a GenomicRanges object as the first argument")
  }

  grl = transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, "gene")
  gene_id_to_gene_symbol = mget(names(grl), org.Hs.egSYMBOL, ifnotfound = NA)

  if (tolower(command@what) == "gene") {
    gr <- unlist(reduce(grl))
  }
  else if (tolower(command@what) == "transcript") {
    gr <- unlist(grl)
  }
  else if (tolower(command@what) == "tss") {
    gr <- unique(resize(unlist(grl), 1, fix = "start"))
  }
  object_nearest = as(distanceToNearest(object, gr), "data.frame")
  rownames(object_nearest) = names(object)[object_nearest$queryHits]
  foll = follow(object, resize(gr, 1))

  mcols(object)[[paste0(command@colName, 'Distance')]] = ifelse(object_nearest$subjectHits == foll, object_nearest$distance, -object_nearest$distance)
  mcols(object)[[paste0(command@colName, 'GeneSymbol')]] = unlist(gene_id_to_gene_symbol[names(gr[object_nearest$subjectHits])])
  mcols(object)[[paste0(command@colName, 'GeneId')]] = as.numeric(names(gr[object_nearest$subjectHits]))

  if (tolower(command@what) %in% c("transcript", "tss")) {
    mcols(object)[[paste0(command@colName, 'Transcript')]] = gr$tx_name[object_nearest$subjectHits]
  }
  return(object)
}

#'
#' NearestGenCommand implementation of execute
#'
#' NearestGenCommand labels each input genomic region with information regarding to the nearest
#' gene. Gene information is obtained from the TxDb.Hsapiens.UCSC.hg19.knownGene transcripts
#' database. A gene region is defined as the union of all its transcript regions. The gene symbol
#' is obtained from the org.Hs.eg.db package.
#'
#' @importFrom AnnotationDbi mget
#' @importFrom TxDb.Hsapiens.UCSC.hg19.knownGene TxDb.Hsapiens.UCSC.hg19.knownGene
#' @importFrom org.Hs.eg.db org.Hs.egSYMBOL
#' @param command A NearestGenCommand.
#' @param object A GRanges object containing the genomic regions to annotate.
#'
setMethod('execute', c('NearestGenCommand', 'GRanges'), .executeNearestGenCommand)




#'
#' AnnotationCommandList
#'
#' This AnnotationCommand contains a list of several AnnotationCommands, and
#' executes them in order
#'
#' @export
#' @slot commandList A list containing AnnotationCommand objects.
#'
setClass('AnnotationCommandList',
         representation(commandList='list'),
         prototype(commandList=list()),
         contains='AnnotationCommand',
         validity=function(object) {
           if (object@colName != '') {
             stop('colName slot must be the empty string.')
           } else if (!is.list(object@commandList)) {
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
#' AnnotationCommand commands
#'
#' @export
#' @param ... AnnotationCommand objects in the order to be processed.
#'
annotationCommandList <- function(...) {
  commandList <- list(...)
  return(new('AnnotationCommandList', colName='', commandList=commandList))
}

#'
#' getCommandList generic
#'
#' Generic definition of the getCommandList getter methods.
#'
#' @param object An object that probably contains a commandList slot.
#'
setGeneric('getCommandList', function(object) standardGeneric('getCommandList'))

#'
#' AnnotationCommandList base implementation of getCommandList
#'
#' Getter method for commandList slot for all AnnotationCommandList objects.
#'
#' @param object An AnnotationCommandList object.
#'
setMethod('getCommandList', 'AnnotationCommandList',
          function(object) {
            return(object@commandList)
          })

#
# Internal AnnotationCommandList implementation of execute
#
.executeAnnotationCommandList <- function(command, object) {
  object <- callNextMethod()
  newRanges <- Reduce(function(xx, yy) execute(yy, xx), command@commandList, object)
  return(newRanges)
}

#'
#' AnnotationCommandList implementation of execute
#'
#' An AnnotationCommandList is just a container for several AnnotationCommand objects. When
#' executed, the objects in the internal list are executed in the same order they were provided to
#' the constructor. The colName slot of AnnotationCommandList has to be empty, because it has no
#' useful meaning at the moment.
#'
#' @param command An AnnotationCommandList.
#' @param object A GRanges object containing the genomic regions to annotate.
#'
setMethod('execute', c('AnnotationCommandList', 'GRanges'), .executeAnnotationCommandList)

