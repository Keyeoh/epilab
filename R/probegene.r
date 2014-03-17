#'
#' Get gene identifiers for a given set of target Ids
#'
#' This function takes as input a given set of Illumina 450k target IDs and produces a list of gene
#' Entrez IDs, according to some parameters governing the getProbeGeneRelationship inner function.
#' 
#' @param tids A character vector containing the input Target IDs.
#' @param ... Parameters to be passed to getProbeGeneRelationship.
#' @return A character vector containing the gene identifier.
#' @export
#'
getGeneEntrezIds <- function(tids, ...) {
  probeGene <- getProbeGeneRelationship(tids, ...)

  geneEntrezIds <- unique(probeGene$geneId)

  return(geneEntrezIds)
}

#'
#' Get gene symbols for a given set of Entrez IDs
#'
#' This function takes as input a set of Entrez IDs and generates a list of gene symbols using the
#' annotation information from the org.Hs.eg.db Bioconductor package.
#' 
#' @param entrezIds A character vector containing the Entrez IDs.
#' @export
#' @importFrom org.Hs.eg.db org.Hs.egSYMBOL
#'
getSymbolsFromEntrezIds <- function(entrezIds) {
  symbolList <- mget(entrezIds, org.Hs.egSYMBOL, ifnotfound=NA)
  if (length(symbolList) != 0) {
    return(unlist(symbolList))
  } else {
    return(character(0))
  }
}

#'
#' Get relationship between probes and genes
#'
#' This function computes the overlapping between Illumina450k probes and a set of genes determined
#' by a given method. The 'illumina' method just uses the annotation information in the 
#' IlluminaHumanMethylation450k.db package, while the 'ucsc19' method uses the information from 
#' UCSC hg19 genome version and computes the overlaps at a transcript level. There is also the
#' possibility to define a given promoter size to also compute the overlapping with these special
#' regions. The function uses the getSymbolsFromEntrezIds function from this package in order to 
#' get the genes' names. This function returns a data.frame with all the overlaps corresponding
#' to the input target Id's.
#' 
#' @param tids A character vector containing the input target ID's.
#' @param method Gene overlap computing method.
#' @param promoterSize Size to add to transcripts' upstream promoter regions.
#' @return A data.frame containing four columns: targetId, txId, geneId and geneSymbol.
#' @export
#' @importFrom IlluminaHumanMethylation450k.db IlluminaHumanMethylation450kENTREZID
#' @importFrom FDb.InfiniumMethylation.hg19 get450k
#' @importFrom TxDb.Hsapiens.UCSC.hg19.knownGene TxDb.Hsapiens.UCSC.hg19.knownGene
#'
getProbeGeneRelationship <- function(tids, method=c('illumina', 'ucsc19'), promoterSize=2000) {
  method <- match.arg(method)

  if (!is.character(tids)) {
    stop('Target IDs should be defined as a character vector.')
  }

  if (!is.numeric(promoterSize) || length(promoterSize) != 1 || promoterSize < 0) {
    stop('Promoter size should be a numeric atomic value.')
  }

  if (method == 'illumina') {
    .Deprecated(msg='The Illumina Annotation package has been marked as deprecated in BioC')

    multipleGeneIds <- toggleProbes(IlluminaHumanMethylation450kENTREZID, 'all')
    idTable <- toTable(multipleGeneIds)
    subtableWithTids <- idTable[idTable$probe_id %in% tids, ]

    probeGene <- data.frame(probeId=subtableWithTids$probe_id,
                            geneId=subtableWithTids$gene_id,
                            stringsAsFactors=FALSE)
    probeGene$geneSymbol <- getSymbolsFromEntrezIds(probeGene$geneId)
  } else if (method == 'ucsc19') {
    transcriptsByGene <- transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, 'gene')
    transcriptRanges <- unlist(transcriptsByGene)
    
    ranges450k <- get450k()
    transcriptPromoterRanges <- resize(transcriptRanges, 
                                       width(transcriptRanges) + promoterSize,
                                       fix='end')
    overlaps <- findOverlaps(ranges450k[tids], transcriptPromoterRanges)

    probeGene <- data.frame(probeId=tids[queryHits(overlaps)],
                            geneId=names(transcriptPromoterRanges)[subjectHits(overlaps)],
                            #txId=transcriptPromoterRanges$tx_name[subjectHits(overlaps)],
                            stringsAsFactors=FALSE)
    probeGene$geneSymbol <- getSymbolsFromEntrezIds(probeGene$geneId)
  } else {
    stop('Method not supported.')
  }

  probeGene <- unique(probeGene[order(probeGene$probeId), ])
  rownames(probeGene) <- NULL

  return(probeGene)
}

