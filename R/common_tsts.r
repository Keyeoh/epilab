#'
#' Test ranges against background with respect to a list of ranges.
#'
#' This function traverses the provided GRangesList and, for each element of the list, it performs
#' a simple proportionality test based on a 2x2 contingency table. This is similar to a gene set
#' testing (GO-like), but a lower level and intended for being used as a building block for more
#' complex testing procedures. Currently, Fisher test is the only choice.
#'
#' @param selectedRange A GRanges object containing the sequences of interest.
#' @param backgroundRange A GRanges object containing sequences which are not supposed to resemble
#' the selected subset.
#' @param rangeList A GRangesList object containing the different subsets we are going to test 
#' against.
#' @param \dots Parameters to be passed down to mclapply.
#' @return A data.frame containing the results from all the tests.
#'
#' @export
#'
rangeTestAgainstList <- function(selectedRange, backgroundRange, rangeList, ...) {
 
  if (!is(selectedRange, 'GenomicRanges')) {
    stop('Selected range must be a GenomicRanges object')
  }
  if (!is(backgroundRange, 'GenomicRanges')) {
    stop('Background range must be a GenomicRanges object')
  }
  if (!is(rangeList, 'GRangesList')) {
    stop('Ranges list must be a GRangesList object')
  }
  if (length(rangeList) == 0) {
    stop('Ranges list must contain at list a GRanges object')
  }
  if (unique(lapply(rangeList, class)) != 'GRanges') {
    stop('Ranges list must contain only GRanges objects')
  }
  if (is.null(names(rangeList))) {
    stop('Objects in ranges list must have names')
  }

  testASingleElement <- function(elementName, selectedRange, backgroundRange, rangeList) {

    currentElement <- rangeList[[elementName]]
    currentElement <- reduce(currentElement)

    sigTable <- matrix(0, nrow=2, ncol=2, 
                       dimnames=list(c('mark', 'nomark'), c('selected', 'background')))
    sigTable['mark', 'selected'] <- sum(countOverlaps(selectedRange, currentElement) > 0)
    sigTable['mark', 'background'] <- sum(countOverlaps(backgroundRange, currentElement) > 0)
    sigTable['nomark', 'selected'] <- length(selectedRange) - sigTable['mark', 'selected']
    sigTable['nomark', 'background'] <- length(backgroundRange) - sigTable['mark', 'background']
 
    currentTest <- fisher.test(sigTable)

    currentResult <- list(
                          Id=elementName,
                          Mark_Selected=sigTable['mark', 'selected'],
                          Mark_Background=sigTable['mark', 'background'],
                          NoMark_Selected=sigTable['nomark', 'selected'],
                          NoMark_Background=sigTable['nomark', 'background'],
                          PValue=currentTest$p.value,
                          Method=currentTest$method
                          )
    currentResult$OR <- oddsRatio(sigTable)
    currentResult$Log2OR <- log2(currentResult$OR)
    currentResult$P_Selected <- sigTable['mark', 'selected'] / sum(sigTable[, 'selected'])
    currentResult$P_Background <- sigTable['mark', 'background'] / sum(sigTable[, 'background'])

    return(currentResult)
  }

  listOfResults <- mclapply(names(rangeList), testASingleElement, 
                            selectedRange=selectedRange,
                            backgroundRange=backgroundRange,
                            rangeList=rangeList, ...)

  results <- do.call(rbind, lapply(listOfResults, 
                                   function(xx) as.data.frame(xx[setdiff(names(xx), 'Id')],
                                                              row.names=xx$Id)))

  results$QValue <- p.adjust(results$PValue, method='fdr')
  results$Bonferroni <- p.adjust(results$PValue, method='bonferroni')

  return(results)
}

#'
#' Test Illumina450k target ids against a list of ranges.
#'
#' This function is a wrapper for the rangeTestAgainstList, specifically designed to handle
#' Illumina450k Target IDs. It internally transforms the ids into both selected and background
#' GRanges objects and then call the main workhorse function.
#'
#' @param targetIds A character vector containing the Illumina450k target identifiers.
#' @param rangeList A GRangesList object containing the different subsets we are going to test 
#' against.
#' @param \dots Parameters to be passed down to rangeTestAgainstList.
#' @return A data.frame containing the results from all the tests.
#'
#' @export
#'
tids450kTestAgainstList <- function(targetIds, rangeList, ...) {
  hm450 <- get450k()
  query <- hm450[targetIds]
  background <- hm450[setdiff(names(hm450), targetIds)]
  return(rangeTestAgainstList(query, background, rangeList, ...))
}

#'
#' General categorical testing.
#'
#' This function is a general building block for the common practice of testing if a given subset
#' of elements is differentially enriched or impovered with respect to a given factor. It builds
#' a contingency matrix, performs a statistical test and returns information about the test results
#' and effect sizes.
#'
#' @param targetFactor A factor or character vector (which is then converted to a factor) encoding 
#' the universe for the test along with its labels.
#' @param selectedIndices A numeric, logical or character vector indicating the elements of the 
#' factor that resemble the subset currently being tested.
#' @param testId A character element identifying the current test.
#' @return A list containing information about the result.
#'
#' @export
#'
categoricalTest <- function(targetFactor, selectedIndices, testId=NULL) {
 
  if (any(is.na(selectedIndices))) {
    stop('Indices must not contain NA')
  }

  if (!class(targetFactor) %in% c('character', 'numeric', 'factor')) {
    stop('Target Factor must be of character, numeric or factor type')
  }

  if (!class(selectedIndices) %in% c('numeric', 'logical', 'character', 'integer')) {
    stop('Selected indices must be of numeric, logical or character type')
  }

  if (!is.null(testId) && !(is(testId, 'character') && length(testId) == 1)) {
    stop('Test Identifier must be a character of length 1')
  }

  if (is(selectedIndices, 'numeric')) {
    if (any(selectedIndices > length(targetFactor))) {
      stop('All numeric indices must be lower or equal than the size of the factor')
    }
    inOut <- ifelse(1:length(targetFactor) %in% selectedIndices, 'In', 'Out')
  } else if (is(selectedIndices, 'logical')) {
    if (length(selectedIndices) != length(targetFactor)) {
      stop('Logical indices must be of the same length than the factor')
    }
    inOut <- ifelse(selectedIndices, 'In', 'Out')
  } else if (is(selectedIndices, 'character')) {
    if (length(setdiff(selectedIndices, names(targetFactor))) > 0) {
      stop('Names of character indices must be included in names of the target factor')
    }
    inOut <- ifelse(names(targetFactor) %in% selectedIndices, 'In', 'Out')
  } else {
    stop('Indices must be  of type numeric, logical or character')
  }
  inOut <- factor(inOut, levels=c('In', 'Out'))
  
  if (is(targetFactor, 'character') || is(targetFactor, 'numeric')) {
    targetFactor <- as.factor(targetFactor)
  }

  sigTable <- table(targetFactor, inOut)
  sigTest <- chisq.test(sigTable)

  results <- list()
  if (!is.null(testId)) {
    results$Id <- testId
  }
  for (l in levels(targetFactor)) {
    results[[paste(l, 'In', sep='_')]] <- sigTable[l, 'In']
  }
  for (l in levels(targetFactor)) {
    results[[paste(l, 'Out', sep='_')]] <- sigTable[l, 'Out']
  }
  results$PValue <- sigTest$p.value
  results$Method <- sigTest$method
  for (l in levels(targetFactor)) {
    results[[paste('OR', l, sep='_')]] <- oddsRatioLevel(sigTable, l)
  }
  for (l in levels(targetFactor)) {
    results[[paste('P', l, sep='_')]] <- sigTable[l, 'In'] / sum(sigTable[, 'In'])
  }

  return(results)
}
