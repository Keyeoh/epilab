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
#' @return A data.frame containing the results from all the tests.
#'
#' @export
#'
rangeTestAgainstList <- function(selectedRange, backgroundRange, rangeList) {
  
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
                            rangeList=rangeList)

  results <- do.call(rbind, lapply(listOfResults, 
                                   function(xx) as.data.frame(xx[setdiff(names(xx), 'Id')],
                                                              row.names=xx$Id)))

  return(results)
}
