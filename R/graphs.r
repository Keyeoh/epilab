
#'
#' Bar Group Graph
#'
#' This low-level function builds a basic bar graph from an input list. Each element of the list
#' corresponds to a bar in the graph, which is then segmented according to the factor defined by the
#' element. 
#' 
#' @param dataList A list containing the bar groups. Each of the groups has to be of type vector or
#' factor. All the groups are then concatenated and converted to factors before generating the 
#' graph.
#'
#' @importFrom ggplot2 ggplot aes_string geom_bar theme_bw scale_fill_grey
#'
barGroupGraph <- function(dataList) {
  .bindNameData <- function(name, aList) {
    cbind(name, as.character(aList[[name]]))
  }

  if (!is.list(dataList)) {
    stop('Input must be a list.')
  }

  plotData <- as.data.frame(Reduce(rbind, 
                                   lapply(names(dataList), .bindNameData, aList=dataList)))
  names(plotData) <- c('group', 'status')

  return(ggplot(plotData, aes_string(x='group', fill='status')) +
         geom_bar(position='fill', color='black', width=0.6) +
         theme_bw() +
         scale_fill_grey())
}

#'
#' Bar Group Graph from matrix
#'
#' This low-level function builds a basic bar graph from an input matrix. The columns in the matrix
#' indicate the groups/bars, and the rows correspond to the different levels of the factor being 
#' studied. It is similar to the barGroupGraph function, but aimed to be used when the data list is
#' big.
#' 
#' @param dataMatrix A matrix containing the bar group element counts. There is a colunn for each of
#' the bars/groups and a row for each level of the factor being studied.
#'
#' @importFrom ggplot2 ggplot aes_string geom_bar theme_bw scale_fill_grey
#' @importFrom reshape2 melt 
#'
barGroupGraphFromMatrix <- function(dataMatrix) {

  if (!is.matrix(dataMatrix)) {
    stop('Input must be a matrix')
  }

  plotData <- melt(dataMatrix, varnames=c('status', 'group'))

  return(ggplot(plotData, aes_string(x='group', y='value', fill='status')) +
         geom_bar(position='fill', color='black', width=0.6, stat='identity') +
         theme_bw() +
         scale_fill_grey())
}

