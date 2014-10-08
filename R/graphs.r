
#
# Internal helper function
#
# This function binds a name with the contents of a list element with the same name to create a
# data.frame.
#
.bindNameData <- function(name, aList) {
  values <- aList[[name]]
  if (is.factor(values)) {
      values <- as.character(values)
  }
  data.frame(name, values)
}

#'
#' Bar Group Graph
#'
#' This low-level function builds a basic bar graph from an input list. Each element of the list
#' corresponds to a bar in the graph, which is then segmented according to the factor defined by the
#' element. 
#' 
#' @param dataList A list containing the bar groups. Each of the groups has to be of type character 
#' or factor. All the groups are then concatenated and converted to factors before generating the 
#' graph.
#'
#' @importFrom ggplot2 ggplot aes_string geom_bar theme_bw scale_fill_grey
#' @export
#'
barGroupGraph <- function(dataList) {
  if (!is.list(dataList)) {
    stop('Input must be a list.')
  }

  elementIsCharacter <- sapply(dataList, is.character)
  elementIsFactor <- sapply(dataList, is.factor)

  if (!all(elementIsCharacter | elementIsFactor)) {
    stop('List elements must be of type character or factor.')
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
#' @export
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

#'
#' Violin Group Graph
#'
#' This function builds a basic violin plot from a data list. Each element of the list will be
#' represented by a different violin and named after its list name. The elements of the list have to
#' be numerical vectors.
#'
#' @param dataList A list containing the violin groups. Each of the groups has to be of type 
#' numeric. 
#'
#' @importFrom ggplot2 ggplot aes_string geom_violin theme_bw scale_fill_grey
#' @export
#'
violinGroupGraph <- function(dataList) {
  if (!is.list(dataList)) {
    stop('Input must be a list.')
  }

  elementIsNumeric <- sapply(dataList, is.numeric)

  if (!all(elementIsNumeric)) {
    stop('List elements must be of type numeric.')
  }

  plotData <- as.data.frame(Reduce(rbind, 
                                   lapply(names(dataList), .bindNameData, aList=dataList)))
  names(plotData) <- c('group', 'value')

  return(ggplot(plotData, aes_string(x='group', y='value', fill='group')) +
         geom_violin() +
         geom_boxplot(width=0.2) +
         theme_bw() +
         scale_fill_grey())
}

#'
#' Simple Bar Graph
#'
#' This function builds a simple count bar plot from a named numeric vector. Each element of the
#' vector represents the number of elements in a group (bar). The group is named after the name of
#' the vector element.
#' 
#' @param dataVector A named vector containing the counts for each group.
#'
#' @importFrom ggplot2 ggplot aes_string geom_bar theme_bw scale_fill_grey
#' @export
#'
simpleBarGraph <- function(dataVector) {

  if (!is.numeric(dataVector) || is.null(names(dataVector))) {
    stop('Input data must be a named numeric vector.')
  }

  plotData <- data.frame(group=names(dataVector), count=dataVector)

  return(ggplot(plotData, aes_string(x='group', y='count', fill='group')) +
         geom_bar(width=0.6, stat='identity', col='black') +
         theme_bw() + 
         scale_fill_grey())
}

