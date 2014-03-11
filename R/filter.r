#'
#' FilterCommandIndices
#'
#' Class for storing logical indices over rows and columns at the same time. Used by FilterCommand
#' objects to represent their output.
#'
setClass('FilterCommandIndices',
         representation(rows='logical', cols='logical')
         )

#' 
#' getRows generic
#'
#' Generic definition of the getRows getter methods
#'
#' @param object An object that contains a rows slot.
#'
setGeneric('getRows', function(object) standardGeneric('getRows'))

#'
#' FilterCommandIndices implementation of getRows
#'
#' Getter method for rows slot for all FilterCommandIndices objects.
#'
#' @param object A FilterCommandIndices object.
#'
setMethod('getRows', 'FilterCommandIndices',
          function(object) {
            return(object@rows)
          })

#' 
#' setRows generic
#'
#' Generic definition of the setRows setter methods
#'
#' @param object An object that contains a rows slot.
#' @param value the new rows value
#'
setGeneric('setRows<-', function(object, value) standardGeneric('setRows<-'))

#'
#' FilterCommandIndices implementation of setRows
#'
#' Setter method for rows slot for all FilterCommandIndices objects.
#'
#' @param object A FilterCommandIndices object.
#' @param value the new rows value
#' @name setRows
#'
setReplaceMethod('setRows', 'FilterCommandIndices',
          function(object, value) {
            object@rows <- value
            validObject(object)
            return(object)
          })

#' 
#' getCols generic
#'
#' Generic definition of the getCols getter methods
#'
#' @param object An object that contains a cols slot.
#'
setGeneric('getCols', function(object) standardGeneric('getCols'))

#'
#' FilterCommandIndices implementation of getCols
#'
#' Getter method for cols slot for all FilterCommandIndices objects.
#'
#' @param object A FilterCommandIndices object.
#'
setMethod('getCols', 'FilterCommandIndices',
          function(object) {
            return(object@cols)
          })

#' 
#' setCols generic
#'
#' Generic definition of the setCols setter methods
#'
#' @param object An object that contains a cols slot.
#' @param value the new cols value
#'
setGeneric('setCols<-', function(object, value) standardGeneric('setCols<-'))

#'
#' FilterCommandIndices implementation of setCols
#'
#' Setter method for cols slot for all FilterCommandIndices objects.
#'
#' @param object A FilterCommandIndices object.
#' @param value the new cols value
#' @name setCols
#'
setReplaceMethod('setCols', 'FilterCommandIndices',
          function(object, value) {
            object@cols <- value
            validObject(object)
            return(object)
          })

#'
#' FilterCommandIndices constructor
#'
#' This function builds a FilterCommandIndices object from parameters rows and cols.
#'
#' @param rows A logical vector indicating rows.
#' @param cols A logical vector indicating cols.
#' @export
#' 
filterCommandIndices <- function(rows, cols) {
  return(new('FilterCommandIndices', rows=rows, cols=cols))
}

#'
#' FilterCommand
#'
#' Abstract class that serves as interface for all the filtering commands.
#'
setClass('FilterCommand')

#'
#' FilterCommand base implementation of execute
#'
#' Implementation of the logic that is common to all FilterCommands. For 
#' example, error management and input control. For now, it only prevents the command from being
#' executed on an empty object.
#'
#' @param command A FilterCommand object.
#' @param object Any object for filtering.
#'
setMethod('execute', c('FilterCommand', 'ANY'),
          function(command, object) {
            if (is.null(dim(object))) {
              stop('Cannot execute filter command on object without dimensions.')
            } else if(is.null(object) || nrow(object) == 0 || ncol(object) == 0) {
              stop('Cannot execute filter command on empty object.')
            } else {
              return(object)
            }
          })

#'
#' AtomicFilterCommand
#'
#' Abstract class that serves as interface for all the filtering commands that are not lists. Allows
#' to define the filtering direction.
#'
#' @slot byRow Indicates the direction of filtering, wheter rows or columns should be removed.
#'
setClass('AtomicFilterCommand',
         representation(byRow='logical'),
         prototype(byRow=TRUE),
         contains='FilterCommand',
         validity=function(object) {
           return(length(object@byRow) == 1)
         }
         )

#'
#' getByRow generic
#'
#' Generic definition of the getByRow getter methods.
#'
#' @param object An object that contains a byRow slot.
#'
setGeneric('getByRow', function(object) standardGeneric('getByRow'))

#'
#' AtomicFilterCommand base implementation of getByRow
#'
#' Getter method for byRow slot for all AtomicFilterCommand objects.
#'
#' @param object An AtomicFilterCommand object.
#'
setMethod('getByRow', 'AtomicFilterCommand',
          function(object) {
            return(object@byRow)
          })

#' 
#' setByRow generic
#'
#' Generic definition of the setByRow setter methods.
#'
#' @param object An object that contains a byRow slot.
#' @param value The new byRow value.
#'
setGeneric('setByRow<-', function(object, value) standardGeneric('setByRow<-'))

#' 
#' AtomicFilterCommand base implementation of setByRow
#'
#' Getter method for byRow slot for all AtomicFilterCommand objects.
#'
#' @param object An AtomicFilterCommand object.
#' @param value The new byRow value.
#' @name setByRow
#'
setReplaceMethod('setByRow', 'AtomicFilterCommand',
          function(object, value) {
            object@byRow <- value
            validObject(object)
            return(object)
          })

#'
#' DetPFilterCommand
#'
#' Abstract class that serves as interface for all the filtering commands that use information from
#' an additional matrix. In our case, this is usually a matrix of detection p-values, and it helps
#' us to clean a set according to criteria defined over the input objec and the detection p-values
#' matrix.
#'
#' @slot detectionP The additional matrix used in the filtering criterium.
#'
setClass('DetPFilterCommand',
         representation(detectionP='matrix'),
         prototype(detectionP=NULL),
         contains='AtomicFilterCommand'
         )

#'
#' getDetectionP generic
#'
#' Generic definition of the getDetectionP getter methods.
#'
#' @param object An object that contains a detectionP slot.
#'
setGeneric('getDetectionP', function(object) standardGeneric('getDetectionP'))

#'
#' DetPFilterCommand base implementation of getDetectionP
#'
#' Getter method for detectionP slot for all DetPFilterCommand objects.
#'
#' @param object An DetPFilterCommand object.
#'
setMethod('getDetectionP', 'DetPFilterCommand',
          function(object) {
            return(object@detectionP)
          })

#' 
#' setDetectionP generic
#'
#' Generic definition of the setDetectionP setter methods.
#'
#' @param object An object that contains a detectionP slot.
#' @param value The new detectionP value.
#'
setGeneric('setDetectionP<-', function(object, value) standardGeneric('setDetectionP<-'))

#' 
#' DetPFilterCommand base implementation of setDetectionP
#'
#' Getter method for detectionP slot for all DetPFilterCommand objects.
#'
#' @param object An DetPFilterCommand object.
#' @param value The new detectionP value.
#' @name setDetectionP
#'
setReplaceMethod('setDetectionP', 'DetPFilterCommand',
          function(object, value) {
            object@detectionP <- value
            validObject(object)
            return(object)
          })

#'
#' DetPFilterCommand implementation of execute for ANY
#'
#' This base implementation of execute just checks if the dimensions of the object to be filtered
#' are compatible with those from the detection p-value matrix. 
#'
#' Dimensions do not have to be the same, but the row and column names of the input object should be 
#' included into row and column names of the detectionP matrix.
#'
#' @param command A KOverADetPFilterCommand command.
#' @param object An eSet object.
#'
setMethod('execute', c('DetPFilterCommand', 'ANY'),
          function(command, object) {
            object <- callNextMethod()

            if (!all(rownames(object) %in% rownames(command@detectionP)) ||
                !all(colnames(object) %in% colnames(command@detectionP))) {
              stop("Dimension names from input object must be included in dimension names of \
                   detection p-values matrix.")
            } else {
              return(object)
            }
          })

#'
#' KOverADetPFilterCommand filter command.
#'
#' This FilterCommand filters out the rows/columns of the input object when there are at least K
#' elements over a value A in a given row/column. It is equivalent to the kOverA filter function in
#' the package genefilter.
#'
#' @export
#' 
setClass('KOverADetPFilterCommand',
         representation(k='numeric', a='numeric'),
         prototype(k=2, a=0.01),
         contains='DetPFilterCommand',
         validity=function(object) {
           return(length(object@k) == 1 && 
                  length(object@a) == 1 && 
                  object@k > 0 &&
                  (!object@byRow || object@k <= ncol(object@detectionP)) &&
                  (object@byRow || object@k <= nrow(object@detectionP)) &&
                  object@a >= 0 &&
                  object@a <= 1
                  )
         })

#'
#' getK generic
#'
#' Generic definition of the getK getter methods.
#'
#' @param object An object that contains a k slot.
#'
setGeneric('getK', function(object) standardGeneric('getK'))

#'
#' KOverADetPFilterCommand base implementation of getK
#'
#' Getter method for k slot for all KOverADetPFilterCommand objects.
#'
#' @param object An KOverADetPFilterCommand object.
#'
setMethod('getK', 'KOverADetPFilterCommand',
          function(object) {
            return(object@k)
          })

#' 
#' setK generic
#'
#' Generic definition of the setK setter methods.
#'
#' @param object An object that contains a k slot.
#' @param value The new k value.
#'
setGeneric('setK<-', function(object, value) standardGeneric('setK<-'))

#' 
#' KOverADetPFilterCommand base implementation of setK
#'
#' Getter method for k slot for all KOverADetPFilterCommand objects.
#'
#' @param object An KOverADetPFilterCommand object.
#' @param value The new k value.
#' @name setK
#'
setReplaceMethod('setK', 'KOverADetPFilterCommand',
          function(object, value) {
            object@k <- value
            validObject(object)
            return(object)
          })

#'
#' getA generic
#'
#' Generic definition of the getA getter methods.
#'
#' @param object An object that contains an a slot.
#'
setGeneric('getA', function(object) standardGeneric('getA'))

#'
#' KOverADetPFilterCommand base implementation of getA
#'
#' Getter method for a slot for all KOverADetPFilterCommand objects.
#'
#' @param object A KOverADetPFilterCommand object.
#'
setMethod('getA', 'KOverADetPFilterCommand',
          function(object) {
            return(object@a)
          })

#' 
#' setA generic
#'
#' Generic definition of the setA setter methods.
#'
#' @param object An object that contains a a slot.
#' @param value The new a value.
#'
setGeneric('setA<-', function(object, value) standardGeneric('setA<-'))

#' 
#' KOverADetPFilterCommand base implementation of setA
#'
#' Getter method for a slot for all KOverADetPFilterCommand objects.
#'
#' @param object A KOverADetPFilterCommand object.
#' @param value The new a value.
#' @name setA
#'
setReplaceMethod('setA', 'KOverADetPFilterCommand',
          function(object, value) {
            object@a <- value
            validObject(object)
            return(object)
          })

#'
#' KOverADetPFilterCommand constructor
#'
#' This function builds a KOverADetPFilterCommand from parameters k and a.
#'
#' @param detectionP A matrix of detection p-values used for filtering.
#' @param byRow A logical indicating the direction of filtering.
#' @param k The minimum number of elements to label a row/column as invalid.
#' @param a The p-value threshold.
#' @export
#' 
kOverADetPFilterCommand <- function(detectionP, byRow=TRUE, k=ifelse(byRow, 5, 5000), a=0.01) {
  return(new('KOverADetPFilterCommand', detectionP=detectionP, k=k, a=a, byRow=byRow))
}

#'
#' KOverADetPFilterCommand percentage constructor
#'
#' This function builds a KOverADetPFilterCommand from parameter a and a percentage of rows/columns.
#'
#' @param detectionP A matrix of detection p-values used for filtering.
#' @param byRow A logical indicating the direction of filtering.
#' @param fraction The minimum percentage of elements to label a row/column as invalid.
#' @param a The p-value threshold.
#' @export
#' 
kOverADetPFilterCommandFromFraction <- function(detectionP, byRow=TRUE, fraction=0.1, a=0.01) {
  kFromFraction <- ceiling(ifelse(byRow, ncol(detectionP), nrow(detectionP)) * fraction)
  return(new('KOverADetPFilterCommand', detectionP=detectionP, k=kFromFraction, a=a, byRow=byRow))
}

#'
#' KOverADetPFilterCommand implementation of execute for eSet
#'
#' KOverADetPFilterCommand discards rows or columns from the input object according to the 
#' internal detection p-values matrix of the command. It is equivalent to the kOverA function in the
#' genefilter package. A row/column is discarded if k or more elements have a detection p-value over
#' the parameter a. Direction of filtering is controlled by the byRow parameter. 
#'
#' In order to properly chain a list of commands (see FilterCommandList), the detection p-value
#' matrix is projected according to the row and column names of the object which is being filtered.
#' 
#' @param command A KOverADetPFilterCommand command.
#' @param object An eSet object.
#'
setMethod('execute', c('KOverADetPFilterCommand', 'eSet'),
          function(command, object) {
            object <- callNextMethod()

            detectionPView <- command@detectionP[rownames(object), colnames(object)]

            if (command@byRow) {
              badSums <- rowSums(detectionPView > command@a)
            } else {
              badSums <- colSums(detectionPView > command@a)
            }

            badElements <- badSums >= command@k

            if (command@byRow) {
              return(object[!badElements, ])
            } else {
              return(object[, !badElements])
            }
          })

#'
#' SimpleFilterCommand
#'
#' Abstract class that serves as interface for all the simple filtering commands that only use 
#' information from the object to be filtered.
#'
setClass('SimpleFilterCommand', contains='AtomicFilterCommand')

#'
#' FilterCommandList
#'
#' This FilterCommand contains a list of several AtomicFilterCommands and executes them in order.
#'
#' @export
#' @slot commandList A list containing AtomicFilterCommand objects.
#'
setClass('FilterCommandList',
         representation(commandList='list'),
         prototype(commandList=list()),
         contains='FilterCommand',
         validity=function(object) {
           if (!is.list(object@commandList)) {
             stop('commandList slot must be a list.')
           } else if (!all(sapply(object@commandList, 
                                  function(xx) is(xx, 'AtomicFilterCommand')))) {
             stop('All of commandList members must be AtomicFilterCommand objects.')
           } else {
             return(TRUE)
           }
         }
         )

#'
#' FilterCommandList constructor.
#'
#' This function builds a FilterCommandList from a list of AtomicFilterCommands.
#'
#' @export
#' @param ... AtomicFilterCommand objects in the order to be processed.
#'
filterCommandList <- function(...) {
  commandList <- list(...)
  return(new('FilterCommandList', commandList=commandList))
}

#'
#' FilterCommandList base implementation of getCommandList generic.
#'
#' Getter method for commandList slot for all FilterCommandList objects.
#'
#' @param object A FilterCommandList object.
#'
setMethod('getCommandList', 'FilterCommandList',
          function(object) {
            return(object@commandList)
          })

#
# FilterCommandList internal implementation
#
.executeFilterCommandList <- function(command, object) {
  object <- callNextMethod()
  newObject <- Reduce(function(xx, yy) execute(yy, xx), command@commandList, object)
  return(newObject)
}

#'
#' FilterCommandList implementation of execute
#'
#' A FilterCommandList is just a container for several AtomicFilterCommand objects. When 
#' executed, the objects in the internal list are executed in the same order they were provided to
#' the constructor.
#'
#' @param command A FilterCommandList.
#' @param object An object to be filtered.
#'
setMethod('execute', c('FilterCommandList', 'ANY'), .executeFilterCommandList)


