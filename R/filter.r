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
#' MatrixFilterCommand
#'
#' Abstract class that serves as interface for all the filtering commands that use information from
#' an additional matrix. In our case, this is usually a matrix of detection p-values, and it helps
#' us to clean a set according to criteria defined over the input objec and the detection p-values
#' matrix.
#'
#' @slot m The additional matrix used in the filtering criterium.
#'
setClass('MatrixFilterCommand',
         representation(m='matrix'),
         prototype(m=NULL),
         contains='AtomicFilterCommand'
         )

#'
#' getMatrix generic
#'
#' Generic definition of the getMatrix getter methods.
#'
#' @param object An object that contains a m slot.
#'
setGeneric('getMatrix', function(object) standardGeneric('getMatrix'))

#'
#' MatrixFilterCommand base implementation of getMatrix
#'
#' Getter method for m slot for all MatrixFilterCommand objects.
#'
#' @param object A MatrixFilterCommand object.
#'
setMethod('getMatrix', 'MatrixFilterCommand',
          function(object) {
            return(object@m)
          })

#' 
#' setMatrix generic
#'
#' Generic definition of the setMatrix setter methods.
#'
#' @param object An object that contains a m slot.
#' @param value The new m value.
#'
setGeneric('setMatrix<-', function(object, value) standardGeneric('setMatrix<-'))

#' 
#' MatrixFilterCommand base implementation of setMatrix
#'
#' Getter method for m slot for all MatrixFilterCommand objects.
#'
#' @param object A MatrixFilterCommand object.
#' @param value The new m value.
#' @name setMatrix
#'
setReplaceMethod('setMatrix', 'MatrixFilterCommand',
          function(object, value) {
            object@m <- value
            validObject(object)
            return(object)
          })

#'
#' MatrixFilterCommand implementation of execute for ANY
#'
#' This base implementation of execute just checks if the dimensions of the object to be filtered
#' are compatible with those from the additional matrix. 
#'
#' Dimensions do not have to be the same, but the row and column names of the input object should be 
#' included into row and column names of the additional matrix.
#'
#' @param command A MatrixFilterCommand command.
#' @param object An object to be filtered.
#'
setMethod('execute', c('MatrixFilterCommand', 'ANY'),
          function(command, object) {
            object <- callNextMethod()

            if (!all(rownames(object) %in% rownames(command@m)) ||
                !all(colnames(object) %in% colnames(command@m))) {
              stop("Dimension names from input object must be included in dimension names of \
                   additional matrix.")
            } else {
              return(object)
            }
          })

#'
#' KOverAFilterCommand filter command.
#'
#' This FilterCommand filters out the rows/columns of the input object when there are at least K
#' elements over a value A in a given row/column. It is equivalent to the kOverA filter function in
#' the package genefilter.
#'
#' @export
#' 
setClass('KOverAFilterCommand',
         representation(k='numeric', a='numeric'),
         prototype(k=2, a=0.01),
         contains='MatrixFilterCommand',
         validity=function(object) {
           return(length(object@k) == 1 && 
                  length(object@a) == 1 && 
                  object@k > 0 &&
                  (!object@byRow || object@k <= ncol(object@m)) &&
                  (object@byRow || object@k <= nrow(object@m)) &&
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
#' KOverAFilterCommand base implementation of getK
#'
#' Getter method for k slot for all KOverAFilterCommand objects.
#'
#' @param object An KOverAFilterCommand object.
#'
setMethod('getK', 'KOverAFilterCommand',
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
#' KOverAFilterCommand base implementation of setK
#'
#' Getter method for k slot for all KOverAFilterCommand objects.
#'
#' @param object An KOverAFilterCommand object.
#' @param value The new k value.
#' @name setK
#'
setReplaceMethod('setK', 'KOverAFilterCommand',
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
#' KOverAFilterCommand base implementation of getA
#'
#' Getter method for a slot for all KOverAFilterCommand objects.
#'
#' @param object A KOverAFilterCommand object.
#'
setMethod('getA', 'KOverAFilterCommand',
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
#' KOverAFilterCommand base implementation of setA
#'
#' Getter method for a slot for all KOverAFilterCommand objects.
#'
#' @param object A KOverAFilterCommand object.
#' @param value The new a value.
#' @name setA
#'
setReplaceMethod('setA', 'KOverAFilterCommand',
          function(object, value) {
            object@a <- value
            validObject(object)
            return(object)
          })

#'
#' KOverAFilterCommand constructor
#'
#' This function builds a KOverAFilterCommand from parameters k and a.
#'
#' @param m A matrix used for filtering.
#' @param byRow A logical indicating the direction of filtering.
#' @param k The minimum number of elements to label a row/column as invalid.
#' @param a The p-value threshold.
#' @export
#' 
kOverAFilterCommand <- function(m, byRow=TRUE, k=ifelse(byRow, 5, 5000), a=0.01) {
  return(new('KOverAFilterCommand', m=m, k=k, a=a, byRow=byRow))
}

#'
#' KOverAFilterCommand percentage constructor
#'
#' This function builds a KOverAFilterCommand from parameter a and a percentage of rows/columns.
#'
#' @param m A matrix used for filtering.
#' @param byRow A logical indicating the direction of filtering.
#' @param fraction The minimum percentage of elements to label a row/column as invalid.
#' @param a The p-value threshold.
#' @export
#' 
kOverAFilterCommandFromFraction <- function(m, byRow=TRUE, fraction=0.1, a=0.01) {
  kFromFraction <- ceiling(ifelse(byRow, ncol(m), nrow(m)) * fraction)
  return(new('KOverAFilterCommand', m=m, k=kFromFraction, a=a, byRow=byRow))
}

#
# Common implementation of KOverAFilterCommand execute
#
.executeKOverAFilterCommand <- function(command, object) {
  mView <- command@m[rownames(object), colnames(object)]

  if (command@byRow) {
    badSums <- rowSums(mView > command@a)
  } else {
    badSums <- colSums(mView > command@a)
  }

  badElements <- badSums >= command@k

  if (command@byRow) {
    return(object[!badElements, ])
  } else {
    return(object[, !badElements])
  }
}

#'
#' KOverAFilterCommand implementation of execute for ANY
#'
#' KOverAFilterCommand discards rows or columns from the input object according to the 
#' internal detection p-values matrix of the command. It is equivalent to the kOverA function in the
#' genefilter package. A row/column is discarded if k or more elements have a detection p-value over
#' the parameter a. Direction of filtering is controlled by the byRow parameter. 
#'
#' In order to properly chain a list of commands (see FilterCommandList), the detection p-value
#' matrix is projected according to the row and column names of the object which is being filtered.
#' 
#' @param command A KOverAFilterCommand command.
#' @param object An eSet object.
#'
setMethod('execute', c('KOverAFilterCommand', 'ANY'),
          function(command, object) {
            object <- callNextMethod()

            return(.executeKOverAFilterCommand(command, object))
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
# TODO: Improve the efficiency of this implementation.
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

#'
#' VarFilterCommand filter command
#'
#' This FilterCommand filters out the rows (or columns) of the input object where the variance meets
#' a certain requirement. Two criteria are implemented at the moment: filter the rows (or columns)
#' where the variance is under a given threshold, or filter the rows (or columns) where the variance
#' is under a given quantile of all elements' variances.
#'
#' @slot type Indicates type of filtering. Should be one of 'quantile' or 'absolute'
#' @slot threshold The numerical value of the threshold.
#'
setClass('VarFilterCommand',
         representation(type='character', threshold='numeric'),
         prototype(type='quantile', threshold=0.25),
         contains='MatrixFilterCommand',
         validity=function(object) {
           if (!object@type %in% c('quantile', 'absolute')) {
             return('Type of filtering can only be quantile or absolute.')
           } else if (object@type == 'quantile' && (object@threshold < 0 || object@threshold > 1)) {
             return('Quantile threshold should be in [0, 1].')
           } else if (object@type == 'absolute' && object@threshold < 0) {
             return('Absolute threshold should be a positive number.')
           }  else {
             return(TRUE)
           }
         })

#'
#' getType generic
#'
#' Generic definition of the getType getter methods.
#'
#' @param object An object that contains a type slot.
#'
setGeneric('getType', function(object) standardGeneric('getType'))

#'
#' VarFilterCommand base implementation of getType
#'
#' Getter method for type slot for all VarFilterCommand objects.
#'
#' @param object An VarFilterCommand object.
#'
setMethod('getType', 'VarFilterCommand',
          function(object) {
            return(object@type)
          })

#' 
#' setType generic
#'
#' Generic definition of the setType setter methods.
#'
#' @param object An object that contains a type slot.
#' @param value The new type value.
#'
setGeneric('setType<-', function(object, value) standardGeneric('setType<-'))

#' 
#' VarFilterCommand base implementation of setType
#'
#' Getter method for type slot for all VarFilterCommand objects.
#'
#' @param object An VarFilterCommand object.
#' @param value The new type value.
#' @name setType
#'
setReplaceMethod('setType', 'VarFilterCommand',
          function(object, value) {
            object@type <- value
            validObject(object)
            return(object)
          })

#'
#' getThreshold generic
#'
#' Generic definition of the getThreshold getter methods.
#'
#' @param object An object that contains a threshold slot.
#'
setGeneric('getThreshold', function(object) standardGeneric('getThreshold'))

#'
#' VarFilterCommand base implementation of getThreshold
#'
#' Getter method for threshold slot for all VarFilterCommand objects.
#'
#' @param object An VarFilterCommand object.
#'
setMethod('getThreshold', 'VarFilterCommand',
          function(object) {
            return(object@threshold)
          })

#' 
#' setThreshold generic
#'
#' Generic definition of the setThreshold setter methods.
#'
#' @param object An object that contains a threshold slot.
#' @param value The new threshold value.
#'
setGeneric('setThreshold<-', function(object, value) standardGeneric('setThreshold<-'))

#' 
#' VarFilterCommand base implementation of setThreshold
#'
#' Getter method for threshold slot for all VarFilterCommand objects.
#'
#' @param object An VarFilterCommand object.
#' @param value The new threshold value.
#' @name setThreshold
#'
setReplaceMethod('setThreshold', 'VarFilterCommand',
          function(object, value) {
            object@threshold <- value
            validObject(object)
            return(object)
          })

#'
#' VarFilterCommand constructor
#'
#' This function builds a VarFilterCommand from parameters type and threshold.
#'
#' @param m An additional matrix for evaluating the variance.
#' @param byRow A logical indicating the direction of filtering.
#' @param type Indicates type of filtering. Should be one of 'quantile' or 'absolute'
#' @param threshold The numerical value of the threshold.
#' @export
#' 
varFilterCommand <- function(m, byRow=TRUE, type='quantile', threshold=0.25) {
  return(new('VarFilterCommand', m=m, byRow=byRow, type=type, threshold=threshold))
}

# 
# Common implementation for VarFilterCommand execute
#
.executeVarFilterCommand <- function(command, object) {
  if (command@byRow) {
    variances <- rowVars(command@m)
  } else {
    variances <- colVars(command@m)
  }

  if (command@type == 'absolute') {
    badElements <- variances < command@threshold
  } else if (command@type == 'quantile') {
    badElements <- variances < quantile(variances, command@threshold)
  } else {
    stop('Should not arrive here!')
  }

  if (command@byRow) {
    return(object[!badElements, ])
  } else {
    return(object[, !badElements])
  }
}

#'
#' VarFilterCommand implementation of execute for ANY
#'
#' VarFilterCommand discards rows or columns according to their elements' variance. Depending on the
#' type of filtering, the threshold used can be absolute or expressed as a quantile of the whole
#' set of variances. This is a useful filter for the implementation of non-specific filtering.
#'
#' @param command A VarFilterCommand command.
#' @param object An object to be filtered.
#' @importFrom matrixStats rowVars colVars
#'
setMethod('execute', c('VarFilterCommand', 'ANY'),
          function(command, object) {
            object <- callNextMethod()

            return(.executeVarFilterCommand(command, object))
          })

