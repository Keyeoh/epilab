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
            if(is.null(object) || nrow(object) == 0 || ncol(object) == 0) {
              stop('Cannot execute annotation command on empty object.')
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


