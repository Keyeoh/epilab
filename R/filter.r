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
#' @slot byRow 
#'
setClass('AtomicFilterCommand',
         representation(byRow='logical'),
         prototype(byRow=TRUE),
         contains='FilterCommand',
         validity=function(object) {
           return(length(object@byRow) == 1)
         }
         )

