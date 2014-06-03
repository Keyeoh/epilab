#'
#' WhichOne functional
#'
#' This function just returns a given element of a vector.
#'
#' @param x An indexable object.
#' @param i Integer value indicating a position.
#' @return The i'th element of x.
#'
whichOne <- function(x, i) {
  if (i <= 0) {
    stop('Index must be positive.')
  } else {
    return(x[i])
  }
}

#'
#' First closure
#'
#' Just a closure for returning the first element of a vector.
#'
#' @param x An indexable object.
#' @return The first element of x.
#'
first <- function(x) whichOne(x, 1)

#'
#' Second closure
#'
#' Just a closure for returning the second element of a vector.
#'
#' @param x An indexable object.
#' @return The second element of x.
#'
second <- function(x) whichOne(x, 2)

#'
#' matchLength closure
#'
#' Just a closure for returning the match.length attribute of an object.
#'
#' @param x An object with a match.length attribute.
#' @return The value of the match.length attribute.
#'
matchLength <- function(x) {
  mLength <- attr(x, 'match.length')
  if (is.null(mLength)) {
    stop('Input object must have a match.length attribute.')
  } else {
    return(mLength)
  }
}

#'
#' secondLength closure
#'
#' Just a closure for returning the second element of a match.length attribute of an object.
#'
#' @param x An object with a match.length attribute.
#' @return The value of the second element of the match.length attribute.
#'
secondLength <- function(x) second(matchLength(x))

#'
#' Get all values for a key
#'
#' This function retrieves all values associated to a given key from a character vector representing
#' a text file. This text file is encoded as an identifier and a variable set of keys per line. The
#' keys are in the form key=value; while the identifier is the first element in the line and it is 
#' separated from the rest of elements with a tab character.
#'
#' @param key A character value representing the value we are looking for.
#' @param txts A character vector containing the lines of a text file.
#' @return A character vector with all the values corresponding to the input key.
#'
getValuesForKey <- function(key, txts) {
  if (!is.character(key) || !length(key) == 1) {
    stop('Key must be a length 1 character.')
  }
  matches <- regexec(paste0(key, '=([^;]+);'), txts)
  positions <- sapply(matches, second)
  matchLengths <- sapply(matches, secondLength)
  values <- substring(txts, positions, positions + matchLengths - 1)
  return(values)
}

#'
#' Create initial data structure
#'
#' This function reads the identifiers from a character vector representing a file and populates
#' an initial data.frame with their values.
#'
#' @param txtLines A character vector containing the lines of a text file.
#' @return A data.frame with a single column containing all the identifiers.
#'
createFilesData <- function(txtLines) {
  stemMatches <- regexec('([a-zA-Z0-9]+).broadPeak.gz\t', txtLines)
  stemPositions <- sapply(stemMatches, second)
  stemLengths <- sapply(stemMatches, secondLength)
  stems <- substring(txtLines, stemPositions, stemPositions + stemLengths - 1)
  return(data.frame(stem=stems, stringsAsFactors=FALSE))
}

#'
#' Get all keys from file
#'
#' This function retrieves all the different keys that are present in the input representation. As 
#' the number of keys and their names may be variable, it is necessary first to know which keys are 
#' present in the input in order to search for them afterwards.
#'
#' @param txtLines A character vector containing the lines of a text file.
#' @return A character vector containing all the keys found in the input.
#'
getAllKeys <- function(txtLines) {
  if (!is.character(txtLines)) {
    stop('Input text lines must be a character vector.')
  }
  keyPositions <- gregexpr('[a-zA-Z0-9]+=', txtLines)
  keyLengths <- lapply(keyPositions, matchLength)
  keyList <- lapply(1:length(keyPositions), 
                    function(xx) substring(txtLines[xx], keyPositions[[xx]], 
                                           keyPositions[[xx]] + keyLengths[[xx]] - 2))
  return(allKeys <- Reduce(union, keyList))
}

#'
#' Add values from key to data.frame
#'
#' This function adds a column to an existing data.frame. It names the column after the input key 
#' and then populates it with all the values it can find in the input file.
#'
#' @param aDataFrame The input data.frame.
#' @param aKey The key it is going to search for.
#' @param refLines A character vector representation of the input file.
#' @return A new data.frame with the new column added.
#'
addKey <- function(aDataFrame, aKey, refLines) {
  aDataFrame[[aKey]] <- getValuesForKey(aKey, refLines)
  return(aDataFrame)
}

#'
#' Add values closure
#'
#' This is a factory that returns a function that adds all key values to a data.frame, but for a
#' given and fixed character vector representing an input file.
#'
#' @param refLines A character vector representation of the input file.
#' @return A function for adding keys to a data.frame.
#'
addKeyFactory <- function(refLines) {
  if (!is.character(refLines)) {
    stop('Input must be a character vector.')
  }
  return(function(x, y) addKey(x, y, refLines=refLines))
}

#' 
#' Add all keys from a file
#'
#' This function adds all the keys found in a file to a given data.frame.
#'
#' @param initData The input data.frame.
#' @param keys A character vector with all the keys it is going to search for.
#' @param txtLines A character vector representation of the input file.
#' @return A new data.frame with the new columns added.
#'
addAllKeys <- function(initData, keys, txtLines) {
  return(Reduce(addKeyFactory(txtLines), keys, init=initData))
}

#'
#' Read a variable key, value file
#'
#' The files describing the ENCODE Broad Histone experiments are in a special format, where for each
#' line, there is an identifier followed by a variable set of key, value pairs encoding the 
#' necessary information. This function reads the information in a file and returns a data.frame
#' with a column for each key.
#'
#' @param txtLines A character vector representation of the input file.
#' @return A new data.frame with the new columns added.
#'
#' @export
#'
readBroadHistoneFile <- function(txtLines) {
  filesData <- createFilesData(txtLines)
  allKeys <- getAllKeys(txtLines)
  filesData <- addAllKeys(filesData, allKeys, txtLines)
  return(filesData)
}


