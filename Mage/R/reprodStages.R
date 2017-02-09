#' Reproductive stages
#'
#' Determine pre-reproductive, reproductive and post-reproductive stages
#' from a matrix model which pre-formats them to use them as the "collapse"
#' argument of the function above to collapse matrices
#'
#' FIXME: ROB says possibly put output as a named list instead of vector
#'
#' @export
#' @param matF a matrix, just reproduction
#' @param matFmu a matrix
#' @param post post something
#' @param maxRep maximum xxx
#' @param matrixStages stages
#' @examples
#' # reprodStages(x, y, z, a)
reprodStages <- function(matF, matFmu, post, maxRep, matrixStages = NULL) {
  propStage <- NULL
  if ("prop" %in% matrixStages) {
    propStage <- which(matrixStages == "prop")
  } else {
    propStage <- NA
  }

  matDim <- dim(matF)[1]
  Rep <- which(colSums(matFmu) > 0)
  if (min(Rep) == 1) {
    preRep <- NA
  } else if (!is.na(propStage[1]) & (min(Rep) - max(propStage) == 1)) {
    preRep <- NA
  } else {
    preRep <- min(which(matrixStages == "active")):(min(Rep) - 1)
  }
  if (length(post) == 0 & maxRep == matDim) {
    postRep <- NA
  } else {
    postRep <- (maxRep + 1):matDim
  }

  if (length(propStage) > 1) {
    propStages <- paste(propStage[1],"-",propStage[length(propStage)], sep = "")
  } else{
    propStages <- as.character(propStage)
  }
  if (length(preRep) > 1) {
    preRepStages <- paste(preRep[1], "-", preRep[length(preRep)], sep = "")
  } else {
    preRepStages <- as.character(preRep)
  }
  if (length(Rep) > 1) {
    repStages <- paste(Rep[1], "-", Rep[length(Rep)], sep = "")
  } else {
      repStages <- as.character(Rep)
  }
  if (length(postRep) > 1) {
    postRepStages <- paste(postRep[1], "-", postRep[length(postRep)], sep = "")
  } else {
    postRepStages <- as.character(postRep)
  }

  stages <- c(propStages,preRepStages,repStages,postRepStages)

  return(stages)
}
