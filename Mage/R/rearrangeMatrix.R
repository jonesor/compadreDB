#' Re-arrange stages from prop, preRep, rep and postRep
#'
#' @export
#' @param matU a matrix without reproduction and clonality
#' @param matF a matrix, just reproduction
#' @param matFmu a matrix, that has been rearranged
#' @examples
#' matU <- matrix(c(0.2581, 0.1613, 0.1935, 0.2258, 0.1613, 0.0408, 0.2857,
#'                  0.4286, 0.102, 0.0816, 0.0385, 0.0385, 0.2692, 0.2308,
#'                  0.3462, 0, 0.0625, 0.125, 0.25, 0.5625, 0.1061, 0.1608,
#'                  0.2637, 0.1801, 0.2058),
#'                  nrow = 5, byrow = FALSE)
#' matF <- matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#'                  0, 0, 0, 0, 2.75, 1.75, 0, 0),
#'                  nrow = 5, byrow = FALSE)
#' matFmu <- matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#'                    0, 0, 0, 0, 1.0409262, 0.5040727, 0.016433, 0.06956801),
#'                  nrow = 5, byrow = FALSE)
#' rearrangeMatrix(matU, matF, matFmu)
rearrangeMatrix <- function(matU, matF, matFmu) {
  reArrange <- NULL
  matDim <- dim(matF)[1]
  Rep <- which(colSums(matFmu) > 0)
  #These are stages that are inter-reproductive but are truly non-reproductive
  nonRepInterRep <- Rep[1] - 1 +
    as.numeric(c(which(colSums(as.matrix(matFmu[,Rep[1]:Rep[length(Rep)]])) == 0)))
  if (length(nonRepInterRep) > 0) {
    allElseStages <- 1:matDim
    allElseStages <- allElseStages[-which(allElseStages %in% nonRepInterRep)]
    reArrangeStages <- c(allElseStages, nonRepInterRep)
    reArrangeMatU <- matU[reArrangeStages, reArrangeStages]
    reArrangeMatF <- matF[reArrangeStages, reArrangeStages]
    reArrangeMatFmu <- matFmu[reArrangeStages, reArrangeStages]
    reArrange$matU <- reArrangeMatU
    reArrange$matF <- reArrangeMatF
    reArrange$matFmu <- reArrangeMatFmu
  }
  if (length(nonRepInterRep) == 0) {
    reArrange$matU <- matU
    reArrange$matF <- matF
    reArrange$matFmu <- matFmu
  }
  reArrange$nonRepInterRep <- nonRepInterRep
  reArrange$maxRep <- max(Rep)

  return(reArrange)
}
