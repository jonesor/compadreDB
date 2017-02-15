#' Re-arrange matrix stages so that all inter-reproductive and non-reproductive
#' stages fall in the final rows/columns of the matrix. This is a preparatory
#' step to collapsing the matrix model into combined prop, pre-rep, rep and
#' post-rep stages.
#'
#' @export
#' @param matU survival matrix
#' @param matF fecundity matrix
#' @param matFmu mean fecundity matrix
#' @examples
#' ## FIXME: find an example containing non-reproductive stages
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
  if (!(identical(dim(matU), dim(matF)) && identical(dim(matF), dim(matFmu)))) {
    stop("Expecting matrices with equal dimensions", call. = FALSE)
  }
  reArrange <- NULL
  matDim <- dim(matF)[1]
  Rep <- which(colSums(matFmu) > 0)
  allRep <- Rep[1]:Rep[length(Rep)]
  ## These are stages that are inter-reproductive but are truly non-reproductive:
  nonRepInterRep <- allRep[which(!allRep %in% Rep)]
  if (length(nonRepInterRep) > 0) {
    allElseStages <- which(!1:matDim %in% nonRepInterRep)
    reArrangeStages <- c(allElseStages, nonRepInterRep)
    reArrange$matU <- matU[reArrangeStages, reArrangeStages]
    reArrange$matF <- matF[reArrangeStages, reArrangeStages]
    reArrange$matFmu <- matFmu[reArrangeStages, reArrangeStages]
  } else {
    ## No non-repro or inter-repro stages so no need to rearrange matrices
    reArrange$matU <- matU
    reArrange$matF <- matF
    reArrange$matFmu <- matFmu
  }
  reArrange$nonRepInterRep <- nonRepInterRep
  reArrange$maxRep <- max(Rep)

  return(reArrange)
}
