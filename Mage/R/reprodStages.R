#' Group the stages of a matrix model into propagule, pre-reproductive,
#' reproductive and post-reproductive stages.
#'
#' Assumes that fecundity and mean fecundity matrices have been rearranged so
#' that non-reproductive stages are in the final rows/columns.
#'
#' Output indicates groupings to be used when collapsing the matrix model.
#'
#' Note: dormant stages not handled.
#'
#' FIXME: ROB says possibly put output as a named list instead of vector
#' FIXME: can we combine propagule and pre-reproductive?
#' FIXME: what about dormant stages?
#' FIXME: once above solved, then fix logic internally
#'
#' FIXME: Scott and Tamora added the extra parameter flags for switching whether
#' to include propagule and post-reproductive stages but these are not used as
#' yet. On reflection this is probably going too far and we should let the
#' matrix structure guide the process. However, in animals we may wish to
#' combine rep and post-rep to obtain a two-stage matrix model.
#'
#' @export
#' @param matF (matrix) fecundity matrix, rearranged so that non-reproductive
#' stages are in the final rows/columns
#' @param matFmu (matrix) mean fecundity matrix, rearranged so that
#' non-reproductive stages are in the final rows/columns
#' @param post (integer) indicator for post-reproductive stages
#' @param maxRep (integer) maximum number of reproductive stages
#' @param matrixStages (character) vector of stages, values are "prop"
#' (propagule), "active", and "dorm" (dormant)
#' @param includeProp (logical) include propagule stage. default: \code{TRUE}.
#' if \code{TRUE}, propagule stage (if present) is given back in result. If
#' \code{FALSE}, it's included into the pre-reproductive stage
#' @param includePost (logical) include post-reproductive stage. default:
#' \code{TRUE}. if \code{TRUE}, post-reproductive stage (if present) is given
#' back in result. If \code{FALSE}, it's included into the reproductive
#' stage
#' @examples
#' matF <- matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#'                  0, 0, 0, 0, 2.75, 1.75, 0, 0),
#'                  nrow = 5, byrow = FALSE)
#' matFmu <- matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#'                    0, 0, 0, 0, 0, 1.0409262, 0.5040727, 0.016433,
#'                    0.06956801),
#'                    nrow = 5, byrow = FALSE)
#' post <- numeric(0)
#' maxRep <- 5
#' matrixStages <- c("dorm", "active", "active", "active", "active")
#'
#' reprodStages(matF, matFmu, post, maxRep, matrixStages)
#'
#' # say what outputs you want
#' ## combine postrep and rep
#' reprodStages(matF, matFmu, post, maxRep, matrixStages, includeProp = FALSE,
#'   includePost = FALSE)
#' ## warn about prop?
#' reprodStages(matF, matFmu, post, maxRep, matrixStages, includeProp = FALSE)
#' ## NOT ALLOWED
#' #reprodStages(matF, matFmu, post, maxRep, matrixStages, c('prop', 'rep', 'postrep'))
#' ## NOT ALLOWED
#' #reprodStages(matF, matFmu, post, maxRep, matrixStages, c('prop', 'postrep'))
#' ## NOT ALLOWED
#' #reprodStages(matF, matFmu, post, maxRep, matrixStages, c('prerep', 'postrep'))
#' ## NOT ALLOWED
#' #reprodStages(matF, matFmu, post, maxRep, matrixStages, c('prop', 'prerep'))
reprodStages <- function(matF, matFmu, post, maxRep, matrixStages,
                         includeProp = TRUE, includePost = TRUE) {

  if ("prop" %in% matrixStages) {
    propStage <- which(matrixStages == "prop")
    if (!includeProp) {
      warning("Propagule stage exists, but includeProp is set to FALSE")
    }
  } else {
    propStage <- NA
  }

  # prerep
  matDim <- dim(matF)[1]
  Rep <- which(colSums(matFmu) > 0)
  if (min(Rep) == 1) {
    preRep <- NA
  } else if (!is.na(propStage[1]) && (min(Rep) - max(propStage) == 1)) {
    preRep <- NA
  } else {
    preRep <- min(which(matrixStages == "active")):(min(Rep) - 1)
  }

  # postrep
  if (length(post) == 0 && maxRep == matDim) {
    postRep <- NA
  } else {
    postRep <- (maxRep + 1):matDim
  }

  if (length(propStage) > 1) {
    propStages <- paste(propStage[1], "-", propStage[length(propStage)], sep = "")
  } else {
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

  c(propStages, preRepStages, repStages, postRepStages)
}
