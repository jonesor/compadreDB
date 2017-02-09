#' Standardized vital rates
#'
#' @export
#' @param matU a matrix without reproduction and clonality
#' @param matF a matrix, just reproduction
#' @param matFmu a matrix
#' @param matrix_stages a vector of matrix stage statuses
#' @examples
#' # compadre
#' # load("~/COMPADRE_v.4.0.1.RData")
#' #mats <- compadre$metadata$SpeciesAccepted
#'
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
#'
#' foo(
#'  matU = matU,
#'  matF = matF,
#'  matFmu = matFmu,
#'  matrix_stages = c("dorm", "active", "active", "active", "active")
#' )
standardizedVitalRates <- function(matU, matF, matFmu, matrix_stages) {
  # put non-reproductive to the end of the matrix
  rearr <- rearrangeMatrix(matU, matF, matFmu)
  matFmu <- rearr$matFmu

  # non-reproductive stages
  rearr2 <- rearrangeMatrix(matU, matF, matFmu)
  post <- rearr2$nonRepInterRep
  maxRep <- rearr2$maxRep

  # defines which columns need to be collapsed for each of the four stages
  collapse <- reprodStages(matF, matFmu, post, maxRep, matrix_stages)

  matUcollapse <- collapseMatrix(matU, matF, collapse = collapse)$matU
  matUcollapse[is.na(collapse), is.na(collapse)] <- NA

  matFcollapse <- collapseMatrix(matU, matF, collapse = collapse)$matF
  matFcollapse[is.na(collapse), is.na(collapse)] <- NA

  extractVitalRates(matU = matUcollapse, matF = matFcollapse, collapse)
}
