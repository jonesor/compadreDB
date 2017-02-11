#' Collapse matrices
#'
#' @export
#' @param matU a matrix without reproduction and clonality
#' @param matF a matrix, just reproduction
#' @param collapse (character) vector indicating which stages must be collapsed
#  and retained as are
#' @examples
#' matU <- matrix(c(0.2581, 0.1613, 0.1935, 0.2258, 0.1613, 0.0408, 0.2857,
#'                  0.4286, 0.102, 0.0816, 0.0385, 0.0385, 0.2692, 0.2308,
#'                  0.3462, 0, 0.0625, 0.125, 0.25, 0.5625, 0.1061, 0.1608,
#'                  0.2637, 0.1801, 0.2058),
#'                  nrow = 5, byrow = FALSE)
#' matF <- matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#'                  0, 0, 0, 0, 2.75, 1.75, 0, 0),
#'                  nrow = 5, byrow = FALSE)
#' collapse1 <- c("1-2","3-4","5")
#' (out <- collapseMatrix(matU, matF, collapse2))
#' eigen(out$matA)
#' eigen(out$matU)
#' eigen(out$matF)
#'
#' #collapse2 <- c("1-2","3-4-5")
#' #collapse3 <- c("1-2-3-4-5")
collapseMatrix <- function(matU, matF, collapse) {
  matA <- matU + matF
  collapseUnique <- collapse
  originalDim <- dim(matA)[1]
  collapseDim <- length(collapseUnique)
  P <- matrix(0, nrow = collapseDim , ncol = originalDim)

  splitCollapseUnique <- strsplit(collapse, "-")
  for (i in 1:collapseDim) {
    columns <- as.numeric(splitCollapseUnique[[i]])
    if (!is.na(columns[1])) {
      P[i,(columns[1]:columns[length(columns)])] <- 1
    }
  }

  Q <- t(P)
  w <- Re(eigen(matA)$vectors[,which(Re(eigen(matA)$values) == max(Re(eigen(matA)$values)))])
  w <- w/sum(w)

  columns <- which(colSums(Q) > 1)
  for (j in columns) {
    rows <- which(Q[,j] == 1)
    for (i in rows) {
      Q[i,j] <- w[i]/sum(w[rows])
    }
  }
  collapseA <- P %*% matA %*% Q
  collapseU <- P %*% matU %*% Q
  collapseF <- P %*% matF %*% Q
  list(matA = collapseA, matU = collapseU, matF = collapseF)
}
