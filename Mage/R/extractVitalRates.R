#' Function to extract the vital rates of a 3x3 matrix that has been
#' collapsed into pre-, reprod, and post-reprod
#'
#' @export
#' @param matU a matrix without reproduction and clonality
#' @param matF a matrix, just reproduction
#' @param collapse stuff
#' @examples
#' matU <- matrix(c(0.2581, 0.1613, 0.1935, 0.2258, 0.1613, 0.0408, 0.2857,
#'                  0.4286, 0.102, 0.0816, 0.0385, 0.0385, 0.2692, 0.2308,
#'                  0.3462, 0, 0.0625, 0.125, 0.25, 0.5625, 0.1061, 0.1608,
#'                  0.2637, 0.1801, 0.2058),
#'                  nrow = 5, byrow = FALSE)
#' matF <- matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#'                  0, 0, 0, 0, 2.75, 1.75, 0, 0),
#'                  nrow = 5, byrow = FALSE)
#' collapse <- c("1-2","3-4-5")
#' xx <- collapseMatrix(matU, matF, collapse)
#' extractVitalRates(xx$matU, xx$matF, collapse)
extractVitalRates <- function(matU, matF, collapse) {
  if (dim(matU)[1] != 4) {
    stop("This matrix is not 4x4!", call. = FALSE)
  }
  vitalRates <- rep(NA, 10)
  surv <- colSums(matU,na.rm = TRUE)
  surv[which(is.na(collapse))] <- NA
  s1 <- surv[1]
  s2 <- surv[2]
  s3 <- surv[3]
  s4 <- surv[4]
  matUIndep <- matU
  for (i in 1:dim(matU)[1]) {
    matUIndep[,i] <- matU[,i] / surv[i]
  }
  g21 <- matUIndep[2,1]
  g31 <- matUIndep[3,1]
  g32 <- matUIndep[3,2]
  g43 <- matUIndep[4,3]
  g34 <- matUIndep[3,4]

  f13 <- matF[1,3]
  f23 <- matF[2,3]
  f33 <- matF[3,3]

  vitalRates <- data.frame(s1,s2,s3,s4,g21,g31,g32,g43,g34,f13,f23,f33,
                           stringsAsFactors = FALSE)
  stats::setNames(vitalRates, c("s1","s2","s3","s4","g21","g31","g32","g43",
                            "g34","f13","f23","f33"))
}
