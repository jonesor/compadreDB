#' @export
#' @import popbio

qsdConverge <- function(matU, conv = 0.05, startLife = 1, nSteps = 1000){
  
  #Function to determine the cutoff age at quasi-convergence for lx and mx (Code adapted from H. Caswell's matlab code):
  requireNamespace("popbio")
  
  uDim = dim(matU)
  eig = eigen.analysis(matU)
  qsd = eig$stable.stage
  qsd = as.numeric(t(matrix(qsd / sum(qsd))))
  
  #Set up a cohort
  nzero = rep(0, uDim[1]) #Set a population vector of zeros
  nzero[startLife] = 1 #Set the first stage to = 1
  n = nzero #Rename for convenience
  
  #Iterate the cohort (n= cohort population vector, p = proportional structure)
  dist = p = NULL
  survMatrix1 <- matU
  for (j in 1:nSteps){ #j represent years of iteration
    p = n / sum(n) #Get the proportional distribution
    dist[j] = 0.5 * (sum(abs(p - qsd)))
    n = survMatrix1 %*% n #Multiply the u and n matrices to iterate
  }
  #Find the ages for convergence to conv. (default = 0.05).
  #i.e. within 5% of the QSD.
  convage = min(which(dist < conv))
  return(convage) 
}
