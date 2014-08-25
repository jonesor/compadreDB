lifeTimeClonallEvents <- function(matU, matC, conv = 0.05, startLife = 1, nSteps = 1000){
  #Function to determine probability of reaching clonal reprouction, age at clonal propagation and clonal propagation lifespan (Code adapted from H. Caswell's matlab code):
  
  uDim = dim(matU)
  eig = eigen.analysis(matU)
  
  
}