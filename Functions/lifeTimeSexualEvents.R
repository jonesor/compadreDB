lifeTimeSexualEvents <- function(matU, matF, conv = 0.05, startLife = 1, nSteps = 1000){
  #Function to determine probability of reaching reproduction, age at maturity and reproductive lifespan (Code adapted from H. Caswell's matlab code):
  
  uDim = dim(matU)
  eig = eigen.analysis(matU)
  
  
}