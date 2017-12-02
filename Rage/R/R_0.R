#' @export

R_0 <- function(matU, matF, matC=F){
  
  #Demetrius entropy (S):
  
  R0=NULL
  matDim <- dim(matA)[1]
  N <- solve(diag(matDim)-matU)  #Fundamental matrix, which states the amount of time units spent in each stage on average
  
  if (sum(matF)>0){
    R0_matF <- matF%*%N
    R0$Fec <- R0_matF[1,1]
  }
  if (sum(matC)>0){
    R0_matC <- matC%*%N
    R0$Clo <- R0_matC[1,1]
  }
  if (sum(matF)>0 & sum(matC)>0){
    matFC <- matF + matC
    R0_matFC <- matFC%*%N
    R0$FecClo <- R0_matFC[1,1]
  }
  
return(R0)
}