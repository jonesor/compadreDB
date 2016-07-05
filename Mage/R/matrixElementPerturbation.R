#' @export

matrixElementPerturbation <- function(matU, matF, matC=NULL, pert = 0.001){
  matA <- matU + matF + matC
  aDim <- nrow(matA)
  fakeA <- matA
  sensA <- elasA <- matrix(NA,aDim,aDim)
  lambda <- Re(eigen(matA)$values[1])
  
  propU <- matU / matA
  propU[is.nan(propU)] <- NA
  propProg <- propRetrog <- propU
  propProg[upper.tri(propU, diag = TRUE)] <- NA
  propRetrog[lower.tri(propU, diag = TRUE)] <- NA
  propStasis <- matrix(diag(aDim) * diag(propU), aDim, aDim)
  propF <- matF / matA
  propF[is.nan(propF)] <- NA
  propC <- matC / matA
  propC[is.nan(propC)] <- NA
  
  for (i in 1:aDim){
    for (j in 1:aDim){
      fakeA <- matA
      fakeA[i, j] <- fakeA[i,j] + pert
      lambdaPert <- eigen(fakeA)$values[1]
      sensA[i, j] <- (lambda - lambdaPert) / (matA[i, j] - fakeA[i, j])
    }
  }
  
  sensA <- Re(sensA)
  elasA <- sensA * matA / lambda
  
  out <- data.frame("SStasis" = NA, "SProgression" = NA, "SRetrogression" = NA,
                    "SFecundity" = NA, "SClonality" = NA, "EStasis" = NA,
                    "EProgression" = NA, "ERetrogression"=NA,
                    "EFecundity" = NA, "EClonality" = NA)
  
  out$SStasis <- sum(sensA * propStasis, na.rm = TRUE)
  out$SRetrogression <- sum(sensA * propRetrog, na.rm = TRUE)
  out$SProgression <- sum(sensA * propProg, na.rm = TRUE)
  out$SFecundity <- sum(sensA * propF, na.rm=TRUE)
  out$SClonality <- sum(sensA * propC, na.rm=TRUE)
  out$EStasis <- sum(elasA * propStasis,na.rm=TRUE)
  out$EProgression <- sum(elasA * propProg, na.rm=TRUE)
  out$ERetrogression <- sum(elasA * propRetrog, na.rm=TRUE)
  out$EFecundity <- sum(elasA * propF, na.rm=TRUE)
  out$EClonality <- sum(elasA * propC, na.rm=TRUE)
  
  return(out) 
}