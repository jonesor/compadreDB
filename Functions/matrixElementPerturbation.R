matrixElementPerturbation <- function(matU, matF, matC=NULL,pert=0.001){
  #Function to determine the cutoff age at quasi-convergence for lx and mx (Code adapted from H. Caswell's matlab code):
  
  matA=matU+matF+matC
  aDim=dim(matA)[1]
  fakeA=matA
  sensA=elasA=matrix(NA,aDim,aDim)
  lambda=Re(eigen(matA)$values[1])

  propU=matU/matA
    propU[is.nan(propU)]=NA
    propProg=propRetrog=propU
    propProg[upper.tri(propU,diag=T)]=NA
    propRetrog[lower.tri(propU,diag=T)]=NA
    propStasis=matrix(Diagonal(aDim)*diag(propU),aDim,aDim)
  propF=matF/matA
    propF[is.nan(propF)]=NA
  propC=matC/matA
    propC[is.nan(propC)]=NA

  for (i in 1:aDim){
    for (j in 1:aDim){
       fakeA=matA
       fakeA[i,j]=fakeA[i,j]+pert
       lambdaPert=eigen(fakeA)$values[1]
       sensA[i,j]=(lambda-lambdaPert)/(matA[i,j]-fakeA[i,j])
    }
  }

  elasA=sensA*matA/lambda
  
  out = data.frame("SStasis"=NA,"SProgression"=NA,"SRetrogression"=NA,"SFecundity"=NA,"SClonality"=NA,
                  "EStasis"=NA,"EProgression"=NA,"ERetrogression"=NA,"EFecundity"=NA,"EClonality"=NA)
    
    out$SStasis=Re(sum(sensA*propStasis,na.rm=T))
    out$SRetrogression=Re(sum(sensA*propRetrog,na.rm=T))
    out$SProgression=Re(sum(sensA*propProg,na.rm=T))
    out$SFecundity=Re(sum(sensA*propF,na.rm=T))
    out$SClonality=Re(sum(sensA*propC,na.rm=T))
    out$EStasis=Re(sum(elasA*propStasis,na.rm=T))
    out$EProgression=Re(sum(elasA*propProg,na.rm=T))
    out$ERetrogression=Re(sum(elasA*propRetrog,na.rm=T))
    out$EFecundity=Re(sum(elasA*propF,na.rm=T))
    out$EClonality=Re(sum(elasA*propC,na.rm=T))

  return(out) 
}