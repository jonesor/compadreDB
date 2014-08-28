vitalRatePerturbation <- function(matU, matF, matC=NULL,pert=0.001){
  #Function to calculate vital rate level sensitivities and elasticities
  
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

  sensA=Re(sensA)
  elasA=sensA*matA/lambda

  #Extracting survival vital rate
    u=colSums(matU)
    uDistrib=matrix(NA,ncol=aDim,nrow=aDim)
    for (j in which(u>0)) uDistrib[,j]=matU[,j]/u[j]
  #Extracting fecundity vital rates:
    f=colSums(matF)
    fDistrib=matrix(NA,ncol=aDim,nrow=aDim)
    for (j in which(f>0)) fDistrib[,j]=matF[,j]/f[j]
  #Extracting clonality vital rates:
    c=colSums(matC)
    cDistrib=matrix(NA,ncol=aDim,nrow=aDim)
    for (j in which(c>0)) cDistrib[,j]=matC[,j]/c[j]
    
  SuDistrib=sensA*uDistrib
  SfDistrib=sensA*fDistrib
  ScDistrib=sensA*cDistrib

out = data.frame("SSurvival"=NA,"SGrowth"=NA,"SShrinkage"=NA,"SReproduction"=NA,"SClonality"=NA,
                  "ESurvival"=NA,"EGrowth"=NA,"EShrinkage"=NA,"EReproduction"=NA,"EClonality"=NA)
  

  #Still to be done
  #out$SSurvival=sum(uDistrib,na.rm=T)
  out$SGrowth=sum(sensA*uDistrib*propProg,na.rm=T)
  out$SShrinkage=sum(sensA*uDistrib*propRetrog,na.rm=T)
  out$SReproduction=sum(sensA*fDistrib*propF,na.rm=T)
  out$SClonality=sum(sensA*cDistrib*propC,na.rm=T)

  EuDistrib=sensA*uDistrib*matrix(u,nrow=aDim,ncol=aDim,byrow=T)/lambda
  EfDistrib=sensA*fDistrib*matrix(f,nrow=aDim,ncol=aDim,byrow=T)/lambda
  EcDistrib=sensA*cDistrib*matrix(c,nrow=aDim,ncol=aDim,byrow=T)/lambda

  #out$ESurvival=sum(uDistrib,na.rm=T)
  out$EGrowth=sum(EuDistrib*propProg,na.rm=T)
  out$EShrinkage=sum(EuDistrib*propRetrog,na.rm=T)
  out$EReproduction=sum(EfDistrib*propF,na.rm=T)
  out$EClonality=sum(EcDistrib*propC,na.rm=T)

 return(out) 
}