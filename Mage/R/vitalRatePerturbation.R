#' A function to perform perturbation of vital rates of the matrix model
#' 
#' A function to perform perturbation of vital rates of the matrix model
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param matU %% ~~Describe \code{matU} here~~
#' @param matF %% ~~Describe \code{matF} here~~
#' @param matC %% ~~Describe \code{matC} here~~
#' @param pert %% ~~Describe \code{pert} here~~
#' @return %% ~Describe the value returned 
#' @note %% ~~further notes~~
#' @author Roberto Salguero-Gomez <r.salguero@@sheffield.ac.uk>
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references %% ~put references to the literature/web site here ~
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' ##---- Should be DIRECTLY executable !! ----
#' ##-- ==>  Define data, use random,
#' ##--	or do  help(data=index)  for the standard data sets.
#' 
#' ## The function is currently defined as
#' 
#' 
#' @export vitalRatePerturbation
vitalRatePerturbation <- function(matU, matF, matC=NULL,pert=0.001){
  #Function to calculate vital rate level sensitivities and elasticities
  
  matA=matU+matF+matC
  aDim=dim(matA)[1]
  fakeA=matA
  sensA=elasA=matrix(NA,aDim,aDim)
  lambda=Re(eigen(matA)$values[1])

  propU=matU/matA
    propU[is.nan(propU)]=0
    propProg=propRetrog=propU
    propProg[upper.tri(propU,diag=T)]=0
    propRetrog[lower.tri(propU,diag=T)]=0
    propStasis=matrix(diag(aDim)*diag(propU),aDim,aDim)
  propF=matF/matA
    propF[is.nan(propF)]=0
  propC=matC/matA
    propC[is.nan(propC)]=0

  #for (i in 1:aDim){
  #  for (j in 1:aDim){
  #     fakeA=matA
  #     fakeA[i,j]=fakeA[i,j]+pert
  #     lambdaPert=eigen(fakeA)$values[1]
  #     sensA[i,j]=(lambda-lambdaPert)/(matA[i,j]-fakeA[i,j])
  #  }
  #}
  #sensA=Re(sensA)

	sensA=eigen.analysis(matA,zero=F)$sensitivities

    #Survival-independent A matrix
    uIndep=matrix(NA,aDim,aDim)
    u=colSums(matU)
    for (j in which(u>0)) uIndep[,j]=matA[,j]/u[j]
   	
   	sensSigmaA=uIndep*sensA
    
    #Little fix for semelparous species
      uPrime=u
      #uPrime[u==0]=0.001
    elasSigmaA=t(t(sensSigmaA)*uPrime)/lambda
    
  	elasA=sensA*matA/lambda

  #Extracting survival vital rate
    uDistrib=matrix(0,ncol=aDim,nrow=aDim)
    for (j in which(u>0)) uDistrib[,j]=matU[,j]/u[j]
  #Extracting fecundity vital rates:
    f=colSums(matF)
    fDistrib=matrix(0,ncol=aDim,nrow=aDim)
    for (j in which(f>0)) fDistrib[,j]=matF[,j]/f[j]
  #Extracting clonality vital rates:
    c=colSums(matC)
    cDistrib=matrix(0,ncol=aDim,nrow=aDim)
    for (j in which(c>0)) cDistrib[,j]=matC[,j]/c[j]
    

  SuDistrib=sensA*uDistrib
  SfDistrib=sensA*fDistrib
  ScDistrib=sensA*cDistrib


out = data.frame("SSurvival"=NA,"SGrowth"=NA,"SShrinkage"=NA,"SReproduction"=NA,"SClonality"=NA,
                  "ESurvival"=NA,"EGrowth"=NA,"EShrinkage"=NA,"EReproduction"=NA,"EClonality"=NA)
  

  #Still to be done
  out$SSurvival=sum(sensSigmaA,na.rm=T)
  out$SGrowth=sum(sensA*uDistrib*propProg,na.rm=T)
  out$SShrinkage=sum(sensA*uDistrib*propRetrog,na.rm=T)
  out$SReproduction=sum(sensA*fDistrib*propF,na.rm=T)
  out$SClonality=sum(sensA*cDistrib*propC,na.rm=T)

  EuDistrib=sensA*uDistrib*matrix(u,nrow=aDim,ncol=aDim,byrow=T)/lambda
  EfDistrib=sensA*fDistrib*matrix(f,nrow=aDim,ncol=aDim,byrow=T)/lambda
  EcDistrib=sensA*cDistrib*matrix(c,nrow=aDim,ncol=aDim,byrow=T)/lambda

  out$ESurvival=sum(elasSigmaA,na.rm=T)
  out$EGrowth=sum(EuDistrib*propProg,na.rm=T)
  out$EShrinkage=sum(EuDistrib*propRetrog,na.rm=T)
  out$EReproduction=sum(EfDistrib*propF,na.rm=T)
  out$EClonality=sum(EcDistrib*propC,na.rm=T)

 return(out) 
}
