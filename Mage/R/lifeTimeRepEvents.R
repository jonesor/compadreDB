#' @export
#' @import MASS


lifeTimeRepEvents <- function(matU, matF, matC = F, startLife = 1){
  #Function to determine probability of reaching reproduction, age at 
  #maturity and reproductive lifespan (Code adapted from H. Caswell's 
  #matlab code, and Morris & Doak):

  if(missing(matU)){stop('matU missing')}
  if(missing(matF) & missing(matC)){stop('matF or matC missing. You must provide at least one')}
  if(sum(matF + matC,na.rm=T)==0){stop('matF and matC contains only 0 values')}
  #if(sum(matC,na.rm=T)==0){stop('matC contains only 0 values')}
  
  tryCatch({
    
  matDim <- dim(matU)[1]
  surv <- colSums(matU)
  
  out = NULL
  
  if (sum(matF)>0){
    
    fecLifeStages <- colSums(matF)
    fecLifeStages[which(fecLifeStages>0)] <- 1

    #Probability of survival to first sexual reprod event
	  Uprime <- matU
	  Uprime[ , which(fecLifeStages == 1)] <- 0
	  Mprime <- matrix(0, 2, matDim)
	  for (p in 1:matDim) {
	    if (fecLifeStages[p] == 1){ Mprime[2,p] <- 1 }else{Mprime[1, p] = 1 - surv[p]}
	  }
	  Bprime <- Mprime %*% (MASS::ginv(diag(matDim) - Uprime))
	  out$pFec <- Bprime[2, startLife]

    #Age at first sexual reproduction (LaFec; Caswell 2001, p 124)
	  D <- diag(c(Bprime[2,]))
	  Uprimecond <- D %*% Uprime %*% MASS::ginv(D)
	  expTimeReprod <- colSums(MASS::ginv(diag(matDim) - Uprimecond))
	  out$LaFec <- LaFec <- expTimeReprod[startLife]
    
    #Mean life expectancy conditional on entering the life cycle in the first reproductive stage
    firstFecLifeStage <- min(which(fecLifeStages == 1))
	  N <- solve(diag(matDim[1]) - matU)
	  out$meanLifeExpectancyFec <- colSums(N)[firstFecLifeStage]
	
    #Life expectancy from mean maturity
    out$remainingMatureLifeExpectancyFec <- colSums(N)[startLife] - LaFec
  }
  
  
  
  if (sum(matC)>0){
    cloLifeStages <- colSums(matC)
    cloLifeStages[which(cloLifeStages>0)] <- 1
    
    #Probability of survival to first clonal reprod event
    Uprime <- matU
    Uprime[ , which(cloLifeStages == 1)] <- 0
    Mprime <- matrix(0, 2, matDim)
    for (p in 1:matDim) {
      if (cloLifeStages[p] == 1){ Mprime[2,p] <- 1 }else{Mprime[1, p] = 1 - surv[p]}
    }
    Bprime <- Mprime %*% (MASS::ginv(diag(matDim) - Uprime))
    out$pClo <- Bprime[2, startLife]
    
    #Age at first clonal reproduction (LaClo; Caswell 2001, p 124)
    D <- diag(c(Bprime[2,]))
    Uprimecond <- D %*% Uprime %*% MASS::ginv(D)
    expTimeReprod <- colSums(MASS::ginv(diag(matDim) - Uprimecond))
    out$LaClo <- LaClo <- expTimeReprod[startLife]
    
    #Mean life expectancy conditional on entering the life cycle in the first reproductive stage
    firstCloLifeStage <- min(which(cloLifeStages == 1))
    N <- solve(diag(matDim[1]) - matU)
    out$meanLifeExpectancyClo <- colSums(N)[firstCloLifeStage]
    
    #Life expectancy from mean maturity
    out$remainingMatureLifeExpectancyClo <- colSums(N)[startLife] - LaClo
  
    
  }
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	return(out)
 }
