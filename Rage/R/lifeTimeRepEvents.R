#' A function to calculate the timing of lifetime reproductive events.
#' 
#' A function to calculate the timing of lifetime reproductive events such as
#' the probability of achieving maturity, age at first reproduction, mean life
#' expectancy conditional on maturity, and life expectancy for mature
#' individuals
#' 
#' This function applies Markov chain approaches to decompose various moments
#' of along the age-based reproduction of individuals in a matrix population
#' model.
#' 
#' @param matU A matrix containing only survival-dependent processes ( growth,
#' stasis, shrinkage).
#' @param matF A matrix containing only sexual reproduction, with zeros
#' elsewhere.
#' @param matC A matrix containing only clonal reproduction, with zeros
#' elsewhere.
#' @param startLife The first stage at which the author consider the beginning
#' of life.
#' @return This function applies Markov chain approaches to decompose various
#' moments of along the age-based reproduction of individuals in a matrix
#' population model. When both a 'matF' and a 'matC' are provided, the
#' following outputs are calculated for both independently, and these are
#' differentiated with the suffix "Fec" or "Clo", respectively:
#' 
#' - 'p': probability of achiving maturity, sexual or clonal.
#' 
#' - 'La': mean age at maturity (in the same units as the matrix population
#' model).
#' 
#' - 'meanLifeExpectancy': mean life expectancy conditional on entering the
#' life cycle in the first reproductive stage
#' 
#' - 'remainingMatureLifeExpectancy': Life expectancy from mean maturity. This
#' is mean life expectancy - mean age at maturity ('La' above). This value can
#' be negative because both mean life expectancy and mean age at maturity are
#' means of their respective distributions.
#' @note %% ~~further notes~~
#' @author Roberto Salguero-Gomez <rob.salguero@@zoo.ox.ac.uk>
#' 
#' Hal Caswell <hcaswell@@whoi.edu>
#' 
#' Owen R. Jones <jones@@biology.sdu.dk>
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references Caswell, H. (2001) Matrix Population Models: Construction,
#' Analysis, and Interpretation. Sinauer Associates; 2nd edition. ISBN:
#' 978-0878930968
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' matU <- matrix (c(0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0.3, 0, 0, 0, 0, 0.1, 0.1), nrow = 4, byrow = T)
#' matF <- matrix (c(0, 0, 5, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 4, byrow = T)
#' matC <- matrix (c(0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0, 0), nrow = 4, byrow = T)
#' 
#' lifeTimeRepEvents(matU, matF, matU, startLife = 1)
#' 
#' @export lifeTimeRepEvents
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
