#' @export



#' A function to calculate measures of longevity.
#' 
#' A function to calculate the mean life expectancy and maximum longevity of
#' individuals in a matrix population model
#' 
#' This function applies Markov chain approaches to obtain mean life
#' expectancy, and a loop to calculate maximum longevity.
#' 
#' @param matU A matrix containing only survival-dependent processes (growth,
#' stasis, shrinkage).
#' @param startLife The first stage at which the author consider the beginning
#' of life.
#' @param initPop Initial population size.
#' @param run Number of iterations that the maximum longevity will be
#' calculated for.
#' @return This function applies Markov chain approaches to obtain mean life
#' expectancy, and a loop to calculate maximum longevity. Outputs are:
#' 
#' - 'eta': mean life expectancy conditional on entering the life cycle on life
#' stage described by 'startLife'.
#' 
#' - 'Max': maximum longevity observed by iterating a population vector with
#' 'initPop' individuals in stage 'startLife' up to 'run' times.
#' @note %% ~~further notes~~
#' @author Roberto Salguero-Gomez <rob.salguero@@zoo.ox.ac.uk>
#' 
#' Hal Caswell <hcaswell@@whoi.edu>
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references Caswell, H. (2001) Matrix Population Models: Construction,
#' Analysis, and Interpretation. Sinauer Associates; 2nd edition. ISBN:
#' 978-0878930968
#' 
#' Doak & Morris (2002) BOOK INSERT CITATION
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' matU <- matrix (c(0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0.3, 0, 0, 0, 0, 0.1, 0.1), nrow = 4, byrow = T)
#' 
#' longevity(matU, startLife = 1, initPop = 100, run = 1000)
#' 
longevity <- function(matU, startLife = 1, initPop = 100, run = 1000){
  #Function to calculate mean life expectancy and maximum longevity from
  # H. Caswell's matlab code, and Morris & Doak:

  if(missing(matU)){stop('matU missing')}
  
  out = NULL
  
  matDim <- dim(matU)[1]
  
  #Mean life expectancy
  N <- solve(diag(matDim) - matU)
  out$eta <- colSums(N)[startLife]
  
  #Maximum longevity up to 'run' interations.
  popVector <- c(rep(0,matDim))
  popVector[startLife] <- initPop
  lifespanLeftover=matrix(0,run,1)
  for (n in 1:1000)	{
    lifespanLeftover[n]=sum(popVector)
    popVector=matU%*%popVector
  }
  
  out$Max <- max(which(lifespanLeftover>1))
  if (out$Max==Inf) {out$Max=run}
  
	return(out)
 }
