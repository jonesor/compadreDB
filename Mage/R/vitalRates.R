#' A function to derive vital rates from the matrix population model.
#' 
#' A function to derive vital rates from the matrix population model for the
#' separate demographic processes.
#' 
#' This function decomposes vital rates of survival, progression,
#' retrogression, sexual reproduction and clonal reproduction according to
#' various ways of weighted means and organization of stages along the life
#' cycle represented in the matrix population model.
#' 
#' @param matU A matrix containing only survival-dependent processes ( growth,
#' stasis, shrinkage).
#' @param matF A matrix containing only sexual reproduction, with zeros
#' elsewhere.
#' @param matC A matrix containing only clonal reproduction, with zeros
#' elsewhere.
#' @param splitStages Splits vital rates according to some pre-determined
#' criteria (below).
#' @param weighted Allows to weight mean vital rates according to various
#' criteria (below).
#' @return - 'Weighted': This argument allows to weight mean values of vital
#' rates (survival 'surv', progression 'prog', retrogression 'retr', sexual
#' reproduction 'fec' and clonal reproduction 'clo') with an equal contribution
#' for all stages (default), by the stable st/age distribution ('SSD'), or by a
#' given population vector chosen by the user, so long as it is congruent with
#' the dimensions of the chosen 'matU', 'matF', and 'matC'.
#' 
#' - 'splitStages': This argument allows to split the values of vital rates
#' according to recognizable stages in the matrix. When 'all', all vital rates
#' are averaged all existing stages, if 'ontogeny', they are averaged as
#' juveniles ('Juv') and adults ('Adu'), and if by 'MatrixClassOrganized', it
#' takes a vector with the pre-established stages of
#' 'compadre$matrixClass[[i]]$MatrixClassOrganized' or
#' 'compadre$matrixClass[[i]]$MatrixClassOrganized', where 'i' is the index of
#' the chosen study in 'COMPADRE' or 'COMADRE'.
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
#' #Vital rate outputs without weights:
#' vitalRates(matU, matF, matC, splitStages = 'all', weighted = FALSE)
#' vitalRates(matU, matF, matC, splitStages = 'ontogeny', weighted = FALSE)
#' vitalRates(matU, matF, matC, splitStages = c('prop', 'active', 'active', 'active'), weighted = FALSE)
#' 
#' 
#' 
#' #Vital rate outputs weighted by the stable stage distribution of 'matA':
#' vitalRates(matU, matF, matC, splitStages = 'all', weighted = 'SSD')
#' vitalRates(matU, matF, matC, splitStages = 'ontogeny', weighted = 'SSD')
#' vitalRates(matU, matF, matC, splitStages = c('prop', 'active', 'active', 'active'), weighted = 'SSD')
#' 
#' #Vital rate outputs weighted by a chosen population vector of initial conditions:
#' initialConditions <- c(100, 10, 0, 1)
#' 
#' vitalRates(matU, matF, matC, splitStages = 'all', weighted = initialConditions)
#' vitalRates(matU, matF, matC, splitStages = 'ontogeny', weighted = initialConditions)
#' vitalRates(matU, matF, matC, splitStages = c('prop', 'active', 'active', 'active'), weighted = initialConditions)
#' 
#' @export
vitalRates <- function(matU, matF, matC = FALSE, splitStages = FALSE, weighted = FALSE){
  #Function to quantify vital rates values
  
  if(missing(matU)){stop('matU missing')}
  if(missing(matF) & missing(matC)){warning('matF or matC missing. These have been coerced to matrices of zero')}
  #if (sum(weighted)>0 & length(weighted) != dim(matU)[i]){stop('Population vector does not agree with matrix dimension')}
  
  matDim <- dim(matU)[1]
  matA <- matU + matF + matC
  
  surv <- colSums(matU)
  fec <- colSums(matF)
  clo <- colSums(matC)
  
  matUIndep <- matrix(NA, matDim, matDim)
  for (i in 1:matDim) {matUIndep[,i] = matU[,i]/surv[i]}
  prog <- retr <- matUIndep
  prog[is.nan(prog)] <- 0
  retr[is.nan(retr)] <- 0
  
  prog[which(upper.tri(matUIndep, diag = TRUE))]=0
  retr[which(lower.tri(matUIndep, diag = TRUE))]=0
  
  prog <- colSums(prog)
  retr <- colSums(retr)
  
  if (weighted[1] == FALSE) {
    weight <- rep(1,matDim)
  }
  
  if (weighted[1] == 'SSD') {
    weight <- Re(eigen(matA)$vectors[,which.max(Re(eigen(matA)$values))])
  }
  
  weight <- weight/sum(weight)
  
  surv1 <- surv * weight
  fec1  <- fec  * weight
  clo1  <- clo  * weight
  prog1 <- prog * weight
  retr1 <- retr * weight
  
  out = NULL
  
  if (splitStages[1] == 'all'){
    out$surv <- sum(surv1)
    out$retr <- sum(retr1)
    out$prog <- sum(prog1)
    out$fec  <- sum(fec1)
    out$clo  <- sum(clo1)
  }
  
  if (splitStages[1] == 'ontogeny'){
    adu <- which(colSums(matF)>0)  #This classification does not accoutn for non- and post-reproductive
    juv <- which(colSums(matF)==0)
    
    out$survJuv <- mean(surv1[juv], na.rm=T)
    out$retrJuv <- mean(retr1[juv], na.rm=T)
    out$progJuv <- mean(prog1[juv], na.rm=T)
    out$cloJuv  <- mean(clo1[juv], na.rm=T)
    
    out$survAdu <- mean(surv1[adu], na.rm=T)
    out$retrAdu <- mean(retr1[adu], na.rm=T)
    out$progAdu <- mean(prog1[adu], na.rm=T)
    out$fecAdu  <- mean(fec1[adu], na.rm=T)
    out$cloAdu  <- mean(clo1[adu], na.rm=T)
    }
  
  if (splitStages[1] %in% c('prop','active','dorm')){
    prop <- which(splitStages=="prop")
    active <- which(splitStages=="active")
    dorm <- which(splitStages=="dorm")
    
    out$survProp <- mean(surv1[prop], na.rm=T)
    out$progProp <- mean(prog1[prop], na.rm=T)
    
    out$survActive <- mean(surv1[active], na.rm=T)
    out$retrActive <- mean(retr1[active], na.rm=T)
    out$progActive <- mean(prog1[active], na.rm=T)
    out$fecActive  <- mean(fec1[active], na.rm=T)
    out$cloActive  <- mean(clo1[active], na.rm=T)
    
    out$survDorm <- mean(surv1[dorm], na.rm=T)
    out$retrDorm <- mean(retr1[dorm], na.rm=T)
    out$progDorm <- mean(prog1[dorm], na.rm=T)
  }
  
	return(out)
 }
