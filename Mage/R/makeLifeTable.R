#' A function to produce a life table from a matrix population model
#' 
#' This function uses age-from-stage decompositions to calculate a life table
#' from a matrix population model
#' 
#' A detailed description of these methods can be found in section 5.3 of
#' Caswell (2001) and the supplementary information of Jones et al. (2014).
#' 
#' @param matU The U matrix - survival-dependent transitions (e.g. changes in
#' size/ageing/development)
#' @param matF The F matrix - sexual reproduction
#' @param matC The C matrix - clonal reproduction
#' @param startLife The stage of of the life cycle that structures the matrix
#' population model where life is considered to begin. This is usually 1 (i.e.
#' the first stage), but in species with a permanent propagule bank (e.g.
#' seedbank, or similar), users may wish to use the first non-permanent
#' propagule bank stage along the life cycle.
#' @param nSteps Number of time steps for which the life table be constructed.
#' This is on the same units as the matrix population model - see
#' MatrixPeriodicity in metadata of COMPADRE/COMADRE.
#' @return A `data.frame` with between 2 and 4 columns. If just `matU` is
#' provided, the `data.frame` has columns `x` (age), `lx` (survivorship).  If
#' `matF` is provided (in addition to `matU`), `mx` (age-specific sexual
#' reproduction) is included in the output `data.frame`. Likewise, if `matC` is
#' provided, an addition column, `cx` (age-specific clonal reproduction) is
#' included.
#' @note %% ~~further notes~~
#' @author Roberto Salguero-Gómez <rob.salguero@@sheffield.ac.uk>
#' 
#' Hal Caswell <h.caswell@@uva.nl>
#' 
#' Owen R. Jones <jones@@biology.sdu.dk>
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references
#' 
#' Caswell, H. (2001) Matrix Population Models: Construction, Analysis, and
#' Interpretation. Sinauer Associates; 2nd edition. ISBN: 978-0878930968
#' 
#' Caswell, H. (2006) Applications of Markov chains in demography. pp. 319-334
#' in A.N. Langville and W.J. Stewart (editors) MAM2006: Markov Anniversary
#' Meeting. Boson Books, Raleigh, North Caroline, USA
#' 
#' Jones, O.R. et al. (2014) Diversity of ageing across the tree of life.
#' Nature, 505(7482), 169–173
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' 
#' matU <- matrix(c(0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0.3, 0, 0, 0, 0, 0.1, 0.1), nrow = 4, byrow = T)
#' matF <- matrix(c(0, 0, 5, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.), nrow = 4, byrow = T)
#' matC <- matrix(c(0, 0, 0, 0, 0, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 4, byrow = T)
#' 
#' makeLifeTable(matU, matF, matC, startLife = 1, nSteps = 100)
#' 
#' @export makeLifeTable
makeLifeTable <- function(matU, matF = NULL, matC = NULL, startLife = 1, nSteps = 1000){
  
  matDim = ncol(matU)
  
  #Age-specific survivorship (lx) (See top function on page 120 in Caswell 2001):
  matUtemp = matU
  survivorship = array(NA, dim = c(nSteps, matDim))
  for (o in 1:nSteps){
    survivorship[o, ] = colSums(matUtemp %*% matU)
    matUtemp = matUtemp %*% matU
  }
  
  lx = survivorship[, startLife]
  lx = c(1, lx[1:(length(lx) - 1)])

  #Make room for dx and qx under assumption of 0.5 in age gap distributions
  
  #Start to assemble output object
  out = data.frame(x = 0:(length(lx)-1),lx = lx)
  
  if(!missing(matF)){
    if(sum(matF,na.rm=T)==0){
      warning("matF contains only 0 values")
    }
  #Age-specific fertility (mx, Caswell 2001, p. 120)
  ageFertility = array(0, dim = c(nSteps, matDim))
  fertMatrix = array(0, dim = c(nSteps, matDim))
  matUtemp2 = matU
  e = matrix(rep(1, matDim))
  for (q in 1:nSteps) {
    fertMatrix = matF %*% matUtemp2 * (as.numeric((ginv(diag(t(e) %*% matUtemp2)))))
    ageFertility[q, ] = colSums(fertMatrix)
    matUtemp2 = matUtemp2 %*% matU
  }  
  mx = ageFertility[, startLife]
  mx = c(0, mx[1:(length(mx) - 1)])
  out$mx = mx
  }
  
  if(!missing(matC)){
     if(sum(matC,na.rm=T)==0){
      warning("matC contains only 0 values")
    }
  #Age-specific clonality (cx)
  ageClonality = array(0, dim = c(nSteps, matDim))
  clonMatrix = array(0, dim = c(nSteps, matDim))
  matUtemp2 = matU
  e = matrix(rep(1, matDim))
  for (q in 1:nSteps) {
    clonMatrix = matC %*% matUtemp2 * (as.numeric((ginv(diag(t(e) %*% matUtemp2)))))
    ageClonality[q, ] = colSums(clonMatrix)
    matUtemp2 = matUtemp2 %*% matU
  }  
  cx = ageClonality[, startLife]
  cx = c(0, cx[1:(length(cx) - 1)])
  out$cx = cx
  }
  
  return(out)
  }
