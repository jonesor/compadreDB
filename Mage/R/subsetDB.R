#' @export



#' A function to subset the COMPADRE/COMADRE database
#' 
#' This function allows users to subset the COMPADRE/COMADRE database by
#' logical argument.
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param sub An argument made using logical operators (see `subset`) with
#' which to subset the data base. Any of the variables contained in the
#' metadata part of the COMPADRE/COMADRE database may be used.
#' @param db The COMPADRE or COMADRE database object.
#' @return Returns a subset of the database, with the same structure, but where
#' the records in the metadata match the criteria given in the `sub` argument.
#' @note %% ~~further notes~~
#' @author Owen R. Jones <jones@@biology.sdu.dk>
#' 
#' Bruce Kendall
#' 
#' Rob Salguero-Goómez <rob.salguero@@zoo.ox.ac.uk>
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' \dontrun{
#' ssData <- subsetDB(compadre, MatrixDimension > 3)
#' ssData <- subsetDB(compadre, MatrixDimension > 3 & MatrixComposite == "Mean")
#' ssData <- subsetDB(comadre, Continent == "Africa" & Class == "Mammalia")
#' ssData <- subsetDB(comadre, SurvivalIssue < 1 & Altitude > 1000 & Altitude < 1500)
#' }
#' 
#' @export subsetDB
subsetDB <- function(db,sub){
  e <- substitute(sub)
  r <- eval(e, db$metadata, parent.frame())
  subsetID <- (1:length(r))[r & !is.na(r)]
 
  # First make a copy of the database.
  ssdb <- db

  # Subset the sub-parts of the database
  ssdb$metadata <- ssdb$metadata[subsetID,]
  ssdb$mat <- ssdb$mat[subsetID]
  ssdb$matrixClass <- ssdb$matrixClass[subsetID]

  # Version information is retained, but modified as follows.
  if("version" %in% names(ssdb)){
    ssdb$version$Version <- paste(ssdb$version$Version," - subset created on ",format(Sys.time(), "%b_%d_%Y"),sep="")
    ssdb$version$DateCreated <- paste(ssdb$version$DateCreated," - subset created on ",format(Sys.time(), "%b_%d_%Y"),sep="")
    ssdb$version$NumberAcceptedSpecies <- length(unique(ssdb$metadata$SpeciesAccepted))
    ssdb$version$NumberStudies <- length(unique(ssdb$metadata$SpeciesAuthor))
    ssdb$version$NumberMatrices <- length(ssdb$mat)
  }

  return(ssdb)
}
