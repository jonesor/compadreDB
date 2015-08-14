# Function to subset COMPADRE/COMADRE database
#Author: Owen Jones
#
#Example use
#
# x<-subsetDB(comadre,sub = which(comadre$metadata$Class == "Mammalia"))
#

subsetDB <- function(db=comadre,sub=1:100){
  
  #is there a way to put the subset command into the arguments for this function?
  subsetID <- sub

  # First make a copy of the database
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

