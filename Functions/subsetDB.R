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

  # Version information is retained.
  return(ssdb)
}

