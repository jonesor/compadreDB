subsetDB2 <-
function(sub, db=comadre){

  #obtain subset IDs
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
