
# ---------------------------------------------------------------------------- #
# - FILE NAME:   COMPADRE_Update_Taxonomy_TPL.R				
# - VERSION:     01
# - DATE:        July 20 2014
# - DESCRIPTION: Code to extract the currently accepted species names 
#                (variable SpeciesAccepted) from The 'Plant List' 
#                (www.theplantlist.org), using the R library Taxonstand.
# - AUTHORS:     Rob Salguero-Gomez, Owen Jones & Fernando Colchero
# ---------------------------------------------------------------------------- #

# START SCRIPT:
# Delete memory to avoid potential problems with saved variables
rm(list=ls(all=TRUE))

# Required libraries
library(Taxonstand)

# Set the reading directory to where the COMPADRE data have been saved. For 
# example, if COMPADRE are in ''C:/Documents/COMPADRE/Data/'':
setwd("C:/Documents/")

# Read in the COMPADRE data:
COMPADRE <- read.csv("COMPADRE Data - July 20 2014.csv")

# Extract species names used by authors:
speciesAuthor <- as.character(COMPADRE$SpeciesAuthor)

# Some basic modifications in the name before passing them through TPL
speciesAuthor <- gsub("cf. ", "", speciesAuthor)
speciesAuthor <- gsub("_[0-9]*$", "", speciesAuthor)
speciesAuthor <- gsub("_", " ", speciesAuthor)
speciesAuthor2 <- strsplit(speciesAuthor, " ")
maxLen <- max(sapply(speciesAuthor2, length))
speciesAuthor3 <- t(sapply(speciesAuthor2, 
                                   function(x) 
                                     c(x, rep("", maxLen - length(x)))))
speciesAuthor4 <- cbind(speciesAuthor,speciesAuthor3)

colnames(speciesAuthor4) <- c("Full.name", "Genus", "Species", "Var",
                                      "Intraspecific")
speciesAuthorFinal <- data.frame(speciesAuthor4)
speciesAuthorFinalChecked <- 
  TPL(genus = speciesAuthorFinal$Genus, 
      species = speciesAuthorFinal$Species,
      infrasp = speciesAuthorFinal$Intraspecific, corr = TRUE)

# Suggested changes in taxonomic names
ind <- which(as.character(speciesAuthorFinalChecked$Genus) != 
  as.character(speciesAuthorFinalChecked$New.Genus) | 
    as.character(speciesAuthorFinalChecked$Species) !=  
    as.character(speciesAuthorFinalChecked$New.Species))
speciesSuggestedChanges <- speciesAuthorFinalChecked[ind, ]
speciesSuggestedChanges[,c('Genus', 'Species', 'Taxonomic.status',
                                     'New.Genus', 'New.Species')]

#Create TPL-Accepted Full.name
c.infra <- as.character(speciesAuthorFinalChecked$New.Infraspecific)
c.infra[is.na(speciesAuthorFinalChecked$New.Infraspecific)] <- ""
c.gen <- as.character(speciesAuthorFinalChecked$New.Genus)
c.sp <- as.character(speciesAuthorFinalChecked$New.Species)

trim <- function (x) gsub("^\\s+|\\s+$", "", x)
speciesAccepted <- trim((paste(c.gen,c.sp,c.infra,sep=" ")))

COMPADRE$SpeciesAccepted <- speciesAccepted

# Saving to output
write.csv(COMPADRE, "COMPADRE_DATA_UPDATED_TAXONOMY_BY_USER.csv", row.names = FALSE)
