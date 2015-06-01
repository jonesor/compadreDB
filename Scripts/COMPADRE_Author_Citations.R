# ---------------------------------------------------------------------------- #
# - FILE NAME:   COMPADRE_and_COMADRE_Author_Citations.R    		
# - VERSION:     01
# - DATE:        July 16 2014
# - MODIFIED:    June 01 2015
# - DESCRIPTION: Code to extract full citation fror studies where DOI is 
#                available in COMPADRE from The Plant List 
#                (www.theplantlist.org), using the R library Taxonstand.
# - AUTHORS:     Rob Salguero-Gomez, Owen Jones and Fernando Colchero
# ---------------------------------------------------------------------------- #

# START SCRIPT:
# Delete memory to avoid potential problems with saved variables
rm(list=ls(all=TRUE))

# Required libraries
install.packages("dev_tools")
library(dev_tools)
install_github("rmetadata", "ropensci")
library(rmetadata)
library(ropensci)
library(rplos)

# Set the reading directory to where the COMPADRE or COMADRE data have been saved. For 
# example, if COMPADRE are in ''C:/Documents/COMPADRE/Data/'':
setwd("C:/Documents/COMPADRE/Data/")

# Read in data. In 3 columns - author, journal year (title would be useful but 
# was not recorded)
papers <- read.csv("COMPADRE Data - July 20 2014.csv")

# Paste together a query
papers$qry1 <- paste(papers$Author, papers$Journal, papers$Year)

#Create empty vector for DOI
papers$doi <- NA

# Loop submitting query to crossref.
for (i in 1:nrow(papers)){
  qryOutcome <- crossref_search_free(papers$qry1[i])
  if(class(qryOutcome) != "NULL"){
    if (!"doi" %in% names(qryOutcome)) {
      qryDOI = NA
    } else {
      qryDOI <- qryOutcome$doi
    }
  } else {
    qryDOI <- NA
  }
  papers$doi[i] <- as.character(qryDOI)
}

# Attempt to get the Title, author, journal and year from crossref using the DOI.
# This is to be used to double-check that the data are good (e.g. the paper does not have a bizarre physcis title)
papers$CrossRefTitle <- NA
papers$CrossRefAuthor <- NA
papers$CrossRefJournal <- NA
papers$CrossRefYear <- NA

for (i in 1:nrow(papers)){
  temp <- try(crossref(doi = papers$doi[i]), silent = TRUE)
  if (class(temp) == "try-error") {
    papers$CrossRefTitle[i] <- NA
    papers$CrossRefAuthor[i] <- NA
    papers$CrossRefJournal[i] <- NA
    papers$CrossRefYear[i] <- NA
  } else {
    papers$CrossRefTitle[i] <- temp$title
    papers$CrossRefAuthor[i] <- paste(temp$author, collapse = ", ")
    papers$CrossRefJournal[i] <- temp$journal
    papers$CrossRefYear[i] <- temp$year
  }
}

# Write out results.
write.csv(papers, file = "COMPADRE_DATA_UPDATED_CROSSREF_RESULTS BY USER.csv", row.names = FALSE)
