#Function to compare two COMPADRE or COMADRE data objects
#Author: Owen Jones
#
#Example use
#
#compadreCompare(db1=compadreOld, db2 = compadreNew)
#
#

dbCompare <- function(db1,db2){ 

#Quick summary
cat("Quick Summary\n")
#File 1
uniqueSource <- unique(paste(db1$metadata$Authors,db1$metadata$Journal,db1$metadata$YearPublication))                      
cat(paste("File-1 contains the demographic and associated data for ", 
       length(uniqueSource), " studies, corresponding to ",
       length(unique(db1$metadata$SpeciesAccepted))," accepted species, (",length(unique(db1$metadata$SpeciesAuthor)) ," according to authors) and ",
       nrow(db1$metadata), " matrices.\n\n",sep=""))

#File 2
uniqueSource <- unique(paste(db2$metadata$Authors,db2$metadata$Journal,db2$metadata$YearPublication))                      
cat(paste("File-2 contains the demographic and associated data for ", 
      length(uniqueSource), " studies, corresponding to ",
      length(unique(db2$metadata$SpeciesAccepted))," accepted species, (",length(unique(db2$metadata$SpeciesAuthor)) ," according to authors) and ",
      nrow(db2$metadata), " matrices.\n\n",sep=""))

cat("DETAILED SUMMARY\n")
#Accepted species in File 1 that are not in File 2
sp1 <- unique(db1$metadata$SpeciesAccepted)
sp2 <- unique(db2$metadata$SpeciesAccepted)

cat("Number of accepted species in File 1\n")
print(length(sp1))

cat("Number of accepted species in File 2\n")
print(length(sp2))

cat("Accepted species in File 1 that are not in File 2\n")
print(sp1[which(!sp1%in%sp2)])

cat("Accepted species in File 2 that are not in File 1\n")
print(sp2[which(!sp2%in%sp1)])


asp1 <- unique(db1$metadata$SpeciesAuthor)
asp2 <- unique(db2$metadata$SpeciesAuthor)
cat("Number of author species in File 1\n")
print(length(asp1))

cat("Number of author species in File 2\n")
print(length(asp2))

cat("Author species in File 1 that are not in File 2\n")
print(asp1[which(!asp1%in%asp2)])

cat("Author species in File 2 that are not in File 1\n")
print(asp2[which(!asp2%in%asp1)])
}



