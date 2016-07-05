#' @export

dbCompare <- function(db1, db2, verbose = FALSE){ 

#Quick summary
cat("Quick Summary\n")

#File 1
uniqueSource1 <- unique(paste(db1$metadata$Authors," (",db1$metadata$YearPublication,") ",db1$metadata$Journal,sep=""))                      
db1$metadata$binomial <- paste(db1$metadata$GenusAccepted,db1$metadata$SpeciesEpithetAccepted,sep = " ")

cat(paste("File-1 contains the demographic and associated data from ", 
       length(uniqueSource1), " source papers, corresponding to ",
       length(unique(db1$metadata$binomial))," accepted species, and ",
       nrow(db1$metadata), " matrices.\n\n",sep=""))

#File 2
uniqueSource2 <- unique(paste(db2$metadata$Authors," (",db2$metadata$YearPublication,") ",db2$metadata$Journal,sep=""))                      
db2$metadata$binomial <- paste(db2$metadata$GenusAccepted,db2$metadata$SpeciesEpithetAccepted,sep = " ")

cat(paste("File-2 contains the demographic and associated data for ", 
      length(uniqueSource2), " source papers, corresponding to ",
      length(unique(db2$metadata$binomial))," accepted species, and ",
      nrow(db2$metadata), " matrices.\n\n",sep=""))

if(verbose == TRUE){
cat("Detailed summary\n")

#Accepted species in File 1 that are not in File 2
sp1 <- unique(db1$metadata$binomial)
sp2 <- unique(db2$metadata$binomial)

cat("Number of accepted species in File 1, based on latin binomial\n")
print(length(sp1))

cat("Number of accepted species in File 2, based on latin binomial\n")
print(length(sp2))

cat("Accepted species in File 1 that are not in File 2 (based on latin binomial)\n")
print(sp1[which(!sp1%in%sp2)])

cat("Accepted species in File 2 that are not in File 1 (based on latin binomial)\n")
print(sp2[which(!sp2%in%sp1)])

#Get unique author species for both files
asp1 <- unique(db1$metadata$SpeciesAuthor)
asp2 <- unique(db2$metadata$SpeciesAuthor)

cat("Number of unique study-species combinations in File 1\n")
print(length(asp1))

cat("Number of unique study-species combinations in File 2\n")
print(length(asp2))

#cat("Study-species in File 1 that are not in File 2\n")
#print(asp1[which(!asp1%in%asp2)])

#cat("Study-species in File 2 that are not in File 1\n")
#print(asp2[which(!asp2%in%asp1)])

cat("\n\nSource papers in File 2 that are not in File 1\n")
print(sort(uniqueSource2[which(!uniqueSource2%in%uniqueSource1)]))

cat("\n\nSource papers in File 1 that are not in File 2\n")
print(sort(uniqueSource1[which(!uniqueSource1%in%uniqueSource2)]))


cat("See the User Guide for definitiions\n")
}
}
