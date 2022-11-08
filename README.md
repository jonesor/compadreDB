com(p)adreDB
==========

**NOTE:** Most of the functionality that was stored in this repository is now available in the R packages [Rcompadre](https://github.com/jonesor/Rcompadre) and [Rage](https://github.com/jonesor/Rage). The latest version of the database is at www.compadre-db.org and can be downloaded directly into R using `x <- Rcompadre::cdb_fetch("compadre")` or `x <- Rcompadre::cdb_fetch("comadre")`.


------

This repository contains code for the COMPADRE Plant Matrix Database and COMADRE Animal Matrix Database.
It will evolve into an R package but will initially contain scripts for interacting with the COMPADRE and COMADRE databases which are released as structured R list objects.

The data can be downloaded from www.compadre-db.org, or www.comadre-db.org.
Follow news at our blog https://compadredb.wordpress.com, and our Twitter accounts @compadreDB @comadreDB.


The structure of the compadre and comadre data objects
-----------------------------------------
    compadre/comadre +-- metadata {dataframe with ca 60 columns and one row one for each set of matrices (mat A, matU, matF, matC)
             |
             +-- matrixClass {list with one entry for each set of matrices. Each entry is a data frame with 3 columns: MatrixClassOrganized, MatrixClassAuthor, MatrixClassNumber.}
             |
             |-- mat {list with one entry for each row of metadata}
             |     |
             |     +-- matA {matrix}
             |     +-- matU {matrix}
             |     +-- matF {matrix}
             |     \-- matC {matrix}
             |
              \-- version {a vector with version information}


