compadreDB
==========

This repository contains code for the COMPADRE Plant Matrix Database.
It will evolve into an R package but will initially contain scripts for interacting with the COMPADRE database which is released as a structured R list object.

The data can be downloaded from www.compadre-db.org.
Follow news at our blog https://compadredb.wordpress.com, and our Twitter account @compadreDB.


The structure of the compadre data object
-----------------------------------------
    compadre +-- metadata {dataframe with 57 columns and one row one for each set of matrices (AUFC)
             |
             +-- matrixClass {list with one entry for each set of matrices. Each entry is a data frame with 3 columns: MatrixClassOrganized, MatrixClassAuthor, MatrixClassNumber.}
             |
             \-- mat {list with one entries each row of metadata}
                   |
                   +-- matA {matrix}
                   +-- matU {matrix}
                   +-- matF {matrix}
                   \-- matC {matrix}



