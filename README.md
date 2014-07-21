compadreDB
==========

This repository contains code for the COMPADRE Plant Matrix Database.
It will evolve into an R package but will initially contain scripts for interacting with the COMPADRE database which is released as a structured R list object.

The structure of the compadre data object
-----------------------------------------
    compadre |-- metadata {data.frame with 57 columns of metadata and 5621 rows (one for each set of matrices}
             |
             |-- matrixClass {data.frame with 5621 rows: MatrixClassOrganized, MatrixClassAuthor, MatrixClassNumber}
             |
             |-- mat (5621 entries, one for each row of metadata)
                   |-- matA {matrix}
                   |-- matU {matrix}
                   |-- matF {matrix}
                   |-- matC {matrix}



