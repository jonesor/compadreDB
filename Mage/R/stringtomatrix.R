#' A function to convert a matrix from the character string format to a matrix
#' 
#' The function converts a matrix from the character string format to class matrix.
#' 
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param A A matrix presented in the form of a string, with values separated by spaces, and bounded by square brackets.
#' @return A `matrix`.
#' @author Owen R. Jones <jones@@biology.sdu.dk>
#' @seealso convert2flat
#' @examples
#' 
#' x <- "[0.46145 0.153 0.33615 NA 0.5238 0.59225 NA 0.2118 0.33615 0.24765 0.5238 0.59225 0.0705 0.0633 0.10815 0.0827 0 0 0.0267 0.0947 0.03155 0.15405 0 0 0.06295 0.19195 0.3909 0.24355 0 0 0.02305 0.281 0.1321 0.35145 0 0]"
#' stringtomatrix(x)
#' 
#' @export stringtomatrix

stringtomatrix <- function(A) {
  A <- gsub(pattern = "\\[|\\]", "", A)
  A <- gsub(pattern = ";", " ", A)
  A <- strsplit(x = A, split = " |`NA`")[[1]]
  matrix(as.numeric(A), nrow = sqrt(length(A)), byrow = TRUE)
}