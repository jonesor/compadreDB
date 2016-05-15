\name{convert2flat}
\alias{convert2flat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to convert from the list structured database object to a flat sheet.
}
\description{
The function converts from the list structured database object to a flat sheet by converting each matrix to a string.
}
\usage{
convert2flat(db = compadre, Aonly = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{db}{
The COMPADRE or COMADRE matrix database object
}
  \item{Aonly}{
A logical value (TRUE/FALSE) indicatting whether ONLY the A matrix be included in the output.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A `data.frame` with the same columns as are present in the metadata part of the COMPADRE/COMADRE object, followed by the string-form matrix stage information and the matrices themselves.
}
\references{
Salguero-Gómez, R., Jones, O. R., Archer, C. R., Buckley, Y. M., Che-Castaldo, J., Caswell, H., Hodgson, D., Scheuerlein, A., Conde, D. A., Brinks, E., de Buhr, H., Farack, C., Gottschalk, F., Hartmann, A., Henning, A., Hoppe, G., Römer, G., Runge, J., Ruoff, T., Wille, J., Zeh, S., Davison, R., Vieregg, D., Baudisch, A., Altwegg, R., Colchero, F., Dong, M., de Kroon, H., Lebreton, J.-D., Metcalf, C. J. E., Neel, M. M., Parker, I. M., Takada, T., Valverde, T., Vélez-Espino, L. A., Wardle, G. M., Franco, M., Vaupel, J. W. (2015), The compadre Plant Matrix Database: an open online repository for plant demography. Journal of Ecology, 103: 202–218. doi: 10.1111/1365-2745.12334

Salguero-Gómez, R., Jones, O. R., Archer, C. R., Bein, C., de Buhr, H., Farack, C., Gottschalk, F., Hartmann, A., Henning, A., Hoppe, G., Römer, G., Ruoff, T., Sommer, V., Wille, J., Voigt, J., Zeh, S., Vieregg, D., Buckley, Y. M., Che-Castaldo, J., Hodgson, D., Scheuerlein, A., Caswell, H., Vaupel, J. W. (2016), COMADRE: a global data base of animal demography. Journal of Animal Ecology, 85: 371–384. doi: 10.1111/1365-2656.12482
}
\author{
Owen R. Jones <jones@biology.sdu.dk>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

#newDB<-convert2flat(compadre,Aonly=FALSE)
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
