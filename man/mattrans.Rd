\name{mattrans}
\alias{mattrans}
\title{ Creates of a succession matrix }
\description{
  Creates a succession matrix of order 1 from an object of class \code{sequence}.
}
\usage{
mattrans(x, dico = NULL, print = FALSE)
}
\arguments{
  \item{x}{ Cseq: sequence }
  \item{dico}{ character: list of items codes }
  \item{print}{ logical: print the matrix. Default = FALSE }
}
\details{
 Creates a succession, not a transition matrix. Integer numbers (frequencies) and not conditional probabilities. Function used by \code{compseq}
}
\value{
  a matrix: the succession matrix. 
}
\author{ Jean-Sebastien Pierre \email{Jean-sebastien.pierre@univ-rennes1.fr}\cr
}
\note{ Not supposed to be used alone; if needed, the user is invited to use rather \code{\link{matrcom}}
}
\seealso{  \code{\link{compseq}} ,\code{\link{matrcom}}}
\examples{
data(seriseq)
mattrans(unlist(seriseq))

}
\keyword{ multivariate }
