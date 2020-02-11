\name{matrcom}
\alias{matrcom}
\title{ Calculates a common succession matrix }
\description{
  Computes the succession matrix common to a set of sequences.
}
\usage{
matrcom(x, printdata = FALSE, printdico = FALSE, printmat = FALSE, printcom = TRUE)
}
\arguments{
  \item{x}{ seq: list of lists  of class \code{character}. The sequences to analyse. }
  \item{printdata}{ logical: print all sequences.}
  \item{printdico}{ logical: print the dictionary.}
  \item{printmat}{ logical: print each individual matrix.}
  \item{printcom}{ logical: print the common succession matrix.}
}
\details{
  The resulting matrix (succession: matrix) gives the number of times the element i (row) is followed by the element j (column).
}
\value{
	data.frame: the succession matrix with the label of items as first column.
}
\references{ 
Pierre, J. S. and C. Kasper (1990). The Design of Ethological Flow-Charts on Factorial Analysis Representations - an Application to the Study of the Male Mole-Cricket Sexual Courtship. Biology of Behaviour 15(3-4): 125-151.(in french) \cr
Van der Heijden, P. G. M. 1986. Transition matrices, model fitting and correspondence analysis. In: Data Analysis and Informatics IV (Ed. by E. Diday), pp. 221-226. Elsevier Science Publishers.
}
\author{ Jean-Sebastien Pierre \email{Jean-sebastien.pierre@univ-rennes1.fr}\cr
}
\seealso{ \code{\link{flux}}, \code{\link{ca}} }
\examples{
data(seriseq)
matrcom(seriseq)

}
\keyword{ multivariate }
