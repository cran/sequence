\name{symet}
\alias{symet}
\title{ Symmetrization and diagonal loading of a transition matrix }
\description{
  Adds the matrix to its transpose (X + t(X)) and loads its diagonal by the sum of each row.
}
\usage{
symet(x, sym = TRUE, charge = TRUE, ident = FALSE)
}
\arguments{
  \item{x}{ data.frame: The matrix to process as a data.frame }
  \item{sym}{ logical: symmetrisation.}
  \item{charge}{ logical: diagonal loading.}
  \item{ident}{ logical: first column consists in identifiers.}
}
\value{
	a data.frame: the resulting matrix
}
\references{ 
Pierre, J. S. and C. Kasper (1990). The Design of Ethological Flow-Charts on Factorial Analysis Representations - an Application to the Study of the male Mole-Cricket Sexual Courtship. Biology of Behaviour 15(3-4): 125-151.\cr
Van der Heijden, P. G. M. 1986. Transition matrices, model fitting and correspondence analysis. In: Data Analysis and Informatics IV (Ed. by E. Diday), pp. 221-226. Elsevier Science Publishers.
}
\author{ Jean-Sebastien Pierre\cr 
	\email{Jean-sebastien.pierre@univ-rennes1.fr}\cr
}
\seealso{ \code{\link{flux}},\code{\link{matrcom}} }
\examples{
	data(wiepkema)
	swiepkema<-symet(wiepkema,ident=TRUE)
	print(swiepkema)
}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
