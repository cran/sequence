\name{homoth}
\alias{homoth}
\title{ Homothety }
\description{
Homothety (scaling) of center \code{centre} and of factor \code{scale} on an object \code{form}
}
\usage{
homoth(form = rbind(c(0, 1), c(1, 1)), centre = c(0, 1), scale = 0.5)
}
\arguments{
  \item{form}{ matrix: with two columns x and y. Set of points making the shape }
  \item{centre}{ vector: Coordinates of the homothety center }
  \item{scale}{ numeric: Homothety ratio }
}
\details{
  	Utilitary used by \code{\link{smartArrow}}
}
\value{
	A matrix, transformed from \code{form} after homothety.
}
\author{ Jean-Sebastien Pierre\cr 
	\email{Jean-sebastien.pierre@univ-rennes1.fr}\cr

}
\seealso{ \code{\link{smartArrow}} }
\examples{
homoth(rbind(c(0,sqrt(2)/2),c(0,sqrt(2)/2)),c(0,0),2)

}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
