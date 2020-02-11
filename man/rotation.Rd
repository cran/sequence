\name{rotation}
\alias{rotation}
\title{Plane rotation }
\description{
  Rotation of the object \code{form} around a center \code{center}, of angle \code{angle}. Angle in radians in the direct (trigonometric) direction.
}
\usage{
rotation(form, center, angle)
}

\arguments{
  \item{form}{ matrix: two columns, x and y, giving the coordinates of the points composing the object.}
  \item{center}{ vector: coordinates of the rotation center. }
  \item{angle}{ numeric: rotation angle in radians }
}
\details{
	A geometric function used by \code{\link{smartArrow}}
}
\value{
    matrix: \code{form} object rotated.
}
\author{ Jean-Sebastien Pierre \email{Jean-sebastien.pierre@univ-rennes1.fr}\cr
}
\seealso{ \code{\link{smartArrow}} }
\examples{
rotation(cbind(c(0,1),c(0,0)),c(0,0),pi/3)
}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
