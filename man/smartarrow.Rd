\name{smartArrow}
\alias{smartArrow}
\title{ Draws a 'smart' arrow between two points. }
\description{
  Draws a 'smart' arrow between two points surrounded by a circle of known radius. The arrow is arc shaped and ends with an arrowhead.
}
\usage{
smartArrow(A = c(0, 2), B = c(2, 2), Ra = 0.2, Rb = 0.1, ClegA = "A", 
ClegB = "B", width = 0.1, col = "lightgreen", ccol = "yellow", 
density = NULL, angle = 0, marge = 1.3, debord = 1.4, plot = FALSE, 
trace = TRUE)
}
\arguments{
  \item{A}{numeric vector: coordinates of point A. }
  \item{B}{numeric vector: coordinatess of point B. }
  \item{Ra}{numeric: Radius of circle A }
  \item{Rb}{numeric: Radius of circle B } 
  \item{ClegA}{character: label of circle A}
  \item{ClegB}{character: label of circle B }
  \item{width}{numeric: width of the arrow }
  \item{col}{color to fill the arrow. See\code{\link{polygon}} }
  \item{ccol}{color to fill the circle. See \code{\link{symbols}} }
  \item{density}{numeric: density of shading lines per inch (default=NULL). See \code{polygon}}
  \item{angle}{numeric: angle of shading lines in radians (default = 0).  See \code{polygon}}
  \item{marge}{numeric: Angular distance between each circle and the corresponding end of the arrow (default = 1.3).}
  \item{debord}{numeric: overflow of the arrowhead relatively to its body. (Default = 1.4) }
  \item{plot}{logical: to draw a first arrow on a new plot.(Default = FALSE) }
  \item{trace}{logical: if TRUE (default) the arrow is drawn, else it is simply calculated. }
}
\details{
 SmartArrow can work both as a high level and a low level function, depending on the value of the argument \code{plot}.\cr Used in \code{\link{flux}}. Not supposed to be used directly by the user.
}
\value{
  An object of class \code{smartarrow} having attributes:
\item{arrow }{ The polygon representing the arrow.}
\item{bg }{ Color of background}
\item{A }{ Coordinates of point A}
\item{B }{ Coordinates of point B}
\item{Ra }{ Radius of circle A}
\item{Rb }{ Radius of circle B}
\item{left }{The most extreme left point of the set circles-arrows.}
\item{right }{The most extreme right point of this set}
\item{up}{ The most extreme high point. }
\item{down}{ The most extreme low point }
\item{drawcircles}{ boolean : draw the circles or not.}
}
\author{ Jean-Sebastien Pierre\cr 
	\email{Jean-sebastien.pierre@univ-rennes1.fr}\cr
}
\seealso{ \code{\link{flux}}, \code{\link{polygon}},\code{\link{symbols}} \cr \code{\link{arrows}} ~ for simple arrows}
\examples{


# For a quick illustration
   smartArrow(plot=TRUE)
# an arrow from a randomly drawn point on [0,2]x[0,2]
   smartArrow(B=runif(2,0,2),Rb=runif(1,0,0.5),ClegB="C")
}
\keyword{ multivariate }
