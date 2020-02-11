\name{seriseq}
\alias{seriseq}
\docType{data}
\title{Random series }
\description{
 5 random series made of letters "a" to "f"
}
\usage{data(seriseq)}
\format{
  The format is:
List of 5 lists
\tabular{rlllllll}{
[[1]] \tab chr \tab [1:9] \tab "A1" \tab "d" \tab "e" \tab "a" \tab ...\cr
[[2]] \tab chr \tab [1:8] \tab "A2" \tab "f" \tab "d" \tab "b" \tab ...\cr
[[3]] \tab chr \tab [1:9] \tab "A3" \tab "c" \tab "e" \tab "b" \tab ...\cr
[[4]] \tab chr \tab [1:8] \tab "A4" \tab "e" \tab "d" \tab "b" \tab ...\cr
[[5]] \tab chr \tab [1:8] \tab "A5" \tab "d" \tab "b" \tab "f" \tab ...\cr
 }
}
\details{
  Example of format used by \code{compseq} and built from a text file by the function \code{read.series.seq}. List of lists. Each elementary list represents a particular sequence which must be compared to all others.
}
\seealso{
\code{\link{read.series.seq}}
}
\author{ Jean-Sebastien Pierre \email{Jean-sebastien.pierre@univ-rennes1.fr}}
\source{
  Random simulation made by J.S. Pierre.
}
\examples{
data(seriseq)
compseq(seriseq)
}
\keyword{datasets}
