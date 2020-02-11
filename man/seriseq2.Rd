\name{seriseq2}
\alias{seriseq2}
\docType{data}
\title{random series}
\description{
  25 random series made of letters "a" to "f".
}
\usage{data(seriseq2)}
\format{
List of 25 lists
\tabular{rlllllll}{
[[1]] \tab chr \tab [1:23] \tab "A1" \tab "d" \tab "d" \tab "d" \tab ...\cr
[[2]] \tab chr \tab [1:31] \tab "A2" \tab "b" \tab "b" \tab "c" \tab ...\cr
[[3]] \tab chr \tab [1:27] \tab "A3" \tab "f" \tab "e" \tab "e" \tab ...\cr
[[4]] \tab chr \tab [1:12] \tab "A4" \tab "a" \tab "c" \tab "a" \tab ...\cr
[[5]] \tab chr \tab [1:9] \tab "A5" \tab "a" \tab "b" \tab "e" \tab ...\cr
[[6]] \tab chr \tab [1:11] \tab "A6" \tab "e" \tab "f" \tab "a" \tab ...\cr
[[7]] \tab chr \tab [1:30] \tab "A7" \tab "f" \tab "d" \tab "f" \tab ...\cr
[[8]] \tab chr \tab [1:29] \tab "A8" \tab "b" \tab "a" \tab "e" \tab ...\cr
[[9]] \tab chr \tab [1:24] \tab "A9" \tab "d" \tab "b" \tab "a" \tab ...\cr
[[10]] \tab chr \tab [1:12] \tab "A10" \tab "a" \tab "f" \tab "b" \tab ...\cr
[[11]] \tab chr \tab [1:18] \tab "A11" \tab "d" \tab "d" \tab "e" \tab ...\cr
[[12]] \tab chr \tab [1:22] \tab "A12" \tab "e" \tab "c" \tab "e" \tab ...\cr
[[13]] \tab chr \tab [1:25] \tab "A13" \tab "e" \tab "f" \tab "f" \tab ...\cr
[[14]] \tab chr \tab [1:20] \tab "A14" \tab "c" \tab "f" \tab "a" \tab ...\cr
[[15]] \tab chr \tab [1:29] \tab "A15" \tab "e" \tab "c" \tab "c" \tab ...\cr
[[16]] \tab chr \tab [1:19] \tab "A16" \tab "c" \tab "a" \tab "f" \tab ...\cr
[[17]] \tab chr \tab [1:25] \tab "A17" \tab "f" \tab "e" \tab "d" \tab ...\cr
[[18]] \tab chr \tab [1:18] \tab "A18" \tab "d" \tab "c" \tab "e" \tab ...\cr
[[19]] \tab chr \tab [1:25] \tab "A19" \tab "c" \tab "a" \tab "b" \tab ...\cr
[[20]] \tab chr \tab [1:8] \tab "A20" \tab "f" \tab "d" \tab "a" \tab ...\cr
[[21]] \tab chr \tab [1:14] \tab "A21" \tab "c" \tab "e" \tab "e" \tab ...\cr
[[22]] \tab chr \tab [1:13] \tab "A22" \tab "f" \tab "c" \tab "e" \tab ...\cr
[[23]] \tab chr \tab [1:28] \tab "A23" \tab "b" \tab "e" \tab "d" \tab ...\cr
[[24]] \tab chr \tab [1:18] \tab "A24" \tab "e" \tab "f" \tab "a" \tab ...\cr
[[25]] \tab chr \tab [1:29] \tab "A25" \tab "e" \tab "a" \tab "b" \tab ...\cr
 }
}
\details{
Example of format used by \code{compseq} and built from a text file by the function \code{read.serie.seq}. \cr
List of lists. Each elementary list represents a particular sequence which must be compared to all others.
}
\source{
   Random simulation made by J.S. Pierre
}
\examples{
data(seriseq2)
compseq(seriseq2)
}
\keyword{datasets}
