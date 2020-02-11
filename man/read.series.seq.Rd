\name{read.series.seq}
\alias{read.series.seq}
\title{ Reads a series of sequences. }
\description{
 Input from a text file of a series of behavioral (or any other type) sequences. Builds an object of class \code{sequence}
}
\usage{
read.series.seq(fich = NULL)
}
\arguments{
  \item{fich}{ character: name of the text file to read }
}
\details{
The text file must have one line by sequence. Each line must begin by a character identifier of the sequence, followed by a tab character or a space, followed by the series of item codes, separated also by tabs or spaces.
example:\cr\cr
A1 A B A C D A B C \cr
A2 A D D C A B \cr
A3 B D C A D A B C \cr
A4 C B C A \cr
A5 C D A B C \cr
}
\value{
An object of class \code{sequence}
}
\author{ Jean-Sebastien Pierre\cr 
	\email{Jean-sebastien.pierre@univ-rennes1.fr}\cr
}
\seealso{ \code{\link{compseq}} }
\examples{
sink("series-seq")
cat("A1 A B A C D A B C\n")
cat("A2 A D D C A B\n")
cat("A3 B D C A D A B C\n")
cat("A4 C B C A\n")
cat("A5 C D A B C\n")
sink()
read.series.seq("series-seq")
file.remove("series-seq")
}
\keyword{ multivariate }
