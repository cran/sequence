\name{read.old.seq}
\alias{read.old.seq}
\title{ Function to read the deprecated files *.SEQ }
\description{
 Used to recover the sequence files used in a former DOS version of the program.
}
\usage{
read.old.seq(file = "NULL")
}
\arguments{
  \item{file}{ Fichier texte a lire }
}
\details{
Only used for compatibility with the old format *.SEQ of program \code{graflux}: Items coded by at most 4 characters, One code per line. The separator of sequences is a code of 4 characters beginning by  \%.\cr 
Example:\cr
 
\%001 \cr
AAAA \cr
BBBB \cr
EEEE \cr
DDDD \cr
AAAA \cr
FFFF \cr
CCCC \cr
FFFF \cr
\%002 \cr
FFFF \cr
BBBB \cr
CCCC \cr
.... \cr
}

\value{
An object de classe \code{sequence}, i.e. a list of lists.
}
\references{ 
Pierre, J. S. and C. Kasper (1990). The Design of Ethological Flow-Charts on Factorial Analysis Representations - an Application to the Study of the Male Mole-Cricket Sexual Courtship. Biology of Behaviour 15(3-4): 125-151.(in french) \cr
Van der Heijden, P. G. M. 1986. Transition matrices, model fitting and correspondence analysis. 
In: Data Analysis and Informatics IV (Ed. by E. Diday), pp. 221-226. Elsevier Science Publishers.
}
\author{ Jean-Sebastien Pierre \email{Jean-sebastien.pierre\@univ-rennes1.fr}\cr
}
\seealso{ \code{\link{read.series.seq}} for a more convenient format, \code{\link{compseq}} for the analysis}
\examples{
sink("old-seq")
cat("\%001\\n")
cat("aaa\\n")
cat("bbb\\n")
cat("aaa\\n")
cat("ccc\\n")
cat("\%002\\n")
cat("ccc\\n")
cat("ccc\\n")
cat("aaa\\n")
cat("aaa\\n")
cat("bbb\\n")
sink()
seq1<-read.old.seq("old-seq")
seq1
file.remove("old-seq")
}
\keyword{ multivariate }

