\name{geneseq}
\alias{geneseq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{  Generation aof random sequences }
\description{
  Generates \code{nseq} random sequences of items represented by codes, either by random sampling (order 0), either using a Discrete Time Markov Chain (DTMC) transition matrix (order 1).
}
\usage{
geneseq(nseq, lmin, lmax, order = 0, dico = NULL, mattrans = NULL)
}
\arguments{
  \item{nseq}{ numeric: number of sequences to generate. }
  \item{lmin}{ numeric: minimum length of sequences }
  \item{lmax}{ numeric: maximum length of sequences }
  \item{order}{ numeric: order of the Markov Process, 0 or 1. See details. }
  \item{dico}{ character vector: dictionnary of items. }
  \item{mattrans}{ data.frame: transition matrix with identifiers as first column. }
}
\details{
Generates random sequences of states according to a Markov process of order 0 or 1. \cr 
order 0: the following state is drawn with equal probabilities in the item list, no influence of previous state.\cr
order 1: the following state is drawn with a probability proportional to the product of the transition matrix by the vector of states.\cr
\code{nseq} sequences of random size are generated. The size of each sequence is a random integer drawn under uniform probability between \code{lmin} and \code{lmax} (both included). 
}
\value{
A list of nseq items
  \item{A1 }{list: first sequence}
  \item{A2 }{list: second sequence}
  ... and so on until \code{nseq}.
}
\author{ Jean-Sebastien Pierre\cr 
	\email{Jean-sebastien.pierre@univ-rennes1.fr}\cr
}
\note{ The transition matrix must be square and Markovian. A non-square matrix produces an error and program break. Its columns must all sum to 1. No verification is done for this last point. A matrix in which the columns do not sum to 1 may lead to spurious or unpredictible results.
}
\seealso{ \code{\link{compseq}} }
\examples{
# 1 - generates 5 random sequences of length from 1 to 15 with  
# symbols a,b,c,d.(order = 0)
geneseq(5,1,15,dico=letters[1:4]) 
mat<-data.frame(id=letters[1:3],a=c(0,0.75,0.25),b=c(0.5,0.25,0.25),
c=c(0.1,0.8,0.1))
# 2 - generates 2 random sequences of length 10 with symbols a,b,c,d  
# according to the transition matrix mat
geneseq(2,10,10,order=1,mattrans=mat) 
}
\keyword{ multivariate }
