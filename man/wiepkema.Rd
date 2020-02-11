\name{wiepkema}
\alias{wiepkema}
\docType{data}
\title{ Wiepkema's (1961) historical data}
\description{
  Succession matrix for the sexual courtship of the gorchak fish (\emph{Rhodeus sericeus}, cyprinidae)
}
\usage{data(wiepkema)}
\format{
  A data frame with 12 observations on the following 13 variables.
  \describe{
    \item{\code{id}}{factor: levels \code{chf} \code{chs} \code{ffl} \code{fl} \code{hb} \code{hdp} \code{jk} \code{le} \code{qu} \code{sk} \code{sn} \code{tu}}
    \item{\code{jk}}{numeric: jerking}
    \item{\code{tu}}{numeric: turning}
    \item{\code{hb}}{numeric: head butting}
    \item{\code{chs}}{numeric:  chasing}
    \item{\code{fl}}{numeric: fleeing}
    \item{\code{qu}}{numeric: quivering}
    \item{\code{le}}{numeric: leading}
    \item{\code{hdp}}{numeric: head-down}
    \item{\code{sk}}{numerics: kimming}
    \item{\code{sn}}{numeric: snapping}
    \item{\code{chf}}{numeric: chafing}
    \item{\code{ffl}}{numeric: finflickering}
  }
}
\details{
  Historical data of Wipkema (1961) consisting in the succession of patterns in the sexual courtship of a fish. Treated by this author by Principal Component Analysis.
}
\source{
 Wiepkema, P. R. 1961. An ethological analysis of the reproductive behaviour of the bitterling. Arch Neerl Zool, 14, 103-199.
}
\references{
  Van der Heijden, P. G. M. (1990). Correspondance analysis of transition matrices, with special attention to missing entries and asymmetry. Animal Behaviour, 39, 49-64.
}
\examples{
data(wiepkema)
# Symmetrization and diagonal loading before a CA
swiepkema <- symet(wiepkema,ident=TRUE)
ca(swiepkema)
}
\keyword{datasets}
