\name{USNAM_genoABH}
\alias{USNAM_genoABH}
\docType{data}
\title{
  Reduced genotype data from Maize US-NAM population in ABH format
}
\description{
  Selection of marker and genotypes from the Maize US nested association mapping
  (NAM) population in ABH format.
}
\usage{data(USNAM_genoABH)}
\format{
  \code{data.frame}
}

\details{
  Transformation of the marker scores data \code{\link{USNAM_geno}} in ABH format
  using function \code{\link{cross_ABH}}.
}

\examples{
  data(USNAM_genoABH)
}

\seealso{
\code{\link{cross_ABH}}, \code{\link{USNAM_geno}}
}

\keyword{datasets}