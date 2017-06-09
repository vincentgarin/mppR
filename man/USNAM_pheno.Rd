\name{USNAM_pheno}
\alias{USNAM_pheno}
\docType{data}
\title{
  Reduced phenotype data from Maize US-NAM population
}
\description{
  Reduced phenotype data from the Maize US nested association mapping (NAM)
  population (McMullen et al., 2009).
}
\usage{data(USNAM_pheno)}
\format{
  data.frame
}
\details{
  Two columns data frame with genotype indicator and upper leaf angle (ULA)
  trait values. These genotypes correspond to the 500 offspring
  genotypes of the marker matrix (USNAM_geno). The data are available
  on www.panzea.org.
  
}
\source{
  http://www.panzea.org/lit/data_sets.html
}
\references{
  McMullen, M. D., Kresovich, S., Villeda, H. S., Bradbury, P., Li, H., Sun, Q.,
  ... & Buckler, E. S. (2009). Genetic properties of the maize nested
  association mapping population. Science, 325(5941), 737-740.
  
}
\examples{
data(USNAM_pheno)
}

\seealso{
\code{\link{USNAM_geno}}, \code{\link{USNAM_map}}
}

\keyword{datasets}