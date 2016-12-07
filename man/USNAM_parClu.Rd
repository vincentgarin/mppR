\name{USNAM_parClu}
\alias{USNAM_parClu}
\docType{data}
\title{
  Parental clustering object
}
\description{
  Example of parental clustering object for the selection from the Maize US
  nested association mapping (NAM) population.
}
\usage{data(USNAM_parClu)}

\details{
  Parent clustering matrix specify at each genome position the results
  of a parent clustering into ancestral groups. The matrix rows
  represent the position and the columns correspond to each parent. For example,
  if we have at the ith row (1, 2, 3, 2, 1), this means that parents 1
  and 5 are in the same group, that 2 and 4 are in another one and that the
  third parent was assigned to any group. This kind of object is required to
  compute an ancestral model and must be provided in the \code{par.clu} argument.
  
}

\examples{
  data(USNAM_parClu)
}

\seealso{
\code{\link{parent_cluster}}
}

\keyword{datasets}