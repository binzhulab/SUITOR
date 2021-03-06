% --- Source file: data.Rd ---
\name{data}
\alias{data}
\alias{results}

\title{
  Data for examples
}
\description{
 Example input data and results
}
\details{ 
Contains an example input data object of size 96 by 30 and an example matrix of results.
It is generated by \code{rpois} with mean WH where W (96 by 1) is profile of signature 8 obtained from 
\url{https://cancer.sanger.ac.uk/cosmic/signatures/SBS/SBS8.tt} and 
H (1 by 30) is rounded integer generated by uniform distribution between 20,000 and 40,000.
}
\examples{

 data(data, package="SUITOR")

 # Display a subset of data objects
 data[1:5, 1:5]
 names(results)

}
\keyword{data}
\seealso{ \code{\link{suitor}} }



