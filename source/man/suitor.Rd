
\name{suitor}
\alias{suitor}
\title{suitor}
\description{
   Selecting the number of mutational signatures through cross-validation                          
}
\usage{ suitor(data, op=NULL)
}
\arguments{
  \item{data}{Data frame or matrix containing mutational signatures. 
             This object must contain non-negative values}
  \item{op}{List of options (see details). The default is NULL.}
}
\details{
  The algorithm finds the optimal rank by applying k-fold cross validation.
 
  \bold{Options list op:}
  \tabular{lll}{
    \bold{Name} \tab \bold{Description}                                                        \tab \bold{Default Value} \cr
    em.eps      \tab EM algorithm stopping tolerance                                           \tab 1e-5          \cr
    get.summary \tab 0 or 1 to create summary results                                          \tab 1             \cr  
    k.fold      \tab Number of folds                                                           \tab 10            \cr  
    max.iter    \tab Maximum number of iterations in EM algorithm                              \tab 2000          \cr
    max.rank    \tab Maximum rank                                                              \tab 10            \cr
    min.rank    \tab Minimum rank                                                              \tab 1             \cr
    min.value   \tab Minimum value of matrix before factorizing                                \tab 1e-4          \cr
    n.cores     \tab Number of cores to use for parallel computing                             \tab 1             \cr
    n.seeds     \tab Number of seeds (starting points)                                         \tab 30            \cr
    plot        \tab 0 or 1 to produce an error plot                                           \tab 1             \cr
    print       \tab 0 or 1 to print info                                                      \tab 1             \cr 
    seeds       \tab Vector of seeds (takes precendence over n.seeds)                          \tab NULL          \cr
    type        \tab Socket type in \code{\link[parallel]{makeCluster}} for parallel computing \tab NULL          \cr
    kfold.vec   \tab Vector of values in 1:k.fold when running on a cluster                    \tab NULL
  }

  \bold{Parallel computing} \cr
  If \code{n.cores > 1} and \code{type = NULL}, then \code{type} will be set to "FORK" when running on unix, and
  set to "PSOCK" otherwise. \cr
  \bold{NOTE:} The R_LIBS_USER environment variable may need to be set to the path where the SUITOR R package
       was installed, especially if it was not installed in the default location. Examples of setting the 
       R_LIBS_USER environment variable on unix and windows are below. \cr
       (unix) export R_LIBS_USER=/path/to/SUITOR/package \cr
       (windows) set R_LIBS_USER=/path/to/SUITOR/package \cr

  \bold{Utilizing a cluster} \cr
   When running on a cluster, the option \code{get.summary} should be set to 0.
   For fastest running jobs, set the options \code{min.rank = max.rank},
  \code{kfold.vec} to a single integer in \code{1:k.fold}, 
  and \code{seeds} to a single value.

  
}
\value{ A list containing the objects:
  \itemize{
    \item \code{rank}: The optimal rank
    \item \code{all.results}: Matrix containing training and testing errors for all values
                             of seeds, ranks, folds.
    \item{summary}: Data frame of summarized results for each possible rank created from \code{all.results}.
                    The \code{MSErr} column is defined as
                    sqrt(\{fold1 + ... +foldK\}/\{nrow(data)*ncol(data)\})
  }
}
\references{
}
\author{
Donghyuk Lee <dhyuklee@pusan.ac.kr> and Bin Zhu <bin.zhu@nih.gov>
}
\examples{
  data(SimData, package="SUITOR")

  # Using the default options will take several minutes to run
  #ret <- suitor(SimData)
  
}
\keyword{ NMF, cross-validation, mutational signatures } 
\seealso{ \code{\link{getSummary}}, \code{\link{plotErrors}} }
