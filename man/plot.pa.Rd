\name{plot.pa}
\alias{plot.pa}
\title{Summary statistics of the parameters.}
\description{
This function is used to obtain MCMC summary statistics.
}
\usage{
\method{plot}{pa}(x, treatment=FALSE, differences=FALSE, pval=FALSE, 
        variance=FALSE, covariance=FALSE, by=0.5, ...)
}
\arguments{
\item{x}{Object of class inheriting from "pa".}
\item{treatment}{Logical, if TRUE then plots the treatment surface.}
\item{differences}{Logical, if TRUE then plots the surface of the treatment differences.}
\item{pval}{Logical, if TRUE then plots the surface of the p-value for the treatment differences. H0: No difference in the treatmetns.}
\item{variance}{Logical, if TRUE then plots the surface of the treatment variances.}
\item{covariance}{Logical, if TRUE then plots the surface of the treatment covariances.}
\item{by}{Numeric statement to adjust the labels of the plot.}
\item{...}{Other arguments.}
}
\value{
\item{sig2eps}{Summary statistics for \eqn{\sigma_\epsilon^2}.} 
\item{sig2eta}{Summary statistics for \eqn{\sigma_\eta^2}.} 
\item{phi}{Summary statistics for spatial decay parameter \eqn{\phi}, if estimated using \code{spT.decay}.} 
\item{...}{Summary statistics for other parameters used in the models.} 
}
\seealso{
\code{\link{run}}.
}
\examples{
\dontrun{
##

plot(out, treatment=TRUE) # where out is the output from pa class
plot(out, difference=TRUE) # 
plot(out, pval=TRUE) # 

##
}
}
\keyword{pa}
