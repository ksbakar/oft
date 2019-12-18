\name{run}
\alias{run}
\title{Localised geostatistical implementation for identifying treatment dissimilarities in paddocks}
\description{
This function is used to obtain co-kriging estimations based on a set of local methods. 
}
\usage{
run(data = NULL, run.model = "LF", grid.input = NULL, cov.model = "Exp",
  neighbourhood.size = NULL, neighbourhood.points = 50, pts.for.each.treat = 10,
  subsample.percentage = 10, num.of.replications = 100, tol.site = 1000, 
  random.seed=1234, ...)
}
\arguments{
\item{data}{This argument is used to input data. The data set should be pre-processed into 4 columns. First two columns will represent x and y coordinates of the data points. Third column will represent the values, and the fourth column will represent the treatment index or name.}
\item{run.model}{Can take input: 'LA-R-A': local adaptive based on radius expanding algorithm,'LA-P': local adaptive based on points, 'LF': local fixed, 'SK': subsampling Kmeans clustering.}
\item{grid.input}{Prediction grids/locations.}
\item{cov.model}{Covariance function for the spatial effects. The available options  are "Exp", "Sph" and "Gau". Default option: "Exp".}
\item{neighbourhood.size}{This argument is to input initial neighbourhood radius for algorithm 'LA-R-A'. If NULL then the algorithm will select the initial radius based on the data.}
\item{neighbourhood.points}{This argument is to input neighbourhood points for the 'LF' algorithm.}
\item{pts.for.each.treat}{This argument is to ensure the number of points for each treatment, default is 10.}
\item{subsample.percentage}{This argument is for 'SK' algorithm, for user defined subsample percentages. Default value is 10.}
\item{num.of.replications}{This argument is for 'SK' algorithm, for user defined replications. Default value is 100 replications. }
\item{tol.site}{Minimum number of points used for global model. Default size = 1000.}
\item{random.seed}{User input of seed number. Default 1234.}
\item{...}{Other optional inputs.}
}
\value{
\item{T}{Total points of time, e.g., days within each year.}
}
\references{


}
%\note{
%}
\seealso{
\link{plot.pa}.
}
\examples{
\dontshow{
######################################################

## 


## Large density points example: Langhorne Creek
## ref: 
library(sp); library(fields); library(spTimer)
data(loc_LanghorneCreek)
dim(loc_LanghorneCreek)
loc <- loc_LanghorneCreek[sort(sample(round(nrow(loc_LanghorneCreek)*0.5))),]
grid <- data.frame(spT.grid.coords(c(range(loc_LanghorneCreek[,1])),c(range(loc_LanghorneCreek[,2])),by=c(150,150)))
data_fnc_LanghorneCreek <- function(loc,grid){
  coordinates(grid) <- ~X1+X2
  bdr <- loc[chull(loc),]
  coordinates(bdr) <- ~Easting+Northing
  bdr <- SpatialPolygons(list(Polygons(list(Polygon(bdr)),1)))
  grid <- as.data.frame(grid[!is.na(over(grid,as(bdr,”SpatialPolygons”))),])
  plot(grid,pch="+"); points(loc,pch=19,col=2)
  rmvn <- function(n, mu=0, V = matrix(1)){
    # from spBayes package 
    p <- length(mu)
    if(any(is.na(match(dim(V),p))))
      stop("Dimension problem!")
    D <- chol(V)
    t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
  }  
  D <- as.matrix(rdist(loc))
  K <- exp(-0.01*D)
  set.seed(1234)
  loc$value <- as.vector(rmvn(1, rep(0,nrow(K)), 0.05*K))
  loc$treatment <- round(runif(nrow(loc),min=1,max=3))
  loc[loc$treatment==1,]$value <- 3+loc[loc$treatment==1,]$value
  loc[loc$treatment==2,]$value <- 2+loc[loc$treatment==2,]$value
  loc[loc$treatment==3,]$value <- 5+loc[loc$treatment==3,]$value
  loc
}
dat <- data_fnc_LanghorneCreek(loc,grid)
head(loc_LanghorneCreek)




######################################################
}
}
\keyword{paTools}