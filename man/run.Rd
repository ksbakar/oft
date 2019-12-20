\name{run}
\alias{run}
\title{Localised geostatistical implementation for identifying treatment dissimilarities in paddocks}
\description{
This function is used to obtain co-kriging estimations based on a set of local methods. 
}
\usage{
run(data = NULL, run.model = "LF", grid.input = NULL, cov.model = "Exp",
    neighbourhood.points = 50, random.seed=1234, ...)
}
\arguments{
\item{data}{This argument is used to input data. The data set should be pre-processed into 4 columns. First two columns will represent x and y coordinates of the data points. Third column will represent the values, and the fourth column will represent the treatment index or name.}
\item{run.model}{Can take input: 'AP': local adaptive based on points, 'FP': local fixed points.}
\item{grid.input}{Prediction grids/locations.}
\item{cov.model}{Covariance function for the spatial effects. The available options  are "Exp", "Sph" and "Gau". Default option: "Exp".}
\item{neighbourhood.points}{This argument is to input neighbourhood points for the 'FP' algorithm.}
\item{random.seed}{User input of seed number. Default 1234.}
\item{...}{Other optional inputs.}
}
\value{
\item{kriged.output}{A data table with columns related to predictions: (1) tretments (2) differences (3) variances (4) covariances (5) z-values.}
\item{est.parameters}{The cokriging parameters estimated from the localised models.}
\item{est.sample}{Number of samples used for fitting the localised models.}
\item{neighbourhood.points}{The number of neighbourhood points for each treatment.}
\item{grid.coordinates}{A 2 column matrix of prediction grid points.}
\item{model}{Returns the covariance model type.}
\item{grid.size}{Applicable if grid.input=NULL is used.}
\item{shapefile}{Applicable if grid.input is a shapefile.}
\item{type}{Returns the model type as a name.}
\item{comp.time}{Returns the computation time.}
}
\references{

}
%\note{
%}
\seealso{
\link{plot.pa}.
}
\examples{
\dontrun{
######################################################
## 
## Langhorne Creek Data
## ref:
## Bramley, R., Lanyon, D., and Panten, K. (2005). Whole-of-vineyard experimentation - 
## An improved basis for knowledge generation and decision making. 
## In Proceedings of the 5th European Conference on Precision Agriculture, p.883-890.
##

library(sp); library(fields); library(spTimer); library(ggplot2)
data(LangHorneCreek)
dim(LangHorneCreek)
head(LangHorneCreek)
grid <- data.frame(spT.grid.coords(c(range(LangHorneCreek[,1])),c(range(LangHorneCreek[,2])),by=c(100,100)))
p <- ggplot() + 
  geom_point(data=grid,aes(x=X1,y=X2),pch="+") +
  geom_point(data=LangHorneCreek,aes(x=east,y=north,color=treatment)) +
  labs(x="Easting",y="Northing",color="Treatment")
p

##
## fixed neighbourhood algorithm 
##

out.lf <- run(data=LangHorneCreek,run.model="LF",grid.input=grid)
names(out.lf)
plot(out.lf, treatment=TRUE)
plot(out.lf, differences=TRUE)
plot(out.lf, zval=TRUE)
plot(out.lf, variances=TRUE)
plot(out.lf, covariances=TRUE)

##
## adaptive points algorithm
##

out.p <- run(data=LangHorneCreek,run.model="LA-P",grid.input=grid)
plot(out.p, zval=TRUE)

# estimated samples 
tr <- c("control","mulch","ripped")
par(mfrow=c(1,3))
for(i in 1:3){
dat <- data.frame(sample=round(out.p$est.sample[,i],2),out.p$grid.coordinates)
dat <- dat[,c(2,3,1)]
head(dat)
library(raster); library(viridis)
r <- rasterFromXYZ(dat)
r.range <- c(floor(min(unlist(c(dat[,3])))),ceiling(max(unlist(c(dat[,3])))))
plot(r,col=rev(viridis(100)),breaks=seq(from=r.range[1],to=r.range[2],length.out = 100),
     main=paste0("Estimated sample for ",tr[i]),
     axis.args=list(at=seq(r.range[1], r.range[2], length.out=10),
                    labels=round(seq(r.range[1], r.range[2],length.out=10)),
                    cex.axis=2))
contour(r, add=TRUE, lty=2)
}
par(mfrow=c(1,1))

##
######################################################
}
}
\keyword{paTools}