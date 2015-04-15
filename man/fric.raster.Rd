\name{fric.raster}
\alias{fric.raster}
\docType{data}
\title{
a simulated friction map
}
\description{
This data set is used to demonstrate the use of the \code{\link{landgenreport}} and \code{\link{genleastcost}}  functions. It is a simple raster file, with dimensions of 50 x 50 and now coordinate system.
}
\usage{data(fric.raster)}
\seealso{
 \code{\link{landgenreport}}, \code{\link{genleastcost}}
}
\examples{
require(raster)
data(fric.raster)
plot(fric.raster)
}
\keyword{datasets}
