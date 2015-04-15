\name{spautocor}
\alias{spautocor}
\title{Spatial autocorrelation following Smouse and Pekall 1999}
\description{
Global spatial autocorrelation is a multivariate approach combining all loci into a single analysis. The autocorrelation coefficient r is calculated for each pairwise genetic distance pairs for all specified distance classes. For more information see Smouse and Peakall 1999, Peakall et a. 2003 and Smouse et al. 2008. 


}
\usage{
spautocor(gen.m, eucl.m, shuffle=FALSE, bins=10)
}
\arguments{
  \item{gen.m}{a matrix of individual pairwise genetic distances. Easiest to use gd_smouse or gd_kosman to create such a matrix, but in priniciple any other squared distance matrix can be used. see example}
  \item{eucl.m}{A euclidean distance matrix, based on the coordinates of individuals. see example }
  \item{shuffle}{used internally for the permutation calculation}
  \item{bins}{number of bins for the distance classes. Currently only even bins are supported.
 }
}
\value{
Returns a data frame with r values and number of distances within each distance class.
}
\references{
Smouse PE, Peakall R. 1999. Spatial autocorrelation analysis of individual multiallele and multilocus genetic structure. Heredity 82: 561-573.

Double, MC, et al. 2005. Dispersal, philopatry and infidelity: dissecting local genetic structure in superb fairy-wrens (Malurus cyaneus). Evolution 59, 625-635. 

Peakall, R, et al. 2003. Spatial autocorrelation analysis offers new insights into gene flow in the Australian bush rat, Rattus fuscipes. Evolution 57, 1182-1195. 

Smouse, PE, et al. 2008. A heterogeneity test for fine-scale genetic structure. Molecular Ecology 17, 3389-3400. 

Gonzales, E, et al. 2010. The impact of landscape disturbance on spatial genetic structure in the Guanacaste tree, Enterolobium cyclocarpum(Fabaceae). Journal of Heredity 101, 133-143. 

Beck, N, et al. 2008. Social constraint and an absence of sex-biased dispersal drive fine-scale genetic structure in white-winged choughs. Molecular Ecology 17, 4346-4358. }
\author{
Bernd Gruber, Bernd.Gruber@canberra.edu.au
}
\seealso{
\code{\link{popgenreport}}
}
\examples{
#not run:
#data(bilby)
#popgenreport(bilby, mk.spautocor=TRUE, mk.pdf=FALSE)
#to get a pdf output you need to have a running Latex version installed on your system.
#popgenreport(bilby[1:50], mk.spautocor=TRUE, mk.pdf=TRUE)
}
