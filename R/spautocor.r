#' Spatial autocorrelation following Smouse and Pekall 1999
#' 
#' Global spatial autocorrelation is a multivariate approach combining all loci
#' into a single analysis. The autocorrelation coefficient r is calculated for
#' each pairwise genetic distance pairs for all specified distance classes. For
#' more information see Smouse and Peakall 1999, Peakall et a. 2003 and Smouse
#' et al. 2008.
#' 
#' 
#' @param gen.m a matrix of individual pairwise genetic distances. Easiest to
#' use gd_smouse or gd_kosman to create such a matrix, but in priniciple any
#' other squared distance matrix can be used. see example
#' @param eucl.m A euclidean distance matrix, based on the coordinates of
#' individuals. see example
#' @param shuffle used internally for the permutation calculation
#' @param bins number of bins for the distance classes. Currently only even
#' bins are supported.
#' @return Returns a data frame with r values and number of distances within
#' each distance class.
#' @author Bernd Gruber, Bernd.Gruber@@canberra.edu.au
#' @seealso \code{\link{popgenreport}}
#' @references Smouse PE, Peakall R. 1999. Spatial autocorrelation analysis of
#' individual multiallele and multilocus genetic structure. Heredity 82:
#' 561-573.
#' 
#' Double, MC, et al. 2005. Dispersal, philopatry and infidelity: dissecting
#' local genetic structure in superb fairy-wrens (Malurus cyaneus). Evolution
#' 59, 625-635.
#' 
#' Peakall, R, et al. 2003. Spatial autocorrelation analysis offers new
#' insights into gene flow in the Australian bush rat, Rattus fuscipes.
#' Evolution 57, 1182-1195.
#' 
#' Smouse, PE, et al. 2008. A heterogeneity test for fine-scale genetic
#' structure. Molecular Ecology 17, 3389-3400.
#' 
#' Gonzales, E, et al. 2010. The impact of landscape disturbance on spatial
#' genetic structure in the Guanacaste tree, Enterolobium
#' cyclocarpum(Fabaceae). Journal of Heredity 101, 133-143.
#' 
#' Beck, N, et al. 2008. Social constraint and an absence of sex-biased
#' dispersal drive fine-scale genetic structure in white-winged choughs.
#' Molecular Ecology 17, 4346-4358.
#' @examples
#' 
#' \dontrun{
#' data(bilby)
#' popgenreport(bilby, mk.spautocor=TRUE, mk.pdf=FALSE)
#' #to get a pdf output you need to have a running Latex version installed on your system.
#' #popgenreport(bilby[1:50], mk.spautocor=TRUE, mk.pdf=TRUE)
#' }
#' @export
 

spautocor <- function(gen.m,eucl.m, shuffle=FALSE, bins = 10)
{

gd <- gen.m
ed <- eucl.m


if (shuffle==TRUE)
{

gdd <- as.dist(gd)
 gdsample <- sample(1:length(gdd), length(gdd))
gd[lower.tri(gd)] <- gdd[gdsample]
gd[upper.tri(gd)] <- gdd[gdsample]
diag(gd) <- 0
}




 cdmat <- function(gd)
 {
 dimen <- nrow(gd)
 sgd <- sum(gd, na.rm=TRUE)
 cscd <- matrix(colSums(gd, na.rm=TRUE), dimen, dimen) 
 rscd <- matrix(rowSums(gd, na.rm=TRUE), dimen, dimen, byrow=TRUE) 
 cd <- 0.5*( -gd + 1/dimen*( cscd + rscd) - 1/dimen^2*(sgd ))
 cd 
 }
 cd <- cdmat(gd)
#remove upper triangel to speed things up....
ed[upper.tri(ed)] <-NA
diag(ed) <- NA
r<- NA
distance <- NA
N<- NA

steps <- signif(diff(range(ed, na.rm=TRUE))/bins,4)
for (d in 1:bins )
{
index <- which(ed<=d*steps & ed >(d-1)*steps, arr.ind=TRUE)
cx <- sum(cd[index])
cxii<-sum(diag(cd)[index[,1]])
cxjj<-sum(diag(cd)[index[,2]])
r[d] <-  2 * cx /(cxii+cxjj)

distance[d] <- steps*d
N[d] <- length(index)
}
if (shuffle==FALSE) res <- data.frame(bin = distance,  N=N, r =r)
else res <- data.frame(r=r)

res

}





#b<- redpossums[1:100]
#b@other$xy <- b@other$latlong
#

#xy <- read.csv("D:\\Bernd\\Projects\\aprasia\\apfinal\\apxy.csv")
#
#aprasia@other$xy <- xy
#
#gen.m<-as.matrix(gd_smouse(cats, verbose=FALSE))
#eucl.m <- as.matrix(dist(cats@other$xy))
#reps=1000
#bins=10
#
#splist<- spautocor(gen.m, eucl.m, bins=20)
#
#
#system.time(
#bssplist <- replicate(reps, spautocor(gen.m, eucl.m,shuffle=TRUE, bins=bins))
#)
#
#bs <-matrix(unlist(bssplist), nrow=reps, ncol=bins, byrow=TRUE)
#
#bs.l <- apply(bs,2, quantile, probs=0.025, na.rm=TRUE)
#bs.u <- apply(bs,2, quantile, probs=0.975, na.rm=TRUE)
#
#
#
#matplot(cbind(splist$r,bs.u, bs.l), type="l", lty=c(1,2,2), lwd=c(2,1,1), ylab="Spatial autocorrelation r", axes=FALSE, col=c(1,3,3), xlab="distance")
#axis(2)
#axis(1, at=1:nrow(splist), labels=signif(splist$bin,3))
#axis(1, at=1:nrow(splist), labels=splist$N, line=1, tick=FALSE)
#box()
#mtext("N=",1,line=2, at=0)
#mtext("Bins",1,line=1, at=0)
#
