#' Partial Mantel tests on costdistance matrices
#' 
#' This function implements the Causal modelling approach as suggested by
#' Wassermann et al. 2010 and Cushman et al. 2010. It tests for the effect of
#' landscape features using a cost distance matrix on the genetic structure of
#' subpopulation/individuals.
#' 
#' see \code{\link{landgenreport}}
#' 
#' @param gen.mat pairwise genetic distance matrix
#' @param cost.mats pairwise cost distance matrix
#' @param eucl.mat pairwise Eukclidean distance matrix
#' @param plot switch for control plots of the partial mantel test
#' @param nperm number of permutations for the partial mantel test
#' @return A table with the results of the partial mantel test. Using plot=TRUE
#' results in diagnostic plots for the partial mantel tests.
#' @author Bernd Gruber (bernd.gruber@@canberra.edu.au)
#' @seealso \code{\link{popgenreport}}, \code{\link{genleastcost}},
#' \code{\link{landgenreport}}, \code{\link{lgrMMRR}}
#' @references Wassermann, T.N., Cushman, S. A., Schwartz, M. K. and Wallin, D.
#' O. (2010). Spatial scaling and multi-model inference in landscape genetics:
#' Martes americana in northern Idaho. Landscape Ecology, 25(10), 1601-1612.
#' @examples
#' 
#' \donttest{
#' library(raster)
#' fric.raster <- readRDS(system.file("extdata","fric.raster.rdata", package="PopGenReport"))
#' glc <- genleastcost(landgen, fric.raster, "D", NN=8)
#' wassermann(eucl.mat = glc$eucl.mat, cost.mats = glc$cost.mats, gen.mat = glc$gen.mat)
#' }
#' @importFrom lattice xyplot densityplot
#' @export
wassermann <- function(gen.mat, cost.mats, eucl.mat=NULL, plot=TRUE, nperm=999)
{
mats <- cost.mats
if (!is.null(eucl.mat)) mats$Euclidean <- eucl.mat

n.mats <- length(mats)

if (n.mats<2) stop("There are not enough resistance matrices for this approach")
mantel.tab <- data.frame(model=NA, r=NA, p=NA)
cc<- 1


for (i in 1:(n.mats-1))
{

for (j in 2:n.mats)

{

if (i<j)
{
mant1 <- mantel.partial(gen.mat, mats[[i]], mats[[j]],permutations=nperm)
mant2 <- mantel.partial(gen.mat, mats[[j]], mats[[i]],permutations=nperm)


mantel.tab[cc,] <- c(paste("Gen ~",names(mats)[i]," | ", names(mats)[j],sep=""), round(mant1$statistic,4), round(mant1$signif,4))
mantel.tab[cc+1,] <- c(paste("Gen ~",names(mats)[j]," | ", names(mats)[i],sep=""), round(mant2$statistic,4), round(mant2$signif,4))
cc<-cc+2

if (plot==TRUE)
{
x<- densityplot(permustats(mant1),main=paste("Gen ~",names(mats)[i]," | ", names(mats)[j],sep=""), cex.main=0.6)
print(x)

x<- densityplot(permustats(mant2), main=paste("Gen ~",names(mats)[j]," | ", names(mats)[i],sep=""), cex.main=0.6)
print(x)
}

}
}
}

mantel.tab <- mantel.tab[order(mantel.tab$r, decreasing = TRUE),]

return(list(mantel.tab=mantel.tab))
}
