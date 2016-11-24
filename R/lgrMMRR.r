#' Multiple Matrix Regression with Randomization analysis
#' 
#' performs Multiple Matrix Regression with Randomization analysis This method
#' was implemented by Wang 2013 (MMRR function see references) and also by
#' Sarah Goslee in package ecodist. lgrMMRR is a simple wrapper to have a more
#' user friendly output.
#' 
#' Performs multiple regression on distance matrices following the methods
#' outlined in Legendre et al. 1994 and implemented by Wang 2013.
#' 
#' @param gen.mat a genetic distance matrix (e.g. output from
#' \code{\link{genleastcost}}
#' @param cost.mats a list of cost distance matrices
#' @param eucl.mat pairwise Euclidean distance matrix. If not specificed
#' ignored
#' @param nperm the number of permutations
#' @return a table with the results of the matrix regression analysis.
#' (regression coefficients and associated p-values from the permutation test
#' (using the pseudo-t of Legendre et al. 1994).  and also r.squared from and
#' associated p-value from the permutation test. F.test.
#' 
#' Finally also the F-statistic and p-value for overall F-test for lack of fit.
#' @author Bernd Gruber (bernd.gruber@@canberra.edu.au) using the
#' implementation of Wang 2013.
#' @seealso \link[ecodist]{MRM} in package ecodist, \code{\link{popgenreport}},
#' \code{\link{genleastcost}}, \code{\link{landgenreport}},
#' \code{\link{wassermann}}
#' @references Legendre, P.; Lapointe, F. and Casgrain, P. 1994. Modeling brain
#' evolution from behavior: A permutational regression approach. Evolution 48:
#' 1487-1499.
#' 
#' Lichstein, J. 2007. Multiple regression on distance matrices: A multivariate
#' spatial analysis tool. Plant Ecology 188: 117-131.
#' 
#' Wang,I 2013. Examining the full effects of landscape heterogeneity on
#' spatial genetic variation: a multiple matrix regression approach for
#' quantifying geographic and ecological isolation. Evolution: 67-12:
#' 3403-3411.
#' @examples
#' 
#' \dontrun{%
#' require(raster)
#' data(landgen)
#' data(fric.raster)
#' glc <- genleastcost(landgen, fric.raster, "D", NN=4, path="leastcost")
#' lgrMMRR(glc$gen.mat, glc$cost.mats, glc$eucl.mat, nperm=999)
#' }
#' @export
lgrMMRR <- function(gen.mat, cost.mats, eucl.mat=NULL, nperm=999)
{

##### MMRR function of Wang
# MMRR performs Multiple Matrix Regression with Randomization analysis
# Y is a dependent distance matrix
# X is a list of independent distance matrices (with optional names)
MMRR<-function(Y,X,nperm=999){
	#compute regression coefficients and test statistics
	nrowsY<-nrow(Y)
	y<-unfold(Y)
	if(is.null(names(X)))names(X)<-paste("X",1:length(X),sep="")
        Xmats<-sapply(X,unfold)   
        fit<-lm(y~Xmats)
	coeffs<-fit$coefficients
	summ<-summary(fit)
	r.squared<-summ$r.squared
	tstat<-summ$coefficients[,"t value"]
	Fstat<-summ$fstatistic[1]
	tprob<-rep(1,length(tstat))
	Fprob<-1

	#perform permutations
	for(i in 1:nperm){
		rand<-sample(1:nrowsY)
		Yperm<-Y[rand,rand]
		yperm<-unfold(Yperm)
		fit<-lm(yperm~Xmats)
		summ<-summary(fit)
                Fprob<-Fprob+as.numeric(summ$fstatistic[1]>=Fstat)
                tprob<-tprob+as.numeric(abs(summ$coefficients[,"t value"])>=abs(tstat))
	}

	#return values
	tp<-tprob/(nperm+1)
	Fp<-Fprob/(nperm+1)
	names(r.squared)<-"r.squared"
	names(coeffs)<-c("Intercept",names(X))
	names(tstat)<-paste(c("Intercept",names(X)),"(t)",sep="")
	names(tp)<-paste(c("Intercept",names(X)),"(p)",sep="")
	names(Fstat)<-"F-statistic"
	names(Fp)<-"F p-value"
	return(list(r.squared=r.squared,
		coefficients=coeffs,
		tstatistic=tstat,
		tpvalue=tp,
		Fstatistic=Fstat,
		Fpvalue=Fp))
}

# unfold converts the lower diagonal elements of a matrix into a vector
# unfold is called by MMRR

unfold<-function(X){
	x<-vector()
	for(i in 2:nrow(X)) x<-c(x,X[i,1:i-1])
	return(x)
}



mats <- cost.mats
if (!is.null(eucl.mat)) mats$Euclidean <- eucl.mat


res <- MMRR(gen.mat, mats, nperm=nperm)

mmrr.mat <- data.frame(layer=names(res$coefficients), coefficient= res$coefficients, tstatistic=res$tstatistic, tpvalue = res$tpvalue, Fstat= NA, Fpvalue = NA, r2=NA )

row.names(mmrr.mat)=NULL
mmrr.mat <- mmrr.mat[order(mmrr.mat$tpvalue),]
mmrr.mat$Fstat[1] <- res$Fstatistic
mmrr.mat$Fpvalue[1] <- res$Fpvalue
mmrr.mat$r2[1] <- res$r.squared




return(list(mmrr.tab=mmrr.mat))
}
