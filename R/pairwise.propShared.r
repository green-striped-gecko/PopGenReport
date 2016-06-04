################################################################################
### #function to create cost.mat from landscape and populations coordinates
################################################################################
################################################################################
#'Calculates proportion of shared alleles per pairs of populations
#'
#'@param gi a genind object with at least two populations
#'@return a matrix of proportion of shared alleles between populations
#'@description Calculates proportion of shared alleles per pairs of populations based on the minima of allele frequency for each allel (then summed and averaged over loci). Returns a similarity matrix (upper diagonal of the pairwise matrix.

pairwise.propShared <- function(gi)
{
  n.pops <- length(unique(pop(gi)))
  allPairs <- combn(1:n.pops, 2)
  gen.mat<- matrix(0, nrow=n.pops, ncol=n.pops)
  pops <- seppop(gi)
  pspop <- function(x) 
  {
    pp  <- seppop(x)
    p1 <- pp[[1]]
    p2 <- pp[[2]]
    
    na <- ncol(p1@tab)
    maf <- NA
    m1 <- colMeans(p1@tab[,], na.rm=T)/2
    m2 <- colMeans(p2@tab[,], na.rm=T)/2
    
    m12 <- apply(rbind(m1,m2), 2, min, na.rm=T)
    
    lfl <- NA
    facs <- levels(p1@loc.fac)
    for (i in 1:length(locNames(p1))) 	lfl[i] <- sum(m12[p1@loc.fac==facs[i]])
    mean(lfl, na.rm=T)	
  }
  
  for (i in 1:ncol(allPairs))
  {
    np1 <- allPairs[1,i]
    np2 <- allPairs[2,i]
    
    p12 <- repool(pops[[np1]], pops[[np2]])
    ps <- pspop(p12)
    gen.mat[np1,np2] <- ps
    gen.mat[np2,np1] <- ps
    
  }
  la <- levels(pop(gi))
  colnames(gen.mat) <- rownames(gen.mat) <- la
  return(as.dist(gen.mat))
}
#######################
