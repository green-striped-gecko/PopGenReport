#' Calculates the allelic richness for a genind object
#' 
#' The function calculates the allelic richness for each combination of
#' population and locus for a genind object. To account for differences in
#' sample sizes and genotyping success, rarefication is used in the
#' calculation. The sample size for each combination of population and locus
#' was set equal to the smallest number of alleles seen in a sample across all
#' combinations of population and locus. Allelic richness was calculated using
#' the methods of Mousadik and Petit (1996) which are in turn based upon the
#' work of Hurlbert (1971).
#' 
#' This function is similar to the allelic.richness function in hiefstat. The
#' main differences between the two packages are that allel.rich works on a
#' genind object while allelic richness works on a data frame and allel.rich is
#' capable of determining allelic richness for species with most ploidies while
#' allelic.richness only works for haploid and diploid species.
#' 
#' 
#' @param population a \code{\link{genind}} object (from package adegenet)
#' @param min.alleles the minimum number of alleles that will be sampled for
#' calculating allelic richness. If min.alleles is set to NULL the min.alleles
#' sampled will be determined automatically (see description)
#' @return Returns a list with the following entries:
#' 
#' all.richness is the allelic richness for each combination of population and
#' locus sum.richness is the sum of the allelic richnesses for each population
#' mean.richness is the mean allelic richness across all loci alleles.sampled
#' is the smallest number of individuals sampled across all combinations of
#' population and locus multiplied by the ploidy of the species.  pop.sizes is
#' a matrix with the total number of alleles counted for each combination of
#' population and locus.
#' @author Aaron Adamack, aaron.adamack@@canberra.edu.au
#' @seealso \code{\link{popgenreport}}
#' @references El Mousadik A, Petit RJ. (1996) High level of genetic
#' differentiation for allelic richness among populations of the argan tree
#' [Argania spinosa (L.) Skeels] endemic to Morocco
#' @examples
#' 
#'  #not run:
#'  #data(bilby)
#'  #here we use only the first 50 individuals to speep up the example
#'  #popgenreport(bilby, mk.allel.rich=TRUE, mk.pdf=FALSE)
#'  
#' #to get a pdf output you need to have a running Latex version installed on your system.
#' #popgenreport(bilby, mk.allel.rich=TRUE, mk.pdf=TRUE)
#' 
#' #data(bilby)
#' #allel.rich(bilby)
#' @export


allel.rich<-function(population,min.alleles=NULL)
{
  
  npops<-length(levels(population@pop))
  nloci<-length(locNames(population))
  
  # this splits bilby up into loci
  loci<-seploc(population)
  
  # this further subdivides the loci into populations
  locipop<-lapply(loci,seppop)
  
  popsizes<-matrix(NA,nrow=nloci,ncol=npops)
  for (i in 1:nloci){
    for (j in 1:npops){
      popsizes[i,j]<-sum(!is.na(apply(locipop[[i]][[j]]@tab,1,sum)))
    }
  }
  
  colnames(popsizes)<-unname(popNames(population))
  rownames(popsizes)<-unname(locNames(population))
  
  counter<-0
  for(i in 1:dim(popsizes)[2]){
    numzeroes<-length(which(popsizes[,i]==0))
    if(numzeroes>0) {
      warning("Population ",unname(popNames(population))[i]," has ",numzeroes," locus/loci with no genotypes observed which will cause allel.rich to fail. Please adjust your dataset appropriately")
      counter<-counter+1
    }
  }
  if(counter>0){
    message("Please see the pop.sizes matrix in the output list to identify the combinations of population and locus causing problems")
  }
  
  
  richness<-matrix(NA,nrow=nloci,ncol=npops)
  if(is.null(min.alleles)) {
    g<-min(popsizes)*population@ploidy[1]
  } else if(!is.null(min.alleles)){
    g<-min.alleles
  }
  for (i in 1:nloci){ 
    for (j in 1:npops){
      allelecnt<-apply(locipop[[i]][[j]]@tab,2,sum, na.rm=TRUE)*population@ploidy[1]
      richness[i,j]<-sum(1-choose((sum(allelecnt)-allelecnt),g)/choose(sum(allelecnt),g))
    }
  }
  
  variance<-matrix(NA,nrow=nloci,ncol=npops)
  for (i in 1:nloci){
    for (j in 1:npops){
      allelecnt<-apply(locipop[[i]][[j]]@tab,2,sum, na.rm=TRUE)*population@ploidy[1]
    }
  }                               
  
  
  
  
  colnames(richness)<-popNames(population)
  rownames(richness)<-locNames(population)
  srichness<-apply(richness,2,sum)
  mrichness<-apply(richness,2,mean)
  
  names(srichness)<-popNames(population)
  names(mrichness)<-popNames(population)
  # bump up the allele count to recognize the ploidy of individuals
  popsizes2<-popsizes*population@ploidy[1]
  allelic.richness<-list(all.richness=richness,sum.richness=srichness,mean.richness=mrichness,alleles.sampled=g,pop.sizes=popsizes2)
  
  return(allelic.richness)
}