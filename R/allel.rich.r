allel.rich<-function(population,min.alleles=NULL)
{
  
  npops<-length(levels(population@pop))
  nloci<-length(population@loc.names)
  
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
  
  colnames(popsizes)<-unname(population@pop.names)
  rownames(popsizes)<-unname(population@loc.names)
  
  counter<-0
  for(i in 1:dim(popsizes)[2]){
    numzeroes<-length(which(popsizes[,i]==0))
    if(numzeroes>0) {
      warning("Population ",unname(population@pop.names)[i]," has ",numzeroes," locus/loci with no genotypes observed which will cause allel.rich to fail. Please adjust your dataset appropriately")
      counter<-counter+1
    }
  }
  if(counter>0){
    message("Please see the pop.sizes matrix in the output list to identify the combinations of population and locus causing problems")
  }
  
  
  richness<-matrix(NA,nrow=nloci,ncol=npops)
  if(is.null(min.alleles)) {
    g<-min(popsizes)*population@ploidy
  } else if(!is.null(min.alleles)){
    g<-min.alleles
  }
  for (i in 1:nloci){ 
    for (j in 1:npops){
      allelecnt<-apply(locipop[[i]][[j]]@tab,2,sum, na.rm=TRUE)*population@ploidy
      richness[i,j]<-sum(1-choose((sum(allelecnt)-allelecnt),g)/choose(sum(allelecnt),g))
    }
  }
  
  variance<-matrix(NA,nrow=nloci,ncol=npops)
  for (i in 1:nloci){
    for (j in 1:npops){
      allelecnt<-apply(locipop[[i]][[j]]@tab,2,sum, na.rm=TRUE)*population@ploidy
    }
  }                               
  
  
  
  
  colnames(richness)<-population@pop.names
  rownames(richness)<-population@loc.names
  srichness<-apply(richness,2,sum)
  mrichness<-apply(richness,2,mean)
  
  names(srichness)<-population@pop.names
  names(mrichness)<-population@pop.names
  # bump up the allele count to recognize the ploidy of individuals
  popsizes2<-popsizes*population@ploidy
  allelic.richness<-list(all.richness=richness,sum.richness=srichness,mean.richness=mrichness,alleles.sampled=g,pop.sizes=popsizes2)
  
  return(allelic.richness)
}