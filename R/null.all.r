#' Checks for the presence of and determine the frequency of null alleles
#' 
#' The function null.all determines the frequency of null alleles at each locus
#' of a genind object. As an initial step, the function makes a bootstrap
#' estimate (based on the observed allele frequencies) of the probability of
#' seeing the number of homozygotes observed for each allele. If there are a
#' large number of null alleles present at a locus, it would result in multiple
#' alleles at a locus having an excess of homozygotes. The second step of the
#' function estimates the frequency of null alleles and a bootstrap confidence
#' interval for each locus using the methods of Chakraborty et al. (1994) and
#' Brookfield (1996). If the 95\% confidence interval includes zero, it
#' indicates that the frequency of null alleles at a locus does not
#' significantly differ from zero.
#' 
#' 
#' @param population a \code{\link{genind}} object (from the package adegenet)
#' @return The function returns a list with two main components: homozygotes
#' and null.allele.freq.  Homozygotes contains the output of the first part of
#' the analysis: determining the observed number of homozygotes for allele at
#' each locus (homozygotes$observed), generating a distribution of the expected
#' number of homozygotes for each allele at each locus (homozygotes$bootstrap)
#' based upond the observed allele frequencies at a locus, and producing a
#' summary table given the probability of observing the number of homozygotes
#' (homozygotes$probability.obs).  null.allele.freq list contains summary
#' tables of the null allele frequency estimates based upon the forumulas of
#' Chakraborty et al. (1994) (summary1), and Brookfield (1996) (summary2). For
#' each summary table, the observed frequency is the null allele frequency
#' determined from the observed heterozygosity and homozygosity at a locus. The
#' median, 2.5th, and 97.5th percentiles are from bootstrap estimates of null
#' allele frequencies obtained by resampling the individual genotypes from the
#' original genind object.
#' 
#' Brookfield (1996) provides a brief discussion on which estimator should be
#' used. In summary, it was recommended that Chakraborty et al. (1994)'s method
#' (e.g. summary1) be used if there are individuals with no bands at a locus
#' seen, but they are discounted as possible artefacts. If all individuals have
#' one or more bands at a locus then Brookfield (1996)'s method (e.g. summary2)
#' should be used.
#' @author Aaron Adamack, aaron.adamack@@canberra.edu.au
#' @seealso \code{\link{popgenreport}}
#' @references Brookfield JFY. (1996) A simple new method for estimating null
#' allele frequency from heterozygote deficiency. Molecular Ecology 5:453-455
#' 
#' Chakraborty R, Zhong Y, Jin L, Budowle B. (1994) Nondetectability of
#' restriction fragments and independence of DNA fragment sizes within and
#' between loci in RFLP typing of DNA. American Journal of Human Genetics
#' 55:391-401
#' @examples
#' 
#'  \dontrun{
#'  data(bilby)
#'  #here we use only the first 50 individuals to speep up the example
#'  popgenreport(bilby, mk.null.all=TRUE, mk.pdf=FALSE)
#'  
#' #to get a pdf output you need to have a running Latex version installed on your system.
#' #popgenreport(bilby, mk.null.all=TRUE, mk.pdf=TRUE)
#' }
#' @importFrom plyr alply
#' @importFrom reshape2 melt
#' @export
null.all<-function(population)
{
  # Confirm that the function has been provided with a genind object
  if (!is(population,"genind")){
    stop("You did not provide a valid genind object! Script stopped!")

  }
  
  allploidies<-as.numeric(names(table(population@ploidy)))
  
  if(any(allploidies!=2,na.rm=TRUE)){
    stop("One or more populations has a ploidy other than 2! Script stopped!")

  } 
  
  
  # divide the genind object into individual loci
  split<-seploc(population)
  ninds<-dim(population@tab)[1]
  maxalleles<-max(population@loc.n.all)
  
  per.results<-matrix(NA,nrow=length(split),ncol=maxalleles)
  list.obs.ho.cnt<-list()
  list.exp.ho<-list()
  locus_ho_dist<-matrix(NA,nrow=1000,ncol=length(split))
  obs_tot_count<-rep(NA,length(split))
  
  # for each locus...
  for(i in 1:length(split)){
    # for each allele, count the number of that type seen
    allelecnt<-apply(split[[i]]@tab,2,sum,na.rm=TRUE)
    
    # get the number of observed homozygotes for each allele
    obs_ho<-alply(split[[i]]@tab,2,table)
    obs_ho_cnt<-rep(NA,length(obs_ho))
    for (j in 1:length(obs_ho)){
      find1s<-which(names(obs_ho[[j]])=="2") # if this was changed to ploidy, could be generic...
      if(length(find1s)==0) {
        obs_ho_cnt[j]<-0
      } else {
        obs_ho_cnt[j]<-unname(obs_ho[[j]][find1s])
      }
    }
    
    # get the observed number of homozygotes across alleles at a locus
    obs_tot_count[i]<-sum(obs_ho_cnt,na.rm=TRUE)
    
    # calculate the allele frequencies
    allelefreq<-allelecnt/sum(allelecnt)
    
    # get the expected counts of homozygotes for each allele
    # this is assuming HWE...
    numho<-matrix(NA,nrow=1000,ncol=length(allelefreq))
    
    for (k in 1:999){ # this is looping over replicates....
      tempgenotype<-matrix(sample(1:length(allelecnt),sum(allelecnt),replace=TRUE,prob=allelefreq),ncol=2)# generating sets of genotypes
      allelepairtab<-table(tempgenotype[,1],tempgenotype[,2])
      allelepairlong<-reshape2::melt(allelepairtab)
      for(l in 1:length(allelecnt)){
        left<-which(allelepairlong[,1]==l)
        right<-which(allelepairlong[,2]==l)
        matching<-intersect(left,right)
        if (length(matching)==0) {
          numho[k,l]<-0
        } else {
          numho[k,l]<-allelepairlong[matching,3]
        }      
      }
      locus_ho_dist[k,i]<-sum(numho[k,],na.rm=TRUE)
    }
    numho[1000,]<-obs_ho_cnt
    locus_ho_dist[1000,i]<-obs_tot_count[i]
    per.results[i,1:length(obs_ho_cnt)]<-sapply(1:dim(numho)[2],function(x,obs_ho_cnt,numho) sum(numho[,x]>obs_ho_cnt[x])/1000, obs_ho_cnt,numho)
    list.obs.ho.cnt[[i]]<-obs_ho_cnt
    list.exp.ho[[i]]<-numho
  }
  rownames(per.results)<-unname(locNames(population))
  suffix<-seq(1:maxalleles)
  colnames(per.results)<-paste("Allele-",suffix,sep="")
  
  homozygotes<-list(observed=list.obs.ho.cnt,bootstrap=list.exp.ho,probability.obs=per.results,overall=list(observed=obs_tot_count,distribution=locus_ho_dist))
  
  
  # calculate null allele frequencies...
  
  distr1<-matrix(NA,nrow=1000,ncol=length(split))
  distr2<-matrix(NA,nrow=1000,ncol=length(split))
  morethan1<-unname(population@loc.n.all)>1
  for(k in 1:length(morethan1)){
    if (!morethan1[k]) warning("Locus ",unname(locNames(population)[k])," has only 1 allele and null allele frequency will not be estimated for it")
  }
  for (i in 1:999){
    for (j in 1:length(split)){
      if (morethan1[j]){
        #message("i = ",i," j = ",j)
        # randomly draw individuals from the population that have been sampled
        tempalleles<-split[[j]]@tab[sample(1:ninds,ninds,replace=TRUE),]  
        # counting the number of each allele type
        allelecnt<-apply(tempalleles,2,sum,na.rm=TRUE) 
        allelefreq<-allelecnt/sum(allelecnt)
        exphz<-1-sum(allelefreq^2)
        ho_cnt<-reshape2::melt(table(tempalleles))
        if(length(ho_cnt$value[ho_cnt$temp==2])==0) {
          numho<-0
        } else if(length(ho_cnt$value[ho_cnt$temp==2])>0){
          numho<-ho_cnt$value[ho_cnt$temp==2]
        }
        # for numhz, have to subtract number of homozygotes and missing from ninds
        numhz<-ninds-numho-sum(is.na(tempalleles[,1]))
        obshz<-1-numho/(numho+numhz)
        distr1[i,j]<-(exphz-obshz)/(exphz+obshz)
        distr2[i,j]<-(exphz-obshz)/(1+obshz)  
      }
    }
  }
  for (k in 1:length(split)){
    if(morethan1[k]){
      # last observation is for the actual dataset
      tempalleles<-split[[k]]@tab
      allelecnt<-apply(tempalleles,2,sum,na.rm=TRUE)
      allelefreq<-allelecnt/sum(allelecnt)
      exphz<-1-sum(allelefreq^2)
      ho_cnt<-reshape2::melt(table(tempalleles))
      if(length(ho_cnt$value[ho_cnt$tempalleles==2])==0){
        numho<-0
      } else if(length(ho_cnt$value[ho_cnt$tempalleles==2])>0){
        numho<-ho_cnt$value[ho_cnt$tempalleles==2]
      }
      numhz<-ninds-numho-sum(is.na(tempalleles[,1]))
      obshz<-1-numho/(numho+numhz)
      distr1[1000,k]<-(exphz-obshz)/(exphz+obshz)
      distr2[1000,k]<-(exphz-obshz)/(1+obshz)  
    }
  }
  
  null.allele.boot.dist<-list(method1=distr1,method2=distr2)
  
  method1<-matrix(NA,nrow=4,ncol=length(split))
  method2<-matrix(NA,nrow=4,ncol=length(split))
  
  method1[1,]<-distr1[1000,]
  method2[1,]<-distr2[1000,]
  
  method1[2,]<-apply(distr1,2,median,na.rm=TRUE)
  method2[2,]<-apply(distr2,2,median,na.rm=TRUE)
  
  method1[3,]<-apply(distr1,2,quantile,0.025,na.rm=TRUE)
  method1[4,]<-apply(distr1,2,quantile,0.975,na.rm=TRUE)
  method2[3,]<-apply(distr2,2,quantile,0.025,na.rm=TRUE)
  method2[4,]<-apply(distr2,2,quantile,0.975,na.rm=TRUE)
  
  rownames(method1)<-c("Observed frequency","Median frequency","2.5th percentile","97.5th percentile")
  rownames(method2)<-c("Observed frequency","Median frequency","2.5th percentile","97.5th percentile")
  colnames(method1)<-unname(locNames(population))
  colnames(method2)<-unname(locNames(population))
  
  null.allele.freq<-list(summary1=method1,summary2=method2,bootstrap=null.allele.boot.dist)
  results.null.alleles<-list(homozygotes=homozygotes,null.allele.freq=null.allele.freq)
  return(results.null.alleles)
}