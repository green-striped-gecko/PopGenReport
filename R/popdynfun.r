################################################################################
##### function to convert cost distances into dispersal probabilities
################################################################################
#'Function to calculate dispersal distances based on cost distances
#'
#'@param x Cost (Euclidean) distance matrix
#'@param d0 dispersal distance
#'@param of p of all individuals in a population
#'@return a dispersal probability matrix
#'@description converts cost distances to probabilities: to reach a certain patch[ d0 average distance of p of all individuals, for example d0=100, p =0.5 -> 50\% procent of all migrating individuals go up to 100 m.

 p2p <- function(x, d0, p)
 {
 return (exp(((x/d0)*log(p))))
 }




################################################################################
##### functions  of popdynamics
################################################################################
#mutation
#'Function to execute mutation on a pop data.frame
#'
#'@param x a pop object
#'@param n.allels number of alleles 
#'@param mutrate the mutation rate
#'@param mtype the mutation process (currently on n.all is implemented. (=Choose any allel from 1:n.allels))
#'@param n.cov (number of covariates, default  )
#'@return a dispersal probability matrix
#'@description mutation subprocess on single populations

mutation <- function(x, n.allels, mutrate, mtype="n.all", n.cov)
{

n.indp <- nrow(x)
n.loci <-  ( ncol(x)-n.cov)  / 2
p.mut = mutrate
noffset = n.cov


if (!is.null(n.indp) )    #some populations might have gone extinct...
{
mut.ind <- which(matrix(as.logical( rbinom(n.indp*n.loci,1,p.mut)), ncol=n.loci, nrow=n.indp)==TRUE,arr.ind=T )

if (nrow(mut.ind)>0) #are there any mutations in this population?
{
for (i in 1:nrow(mut.ind))
  {
  
  if (mtype=="n.all")
    {
    mf <- rbinom(1,1,0.5)
    x[mut.ind[i,1], noffset+mut.ind[i,2]*2-mf] <- sample(1:n.allels,1)
    }
  }
}
}
x
}


#reprodution

reproduction <- function(x, type="K.limit",K=n.ind, n.off, n.cov)
{
noffset <- n.cov
n.loci <-  ( ncol(x)-n.cov)  / 2
n.ind <- nrow(x)

index.male <- x$sex=="male"
index.female <- x$sex=="female"
n.female <- sum(index.female)
n.male <-  sum(index.male)
#offsprings  <- data.frame(cbind(pop=NA, sex=NA, age=NA,matrix(NA,ncol=n.loci*2,nrow=1)) )


if (n.male>0 & n.female>0)
{
loci.cols <- (noffset+1):(n.loci*2+noffset)

genes.females <- (x[index.female,loci.cols])

#for (i in 1:(n.off-1)) genes.females <- rbind(genes.females,)

genes.females <- genes.females[rep(1:nrow(genes.females),n.off),]

pick.femall <-sample(c(T,F),n.female*n.off*n.loci,replace=T)
all.fem <-  matrix(as.vector(rbind(pick.femall,!pick.femall)),nrow=n.female*n.off, ncol=n.loci*2, byrow=T)


genes.off.fem <- matrix(t(genes.females)[t(all.fem)], ncol=n.loci, nrow=n.female*n.off, byrow=T)



breeding.males <- sample(which(index.male),n.female*n.off, replace=T)
genes.males <- x[breeding.males,loci.cols]
pick.malall <-sample(c(T,F),n.female*n.off*n.loci,replace=T)

all.mal <-  matrix(as.vector(rbind(pick.malall,!pick.malall)),nrow=n.female*n.off, ncol=n.loci*2, byrow=T)
genes.off.mal <- matrix(t(genes.males)[t(all.mal)], ncol=n.loci, nrow=n.female*n.off, byrow=T)


popn <- x$pop[1]
off.sex <- factor( ifelse(runif(n.female*n.off)<0.5,"male","female"))
off.age <- 1

all.seq<- NA

for (i in seq(1,(2*n.loci-1),2))
{ 

if (i==1) all.seq[i] <- 1 else all.seq[i] <- all.seq[i-2]+1
all.seq[i+1] <- all.seq[i]+n.loci
}

off.loci <- cbind(genes.off.mal, genes.off.fem)[,all.seq]

colnames(off.loci) <-paste("locus", rep(seq(1,n.loci,1),each=2),rep(c("A","B"),n.loci),sep="")

offsprings <- data.frame(pop=popn, sex=off.sex , age=off.age, off.loci)


#combine adults and offsprings
x <- rbind(x, offsprings)
#x<- rbindlist(list(x,offsprings))

#cut off at K

npop <- nrow(x)
x<- x[sample(1:npop,K, replace=F),]

} else x<- NULL    #return NULL if only males or only females....

return(x)
}





#survival   (cut of at K)
#surv <- function(x, s=t.surv) { 
#
#x<- x[runif(dim(x)[1])< s,]
#}
#lapply(pops, surv, s=0.1)

################################################################################
#emigration (for all populations at once xp=pops)
emigration <- function(xp=pops,n.ind = n.ind, perc.mig=p.mig, emi.m=emis, emi.table=NULL)
{


n.pops = length(xp)

migs <- matrix(0,nrow=n.pops, ncol=n.pops)
migrants<- NA


if (is.null(emi.table))
{
for (i in 1:n.pops)
{
pop.size <-nrow(xp[[i]])

if (!is.null(pop.size))
{
#how many migrants
migrants[i] <- sum(ifelse(runif(pop.size)<perc.mig,1,0))
#where do they go to (depends on emis)
fromto <- table(sample(1:n.pops,migrants[i],replace=T, prob=emi.m[,i]))
to <- as.numeric(names(fromto))
for (ii in 1:length(to)) migs[i,to[ii ]] <- migs[i,to[ii ]] + fromto[ii]
}
}
} else migs <- emi.table #end emi.type


#deterministic emigration (needs to be passed by emi.table


#let them move (at the moment an individual can migrate more than once, need to be fixed)

for (from in 1:n.pops)
{

for (to in 1:n.pops)
  {
  #update population size during migration events!!!
  psize <- nrow(xp[[from]])
  m <- migs[from, to]
  if (m>0)
    {
    ind.from <- sample(1:psize,m,replace=F)
    xp[[to]] <- rbind(xp[[to]], xp[[from]][ind.from,])
    xp[[from]] <- xp[[from]][-ind.from,]
    xp[[to]][,"pop"] <- to
    }
  }
}
#output number of migrants to a population
#diag(migs) <- 0
#cat(paste("Migrants: ",colSums(migs),"\n"))

results <- list(xp, migs)
return(results)

}


#combine all pops into a single genind object.

################################################################################
pops2genind <- function(x, prop, locs=NULL, n.cov)
{         


n.pops <- length(x)
n.loci <- (ncol(x[[1]])-n.cov )/2 #only diploid currently

noffset <- n.cov

##### convert pops to data.frame
#combine <-do.call(rbind.data.frame,x)
combine <- do.call(rbind.data.frame,x)
#combine <- data.frame(rbindlist(x))

#convert to genind
allele <-combine[,-(1:noffset)]

pop.size <- table(combine[,"pop"] )

res <- data.frame(matrix(NA, nrow=nrow(allele), ncol=n.loci))#rep(NA,sum(pop.size)))
        for (i in seq(1, 2*n.loci, 2))
            res[, ceiling(i/2)] <- paste(allele[, i], allele[, i + 1], sep = "/")

pops.genind <- df2genind(res,pop=combine$pop, sep="/")

#add locations
if (!is.null(locs))
  {
  x <- rep(locs[,1], pop.size)
  y <- rep(locs[,2], pop.size)
  pops.genind@other$xy <- cbind(x=x,y=y)
  pops.genind@pop.names <- row.names(locs)
  }
return(pops.genind)
}



################################################################################
### sim functions
################################################################################

#     run simulation a number of steps (generations)
run.popgensim <- function(simpops, steps, cost.mat, n.offspring , n.ind,  mig.rate, disp.max, disp.rate, mut.rate, n.cov, rec = "none", emi.table=NULL)
{


pops <- simpops


distances <- cost.mat
disdis <- p2p(distances,d0=disp.max, p=disp.rate)

stdemi.table= disdis/rep(colSums(disdis), each=ncol(disdis))
emi.tableint=NULL


#cl<- makePSOCKcluster(6)

for (year in 1:steps)
{

if (!is.null(emi.table))
  {
  emi.table <- stdemi.table * mig.rate * n.ind +  emi.table
  emi.tableint <- floor(emi.table)
  emi.table <- emi.table - emi.tableint
  }


#emigration
emi<- emigration(pops, n.ind, mig.rate, emi.m=disdis, emi.table=emi.tableint)
pops <- emi[[1]]
#emi[[2]] table of migrants....

#reproduction
 pops <- lapply(pops, reproduction,K=n.ind, n.off=n.offspring, n.cov=n.cov)
 # pops <- parLapply(cl, pops, reproduction, K=n.ind, n.off=n.offspring, n.cov=n.cov)
#mutation
pops <- lapply(pops, mutation, n.allels, mutrate=mut.rate, n.cov=n.cov)
}  #end of every year....   time loop

#stopCluster(cl)

out <- pops #list(pops=pops)
if (rec=="emi")
{
out <- list()
out[[1]] <- pops
out[[2]] <- emi[[2]]
}
return (out)
} #end of simpop function

################################################################################
#initialise populations
init.popgensim <- function(n.pops, n.ind, sex.ratio, n.loci, n.allels, locs, n.cov)
{

noffset <- n.cov #number of covariates before the loci

#initialise populations
empty.loci <- matrix(NA,ncol=n.loci*2,nrow=n.ind)
colnames(empty.loci) <-paste("locus", rep(seq(1,n.loci,1),each=2),rep(c("A","B"),n.loci),sep="")
pops <- list( )
for (i in 1:n.pops)
{
pops[[i]]<- data.frame(cbind(pop=NA, sex=NA, age=NA,empty.loci))
#population number
pops[[i]][,"pop"] <- i
#sex
pops[[i]][,"sex"] <- "male"
pops[[i]][1:(round(n.ind*sex.ratio,0)),"sex"] <- "female"

pops[[i]][,"sex"] <- factor(pops[[i]][,"sex"])

#uniform allel distribution
pops[[i]][,(noffset+1):(n.loci*2+noffset) ]  <- sample(1:n.allels,n.loci*2*n.ind, replace=T)
}

names(pops) <- row.names(locs)   #or from genind object....
#gpops <- pops2genind(x=pops, prop=1 , locs=locs, n.cov=noffset)

return (pops)
}


################################################################################
### #function to create cost.mat from landscape and populations coordinates
################################################################################


costdistances <- function(landscape, locs, method, NN)
{
fric.mat <- transition(landscape,function(x) 1/x[2],NN)

# fric.mat <- transition(fr.raster,function(x) 1/(abs(x[1]-x[2])),8)
#set distances to meters  if no projected already
fric.mat@crs@projargs<- "+proj=merc +units=m"
fric.mat.cor <- geoCorrection(fric.mat)

if (method=="leastcost") cd.mat <-costDistance(fric.mat.cor, locs, locs)
if (method=="rSPDistance") cd.mat <- rSPDistance(fric.mat.cor, locs, locs, theta=1)
if (method=="commute") cd.mat <-as.matrix(commuteDistance(fric.mat.cor,locs))

colnames(cd.mat) <- row.names(locs)
rownames(cd.mat) <- row.names(locs)

return (cd.mat)
}

################################################################################
### barallel pairwise fst
################################################################################


  pairwise.fstb <- function (x, pop = NULL, res.type = c("dist", "matrix"), truenames = TRUE)
{
#    require(doSNOW)
#require(foreach)
#cl<-makeCluster(6) #change the 2 to your number of CPU cores
#registerDoSNOW(cl)
#
    res.type <- match.arg(res.type)
    f1 <- function(pop1, pop2) {
        n1 <- nrow(pop1@tab)
        n2 <- nrow(pop2@tab)
        temp <- repool(pop1, pop2)
        b <- weighted.mean(Hs(temp), c(n1, n2))
        pop(temp) <- NULL
        a <- Hs(temp)
        return((a - b)/a)
    }
    lx <- seppop(x, treatOther = FALSE)
    temp <- pop(x)
    levPop <- levels(temp)
    allPairs <- combn(1:length(levPop), 2)
    if (!is.matrix(allPairs)) {
        allPairs <- matrix(allPairs, nrow = 2)
    }
    vecRes <- numeric()

    #for (i in 1:ncol(allPairs)) {
    vecRes <- foreach(i=1:ncol(allPairs), .combine=cbind, .packages='adegenet') %dopar% {
        #vecRes[i] <- f1(lx[[allPairs[1, i]]], lx[[allPairs[2,
         #   i]]])
         f1(lx[[allPairs[1, i]]], lx[[allPairs[2, i]]])
    }
    squelres <- dist(1:length(levPop))
    res <- as.numeric(vecRes)
    attributes(res) <- attributes(squelres)
    if (res.type == "matrix") {
        res <- as.matrix(res)
        if (truenames) {
            lab <- x@pop.names
        }
        else {
            lab <- names(x@pop.names)
        }
        colnames(res) <- rownames(res) <- lab
    }
#    stopCluster(cl)
    return(res)
}

