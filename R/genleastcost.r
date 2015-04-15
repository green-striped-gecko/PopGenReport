#leatcost function
genleastcost <- function(cats, fr.raster, gen.dist, NN=NULL, pathtype="leastcost", plotpath=TRUE, theta=1)
{
  if (is.null(NN) & pathtype=="leastcost") 
  {
    cat("NN is not specified!\nPlease specify the number of nearest neighbour to use for the least-cost path calculations (NN=4 or NN=8). If linear features are tested you may want to consider NN=4 otherwise NN=8 is the most commonly used and prefered option. In any case check the actual least-cost paths for artefacts by inspecting the plot on least-cost paths.\n")
    return()
  }
  
dist.type<-NA
if (gen.dist=="D" || gen.dist=="Gst.Hedrick" || gen.dist=="Gst.Nei") dist.type<- "pop" 

if (gen.dist=="Kosman" || gen.dist=="Smouse" || gen.dist=="propShared") dist.type<- "ind" 

if (is.na(dist.type)) 
  {cat("No valid genetic distance type was provided. Please check ?landgenreport for valid options\n")
   return(-1)
}

if (dist.type=="pop")
{
#calculate the centers if population meassurment is wanted
c.x <- tapply(cats@other$xy[,1],cats@pop, mean)
c.y <- tapply(cats@other$xy[,2],cats@pop, mean)
cp<-cbind(c.x, c.y)
} else 
{
cp <- cbind(cats@other$xy[,1], cats@other$xy[,2])
}


eucl.mat <- round(as.matrix(dist(cp)),3)

if (dist.type=="pop") 
{

  dimnames(eucl.mat) <- list(cats@pop.names, cats@pop.names) 
  npop <- length(levels(cats@pop))
} else
  
{
  dimnames(eucl.mat) <- list(cats@ind.names, cats@ind.names) 
  npop <- length(cats@ind.names)
}





#check if fr.raster is a stack or not...
mats <- list()
mats.names<- NA
mats.pathlength<- list()
mats.paths<- list()

pathlength.mat <- NULL
paths <- NULL



n.mats <- dim(fr.raster)[3] #number of rasters in the stack



for (ci in 1:n.mats)
{

plot(fr.raster[[ci]], main=paste(names(fr.raster)[ci],":",pathtype,", NN=",NN,sep=""))
 #image(fr.raster, col=fr.raster@legend@colortable, asp=1)

 mapcolor <-   col2rgb(cats@other$mapdotcolor)/255
 if (is.null(mapcolor)) colX= rgb(0,0,1,0.8) else 
 colX = rgb(mapcolor[1,],mapcolor[2,], mapcolor[3,],cats@other$mapdotalpha)
points(cats@other$xy,cex=1, pch=16, col="blue")
# points(cats@other$xy,cex=cats@other$mapdotsize, pch= cats@other$mapdottype, col=colX)
if (dist.type=="pop")  points(cp,cex=cats@other$mapdotsize*1.5 , pch= 16, col="black")


#create friction matrix
fric.mat <- transition(fr.raster[[ci]],function(x) 1/x[2],NN)

# fric.mat <- transition(fr.raster,function(x) 1/(abs(x[1]-x[2])),8)
#set distances to meters  if no projected already
fric.mat@crs@projargs<- "+proj=merc +units=m"
fric.mat.cor <- geoCorrection(fric.mat)


if (pathtype=="leastcost") cd.mat <-costDistance(fric.mat.cor, cp, cp)

if (pathtype=="rSPDistance") cd.mat <- rSPDistance(fric.mat.cor, cp, cp, theta=1)

if (pathtype=="commute") cd.mat <-as.matrix(commuteDistance(fric.mat.cor, cp))


  dimnames(cd.mat) <- dimnames(eucl.mat) 




if (pathtype=="leastcost" & plotpath==TRUE)   #only show paths if leastcost otherwise not possible
{
comb <- t(combn(1:npop,2))


#pathlength matrix
pathlength.mat <- cd.mat
pathlength.mat[,] <- 0
paths<- list()

cols <- rainbow(dim(comb)[1], alpha=0.5)
for (i in 1:dim(comb)[1])
{

if (dist(rbind(cp[comb[i,1],], cp[comb[i,2],]))==0)
{
 ll <- Line(rbind( cp[comb[i,1],], cp[comb[i,2],]))
 S1 <- Lines(list(ll),ID="Null")
 sPath <- SpatialLines(list(S1))
} else 
{
sPath <- shortestPath(fric.mat.cor, cp[comb[i,1],], cp[comb[i,2],], output="SpatialLines")
}

lines(sPath, lwd=1.5, col=cols[i])
paths[[i]] <- sPath
ll <-  round(SpatialLinesLengths(sPath),3)
pathlength.mat[comb[i,1],comb[i,2]] <- ll
pathlength.mat[comb[i,2],comb[i,1]] <- ll
}

}

mats[[ci]] <- cd.mat
mats.names[[ci]] <- names(fr.raster)[ci]
mats.pathlength[[ci]] <- pathlength.mat
mats.paths[[ci]] <- paths



} #end of ci loop

#mats[[n.mats+1]] <- eucl.mat
#mats.names[n.mats+1]<- "Euclidean"
#
names(mats)  <- names(fr.raster)
#put other calculations here....
# Calculate genetic distances across subpopulations

if (gen.dist=="Gst.Nei")
{
gendist.mat<-as.matrix(pairwise_Gst_Nei(cats))
}
if (gen.dist=="Gst.Hedrick")
{
gendist.mat<-as.matrix(pairwise_Gst_Hedrick(cats))
}
if (gen.dist=="D")
{
gendist.mat<-round(as.matrix(pairwise_D(cats)),4)
}

if (gen.dist=="Smouse")
{
gendist.mat <- as.matrix(gd.smouse(cats,verbose=FALSE))
}
if (gen.dist=="Kosman")
{
gendist.mat <-as.matrix(as.dist(gd.kosman(cats)$geneticdist))
}
if (gen.dist=="propShared")
{
  gendist.mat <-as.matrix(as.dist(propShared(cats)))
}

  dimnames(gendist.mat)<-dimnames(eucl.mat)


return(list( gen.mat=gendist.mat, eucl.mat=eucl.mat,cost.matnames=mats.names, cost.mats=mats,pathlength.mats= mats.pathlength,  paths=mats.paths))
}


