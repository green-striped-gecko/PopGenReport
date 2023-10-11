#' Create a landscape genetic report
#' 
#' This function is the landscape genetic version of the
#' \code{\link{popgenreport}} function. It needs to be provided with a genind
#' object with spatial coordinates, a friction map (raster) and a specification
#' which type of genetic distance should be used.  Once all three type of input
#' are provided with the necessary input, a landscape genetic analysis using
#' least cost path analysis is computed (see Cushman et al. 2010, Landguth et
#' al. 2010). Depending on the genetic distance meassurement this is done on a
#' subpopulation basis (D, Gst.Hedrick, Gst.Nei=Fst) or on an individual basis
#' (Kosman, Smouse).
#' 
#' Check the help pages of \code{\link{popgenreport}} how to include
#' coordinates to a genind object. The coordinates need to be projected.
#' Latlongs are not valid, because Euclidean distances are calcuated based on
#' these coordinates. For an example how to convert latlongs into a projected
#' format have a look at the vignette that comes with this package. The
#' friction needs to be a raster and needs to be in the same projection as the
#' genind object. Also the type of genetic distance to be used needs to be
#' specified.
#' 
#' @param cats a \code{genind} object with spatial coordinates in the other
#' slot
#' @param fric.raster friction (resistance) raster, that specifies the
#' landscape where the analysis should be computed on. If fric.raser is a stack
#' a cost distances are calculated for each layer in the stack.
#' @param gen.distance type of genetic distance that should be used. Depending
#' on the genetic distance meassurement this is done on a subpopulation basis
#' (D, Gst.Hedrick, Gst.Nei=Fst) or on an individual basis (Kosman, Smouse,
#' propShared). propShared is the proportion of shared alleles between
#' individuals.
#' @param NN Number of neighbours used when calculating the cost distance
#' (possible values 4,8 or 16). As the default is NULL a value has to be
#' provided if pathtype is 'leastcost'. NN=8 is most commonly used as it avoids
#' a directional bias, but be aware that linear structures may cause artefacts
#' in the least-cost paths in the NN=8 case, therefore we strongly recommend to
#' inspect the actual least-cost paths in the provided output.
#' @param pathtype Type of cost distance to be calculated (based on function in
#' the \code{\link{gdistance}} package. Available distances are 'leastcost',
#' 'commute' or 'rSPDistance'. See functions in the gdistance package for
#' futher explanations.
#' @param plotpath switch if least cost paths should be plotted (works only if
#' pathtype='leastcost'. Be aware this slows down the computation, but it is
#' recommended to check least cost paths visually.
#' @param theta value needed for rSPDistance function. see
#' \code{\link{rSPDistance}} in package \code{gdistance}.
#' @param mk.resistance switch to do the landscape genetic analysis based on
#' resistance matrices, should be set to TRUE
#' @param mapdotcolor see \code{\link{popgenreport}}
#' @param mapdotsize see \code{\link{popgenreport}}
#' @param mapdotalpha see\code{\link{popgenreport}}
#' @param mapdottype see \code{\link{popgenreport}}
#' @param mapzoom see \code{\link{popgenreport}}
#' @param mk.custom switch to add a customised part to the landgenreport
#' @param fname see \code{\link{popgenreport}}
#' @param foldername see \code{\link{popgenreport}}
#' @param path.pgr see \code{\link{popgenreport}}
#' @param mk.Rcode see \code{\link{popgenreport}}
#' @param mk.complete see \code{\link{popgenreport}}
#' @param mk.pdf see \code{\link{popgenreport}}
#' @return Four distance matrices are returned. Pairwise Euclidean distances
#' between subpopulations/individuals, cost distances, path lengths and genetic
#' distances. Also following the approach of Wassermann et al. 2010 a series of
#' partial mantel tests are performed. A multiple regression analysis based on
#' Wang 2013 and Legendre 1994 is returned.The actual least-cost paths can be
#' found under paths
#' @author Bernd Gruber (bernd.gruber@@canberra.edu.au)
#' @seealso \code{\link{popgenreport}}, \code{\link{wassermann}},
#' \code{\link{genleastcost}}, \code{\link{lgrMMRR}}
#' @references Cushman, S., Wasserman, T., Landguth, E. and Shirk, A. (2013).
#' Re-Evaluating Causal Modeling with Mantel Tests in Landscape Genetics.
#' Diversity, 5(1), 51-72.
#' 
#' Landguth, E. L., Cushman, S. A., Schwartz, M. K., McKelvey, K. S., Murphy,
#' M. and Luikart, G. (2010). Quantifying the lag time to detect barriers in
#' landscape genetics. Molecular ecology, 4179-4191.
#' 
#' Wang,I 2013. Examining the full effects of landscape heterogeneity on
#' spatial genetic variation: a multiple matrix regression approach for
#' quantifying geographic and ecological isolation. Evolution: 67-12:
#' 3403-3411.
#' 
#' Wasserman, T. N., Cushman, S. A., Schwartz, M. K. and Wallin, D. O. (2010).
#' Spatial scaling and multi-model inference in landscape genetics: Martes
#' americana in northern Idaho. Landscape Ecology, 25(10), 1601-1612.
#' @examples
#' 
#' \donttest{
#' library(raster)
#' fric.raster <- readRDS(system.file("extdata","fric.raster.rdata", package="PopGenReport"))
#' lc<-landgenreport(cats=landgen, fric.raster=fric.raster,
#' gen.distance="D", NN=4, mk.resistance=TRUE, mk.pdf=FALSE)
#' names(lc$leastcost)
#' }
#' @export
#' @importFrom vegan mantel.partial permustats 
#' @importFrom gdistance transition geoCorrection shortestPath costDistance commuteDistance rSPDistance 
#' @importFrom sp SpatialLinesLengths Line Lines SpatialLines SpatialPolygons Polygons Polygon coordinates CRS

landgenreport <- function(cats,
                          fric.raster,  #friction matrix   
                          gen.distance = "Gst.Nei",  #Gst Hedrick, Gst_Nei, smouse, kosman, D, %allels shared 
                          NN=NULL, 
                          pathtype="leastcost",
                          plotpath=TRUE, 
                          theta=1, #for randomSP distance
                          mk.resistance=TRUE,    #mantel tests 
                        
                            mapdotcolor ="blue",
                            mapdotsize=1,
                            mapdotalpha=0.4,
                            mapdottype=19 ,
                            mapzoom=NULL,                          
# "roadmap","mobile","satellite","terrain","hybrid","mapmaker-roadmap","mapmaker-hybrid"

                          
                          mk.custom = FALSE,
                          

                          fname="LandGenReport",
                          foldername="results",
                          path.pgr=NULL,
                          mk.Rcode=FALSE,       # make the code that was ran available as an R file
                          mk.complete=FALSE,    # create a full report)  
                          mk.pdf=TRUE)
{
  if (!is(cats,"genind")) {stop("You did not provide a valid catsnd object! Script stopped!\n")}
  # Check for combinations of populations and loci with low numbers of individuals and alleles  
  npops<-length(levels(cats@pop))
  nloci<-length(locNames(cats))
  
  # this splits bilby up into loci
  loci<-seploc(cats)
  
  # this further subdivides the loci into populations
  locipop<-lapply(loci,seppop)
  
  popsizes<-matrix(NA,nrow=nloci,ncol=npops)
  for (i in 1:nloci){
    for (j in 1:npops){
      popsizes[i,j]<-sum(!is.na(apply(locipop[[i]][[j]]@tab,1,sum)))
    }
  }
  
  for(i in 1:dim(popsizes)[2]){
    numlow<-length(which(popsizes[,i]<3))
    if(numlow>0) message("Population ",unname(popNames(cats))[i]," has ",numlow," locus/loci with less than 3 genotypes. This may cause errors in some analyses. We advice to combine or drop populations with low numbers of genotypes. ")
  }
  
  
#   #cut down length of loci names to  6  and make sure they are unique
#   #locNames(cats) <- substr(locNames(cats),1,6)   
#   # if (length(unique(locNames(cats)))!= length(locNames(cats))) 
#   #   {
#   #   locNames(cats) <- paste(1:length(locNames(cats)),"-",substr(locNames(cats),1,4), sep="")
#   # 
#   # message("Loci names were not unique and therefore adjusted.\n")
#   #   }
#   # levels(cats@loc.fac) <- locNames(cats)  #make sure levels and factors are the same
# #check if indnames are unique!!!!
# #adjust if necessary and issue a notification
# if (length(unique(indNames(cats)))!=length(indNames(cats))) 
#   {indNames(cats) <- paste(1:length(indNames(cats)),"-",substr(indNames(cats),1,8),sep="")
#   message("Individual names were not unique and therefore adjusted.\n")
#   }
# 
# 
# #check if pop.names are unique!!!!
# #adjust if necessary and issue a notification  
# if (length(unique(popNames(cats)))!=length(popNames(cats))) 
#   {
#   popNames(cats) <- paste(1:length(popNames(cats)),"-",substr(popNames(cats),1,6),sep="")
#   message("Subpopulation names were not unique and therefore adjusted.\n")
#   }


if (is.null(NN)) 
{
  stop("NN is not specified!\nPlease specify the number of nearest neighbour to use for the least-cost path calculations (NN=4 or NN=8). If linear features are tested you may want to consider NN=4 otherwise NN=8 is the most commonly used and prefered option. In any case check the actual least-cost paths for artefacts by inspecting the plot on least-cost paths.\n")
}


 
#set directory where to save a file, defaults to tempdir (follow R policy)
  if (is.null(path.pgr)) 
  {
  path.pgr <- tempdir()
  }

  
  #create a foldername folder if not existing...
  dirfiles <- list.dirs(path=path.pgr, recursive=FALSE)
  if (!(tolower (file.path(path.pgr,foldername))) %in% tolower(dirfiles)) {
    dir.create(file.path(path.pgr,foldername))
    message("There is no ",foldername, " folder. I am trying to create it; \notherwise please create the folder manually. \n")
  }
  owd <-getwd()
  on.exit(setwd(owd))
  setwd(file.path(path.pgr, foldername))
  
 
  # coordinates must be in xy !!!!!
  coords=FALSE
  if (is.null(cats@other$xy)) coords=FALSE 
  if ( nrow(cats@other$xy) == length(indNames(cats)) ) coords=TRUE
 
  # give cats a filename that can be seen in the snw chunks
  cats@other$filename<- fname
  cats@other$foldername<-foldername
  cats@other$path <- path.pgr
  #determine the type of map
  if (coords) 
  {

  cats@other$mapdotcolor =mapdotcolor
  cats@other$mapdotsize=mapdotsize
  cats@other$mapdotalpha=mapdotalpha
  cats@other$mapdottype=mapdottype
  cats@other$mapzoom=mapzoom
  
  }  
###################################
##### create a new environment to run knitr in it
pgr <- new.env()
assign("cats",cats,envir=pgr)
assign("gen.distance",gen.distance,envir=pgr)
assign("fric.raster",fric.raster,envir=pgr)
assign("NN",NN,envir=pgr)
assign("pathtype",pathtype,envir=pgr)
assign("theta",theta,envir=pgr)
assign("plotpath",plotpath,envir=pgr)

###################################




  # save the data in a tempfile
 # save(cats, file=paste(foldername,"\\","tempcats.rdata",sep=""))
  
  #check path to the snw files
path <- NULL
  for(i in seq_along(.libPaths()))
{
  if (file.exists(paste(.libPaths()[i],"/PopGenReport/swchunks/header.snw",sep="")))  
  {
  path <-   paste(.libPaths()[i],"/PopGenReport/swchunks/", sep="" )
  break
  }
  
}
if (is.null(path)) {stop("Could not find snw files in the PopGenReport library folder. Please check if the package is installed correctly [e.g.  installed.packages() ]. \n")}
  #for testing:
  #path <- "d:\\bernd\\R\\popgenreport\\inst\\swchunks\\"
  #path<- "C:\\Aaron files\\popgenreport098\\PopGenReport_0.98\\PopGenReport\\swchunks\\"
  header.file <- readLines(paste(path,"header.snw",sep=""))
  required<- readLines(paste(path,"required.snw",sep=""))
  compl<-c(header.file,required) 
  
  message("Compiling report...\n")
  if (coords==FALSE) message(" - No valid coordinates were provided. \n   Be aware you need to provide a coordinate (or NA) for each individual\n   and the coordinate heading in slot @other has to be 'latlong' or 'xy'.\n   All analyses will be skipped!\n") 
 
if ((mk.resistance==TRUE | mk.complete==TRUE)  & (coords & !is.null(fric.raster))) 
  {
    message("- Landscape genetic analysis using resistance matrices...\n")  
 #   fr.raster<<-fric.raster
 #   gen.dist <<- gen.distance
    pmantel<-  readLines(paste(path,"pmantel.snw",sep=""))
    compl<-c(compl,pmantel)
  } 
  

if (mk.custom==TRUE){
  message("- Run customised snw file, custom.snw ...\n")
  custom<-readLines(paste(path,"custom.snw",sep=""))
  compl<-c(compl,custom)
}



footer.file<-readLines(paste(path,"footer.snw",sep=""))  
compl<-c(compl,footer.file)

#compl <- c(header.file, required, loaddata, mapping, popheterozygosity, footer.file)


rnwfile <- paste(fname,".rnw",sep="")
texfile <-  paste(fname,".tex",sep="") 


zz <- file(file.path(path.pgr,foldername,rnwfile), "w")
writeLines(compl,zz)
close(zz) 


#setwd(paste(path.pgr,foldername, sep="/"))
message(paste("Analysing data ...\n", sep=""))
#Sweave(paste(fname,".rnw",sep=""), output=paste(fname,".tex",sep=""), quiet=FALSE, driver=mydriver)
flush.console()
knit(input=rnwfile, output=texfile, quiet=TRUE, envir=pgr)

if (mk.pdf==TRUE)
{
message(paste("Creating pdf from: ",rnwfile," ...\n",sep=""))
knit2pdf(texfile, texfile)
message(paste("Finished.\nCheck ",fname,".pdf for results.\n", sep=""))
}

if (mk.Rcode) {
  message(paste("Creating R code from: ",rnwfile,"...\n"), sep="")
  rfile <-paste(fname,".R",sep="")
  purl(input=rnwfile, output=rfile)
#  Stangle(paste(fname,".rnw",sep=""))
}
    

message(paste("All files are available in the folder: \n",file.path(path.pgr, foldername),"\n",sep=""))

#reset working directory to previous
setwd(owd)
return(pgr$allresults)
}
