#' This is the main function of the package. It analyses an object of class
#' \code{\link{genind}} and then creates a report containing the results of the
#' analysis. There are several routines that can be optionally included in the
#' analysis and there are multiple output options including a PDF with the
#' report, R-code and an object (\code{fname.results}) containing all of the
#' results, which can be used for further analyses.
#' 
#' This function is used to analyse population genetic data. The main idea is
#' to provide a framework for analysing microsatellite and also SNP genetic
#' data (if not too many loci, say below 1000) using a mix of existing and new
#' functions. The function works on an object of class \code{genind}. There are
#' several ways to convert data into a \code{\link{genind}} object using
#' existing functions provided by the \code{adegenet} package (
#' \code{\link{import2genind}},
#' \code{\link{df2genind}},\code{\link{read.fstat}},
#' \code{\link{read.structure}}, \code{\link{read.genetix}}
#' ,\code{\link{read.genepop}}) or refer to \code{read.genetable} how to import
#' data from an EXCEL (csv) document. The function performs a number of
#' different genetic analyses (e.g. counts of indivuals and alleles across
#' sub-populations, tests for heterozygosity and Hardy-Weinberg Equilibrium,
#' differentiation statistics Fst, G'st, Jost's D, and genetic distance between
#' individuals and populations), with users having the option to select which
#' analysis routines are included in the report. To select a routine, the user
#' simply turns on a switch e.g. mk.map=TRUE returns a map with the sampling
#' location for each individual (if coordinates are provided). \cr Coordinates
#' need to specified within the genind object. As a standard genind object does
#' not require spatial coordinates, we extended it by using the \code{other}
#' slot in the genind object. The easiest way to provide spatial coordinates is
#' to use the read.genetable function and use the \code{lat}, \code{long} or
#' \code{x}, \code{y} arguments for WGS1984 projected data or mercator
#' projected data respectively. To calculate distances the data are internally
#' reprojected using the \code{\link{Mercator}} function in package
#' \code{\link{dismo}}), which is the projection used by google maps. Or you
#' can add data manually to your \code{genind} object using the mentioned (e.g.
#' \code{genindobject@other$latlong <- yourlatlong data} or
#' \code{genindobject@other$xy <- your_xy_data}). If you have your data in a
#' different projection you need to reproject them into either WGS1984 or the
#' google maps Mercator projection. If you use a different projection distance
#' calculation may be wrong and probably the map will not be correct. See the
#' manual for an example how to project and add spatial coordinates to your
#' genetic data.\cr Names for alleles (\code{genindobject@loc.names}) are
#' truncated if longer than six characters. If truncated Captial letters linked
#' by a hyphen are added to guarentee they are unique. You can rename them by
#' providing new names by accessing the \code{genind@loc.names} slot prior to
#' running \code{popgenreport}.\cr Note that the popgenreport function can take
#' a long time to run if the options mk.complete, mk.gd.kosman, or mk.gd.smouse
#' are set to \code{TRUE}. For example, running popgenreport with
#' \code{mk.complete=TRUE} on a dataset with 500 individuals with 36 loci will
#' take 14 to 15 minutes on a PC with a 3.5 Ghz processor and nearly 3 hours
#' for a dataset with ~3200 individuals.
#' 
#' 
#' @param cats this is the \code{\link{genind}} object the analysis will be
#' based on.
#' @param mk.counts switch is to provide overview counts of the number of
#' individuals sampled, numbers of individuals and alleles sampled per
#' sub-population, number of alleles per locus, mean number of alleles per
#' locus and the percentatge of missing data.
#' @param mk.map switch to produce a map with the sampling location of each
#' individual marked. This switch requires individual coordinates (latitudes
#' and longitudes in WGS1984) be provided (under cats@other$latlong or see
#' \code{read.genetable} on how to import them from a table of genetic data).
#' An error message will be generated if you turn this routine on, but do not
#' provide the coordinates in the right format. If the coordinates are provided
#' in a seperate file, they must be attached to the genind object in the slot
#' \cr\code{yourgenindobject@other$latlong <- yourlatlongdata}.\cr
#' \code{yourlatlongdata} needs to be a data frame that has the same number and
#' order of individuals per row as the population genetic data. Note that an
#' internet connection is required to connect to the Google Maps server which
#' provides the basemap for this routine.
#' @param maptype Defines the type of map. Default is 'satellite'. Other
#' options are: 'roadmap', 'mobile', 'terrain', 'hybrid'.
#' @param mapdotcolor Color of dots for each individual on the map. Default is
#' 'blue'.
#' @param mapdotsize Size of dots for each individual. Default is 1.
#' @param mapdotalpha Transparency of dots. 1 is invisible, 0 is no
#' transparency. Default is 0.4.
#' @param mapdottype Defines the type of the symbol. For explanation see pch
#' under \code{\link{par}}. Default is 19 - a filled circle.
#' @param mapzoom Zoom level of the map. If not specified the default zoom of
#' Google maps are used. Please be aware if you set the zoom level to high, the
#' map may not show all sample locations.
#' @param mk.locihz switch to test for population heterozygosity
#' @param mk.hwe switch to test for Hardy-Weinberg equilibrium for each loci
#' and population
#' @param mk.fst switch to calculate Fst values for each loci and pairwise Fst
#' (Nei's 1973) over subpopulations
#' @param mk.gd.smouse Individual pairwise genetic distances based on Smouse
#' and Peakall (1999). Refer to \code{gd_smouse}. Spatial coordinates need to
#' be provided to be able to run this analysis.
#' @param mk.gd.kosman Individual pairwise genetic distances based on Kosman &
#' Leonhard (2005). Refer to \code{gd_kosman}. Spatial coordinates need to be
#' provided to be able to run this analysis.
#' @param mk.pcoa Principal component analysis following Jombart et al. 2009.
#' Spatial coordinates need to be provided to be able to run this analysis.
#' Refer to vignettes within \code{adegenet}.
#' @param mk.spautocor Spatial autocorrelation analysis following Smouse &
#' Peakall 1999. Spatial coordinates need to be provided to be able to run this
#' analysis. Refer to \code{spautocor} for more information.
#' @param mk.allele.dist switch to look at allele distributions by loci and
#' subpopulation
#' @param mk.null.all check for null alleles
#' @param mk.allel.rich calculation of allelic richness
#' @param mk.differ.stats switch to look at population differentiation
#' statistics (Nei's Gst, Hedrick's Gst, and Jost's D)
#' @param mk.custom edit custom.snw to include your own function to a report.
#' @param fname filename for the output files. Defauts to PopGenReport. Note
#' that using a filename which includes a space in the name will result in the
#' filename for each figure being printed out in the PDF report for each
#' figure. Replacing the space with an underscore should prevent this from
#' happening.
#' @param foldername name of folder, where files are stored. Defaults to
#' 'results'
#' @param path.pgr Folder where the output files are stored. Defaults to the
#' temporary directory (\code{tempdir()}). If you want to store the output in
#' another directory, simply provide the path here. e.g.
#' \code{path.pgr=getwd()} saves it in your current working directory.
#' @param mk.Rcode switch to get the full R script that is used to generate the
#' report. A great way to get a very detailed insight on the kind of analysis
#' and also an easy way to generate a script which you can customize for your
#' analytical needs.
#' @param mk.complete switch to create a full report using all of the routines
#' (all switches are set to \code{TRUE}, except \code{mk.subgroups}).
#' @param mk.pdf switch to create a shiny pdf output. You need a working
#' \bold{latex} version running on your system (e.g. MikTex (Windows) or
#' Texmaker (Linux, MacOSX). For more information how to install latex on your
#' system refer to the \url{www.popgenreport.org} and to the manuals of the
#' \code{\link{knitr}} package and its manuals.
#' @return The function returns an object (e.g. res) that has all of the
#' results produced by this function in it. The structure of the object can be
#' accessed via \code{str(res)}. The main slots in this object (if you ran a
#' full report) are:\cr \code{dataoverview, PopHet, Alleledist, Fst,
#' HsHtdifferentiate, HWEresults,} \cr\code{subgroups, GDKosman, GDSmouse}
#' 
#' Additional ouput is provided in the form of a PDF (if mk.pdf=TRUE),which
#' will be saved to the specified subfolder (via foldername) in your current
#' working directory, and maps and figures which will be placed in this folder
#' as well. This folder will be generated automatically in your current working
#' directory. If you do not specify a working directory via \code{path.pgr}
#' then the temporary working directory of R will be used (\code{tempdir()}).
#' If \code{mk.Rcode=T} is set, an R file named fname.R will be saved to your
#' specified subfolder.
#' @author Aaron Adamack & Bernd Gruber, aaron.adamack@@canberra.edu.au,
#' bernd.gruber@@canberra.edu.au
#' @seealso \code{\link{adegenet}}, \code{\link{pegas}}, \code{\link{mmod}}
#' @references Kosman E., Leonard K.J. 2005. Similarity coefficients for
#' molecular markers in studies of genetic relationships between individuals
#' for haploid, diploid, and polyploidy species. Molecular Ecology 14:415-424
#' 
#' Peakall R., Smouse P. 2012. GenAlEx 6.5: Genetic analysis in Excel.
#' Population genetic software for teaching and research - an update.
#' Bioinformatics 28:2537-2539
#' 
#' %%Adamack & Gruber (2012)
#' @examples
#' 
#' #not run:
#' #data(bilby) # a generated data set
#' #res <- popgenreport(bilby, mk.counts=TRUE, mk.map=TRUE, mk.pdf=FALSE)
#' #check results via res or use created tables in the results folder.
#' 
#' ### RUN ONLY with a working Latex version installed
#' # res <- popgenreport(bilby, mk.counts=TRUE, mk.map=TRUE, mk.pdf=TRUE, path.pgr="c:/temp")
#' # for a full report in a single pdf set mk.complete to TRUE
#' # res <- popgenreport(bilby, mk.complete=TRUE)
#' @importFrom ggplot2 ggplot 
#' @importFrom dismo Mercator 
#' @importFrom xtable xtable 
#' @importFrom RgoogleMaps GetMap.bbox PlotOnStaticMap 
#' @importFrom calibrate textxy 
#' @importFrom genetics HWE.test HWE.chisq HWE.exact 
#' @importFrom ade4 dudi.pca s.class add.scatter 
#' @importFrom rgdal project 
#' @import mmod 
#' @import R.utils 
#' @import gap 
#' @import knitr 
#' @import adegenet 
#' @importFrom raster raster extent rasterize crs values buffer values<- crs<- plot
#' @importFrom GGally ggpairs
#' @importFrom graphics abline axis box hist image lines par points text
#' @importFrom grDevices col2rgb colorRampPalette rainbow rgb 
#' @importFrom methods as 
#' @importFrom stats as.dist dist ecdf lm median quantile rbinom resid runif sd weighted.mean 
#' @importFrom utils combn flush.console read.csv 
#' @export

popgenreport <- function(cats=NULL,
                          
                            mk.counts=TRUE,   # this switch is to provide a population overview
                            mk.map=FALSE,        # this switch is for the map
                            maptype="satellite",
                            mapdotcolor ="blue",
                            mapdotsize=1,
                            mapdotalpha=0.4,
                            mapdottype=19 ,
                            mapzoom=NULL,
                          
#c("roadmap",
#"mobile",
#"satellite",
#"terrain",
#"hybrid",
#"mapmaker-roadmap",
#"mapmaker-hybrid")",
#                          
                          
                          mk.locihz=FALSE,     # this switch is to test for population heterozygosity
                          mk.hwe=FALSE,   # this switch is for population wide HWE

                          mk.fst=FALSE,        # this switch is to run FST tests on the full population
                          mk.gd.smouse=FALSE,   # this switch is to run the Smouse and Peakall genetic distances
                          mk.gd.kosman=FALSE,   # this switch is to run the Kosman and Leonard genetic distances

                          mk.pcoa=FALSE,
                          mk.spautocor=FALSE,
                          mk.allele.dist=FALSE, # this switch it to look at allele distributions by loci and pop
                          mk.null.all=FALSE,
                          mk.allel.rich=FALSE,

                          mk.differ.stats=FALSE ,     # this switch is to look at population differentiation statistics (Fst, Gst, etc)
                          mk.custom = FALSE,
                          fname="PopGenReport",
                          foldername="results",
                          path.pgr=NULL,
                          mk.Rcode=FALSE,       # make the code that was ran available as an R file
                          mk.complete=FALSE,    # create a full report)  
                          mk.pdf=TRUE)
{
  if (class(cats)!="genind") {stop("You did not provide a valid genind object! Script stopped!\n")}
  
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
  
### adjusting causes error with Adegenet 2.0 and it is controlled for when creating the object!!  
  #cut down length of loci names to  6  and make sure they are unique
  #locNames(cats) <- substr(locNames(cats),1,6)   
  #if (length(unique(locNames(cats)))!= length(locNames(cats))) 
  #{locNames(cats) <- paste(1:length(locNames(cats)),"-",substr(locNames(cats),1,4), sep="")
  #cat("Loci names were not unique and therefore adjusted.\n")
  #}
  #levels(cats@loc.fac) <- locNames(cats)  #make sure levels and factors are the same
#check if pop.names, indnames are unique!!!!
#adjust if necessary and issue a notification
# if (length(unique(indNames(cats)))!=length(indNames(cats))) 
#   {indNames(cats) <- paste(1:length(indNames(cats)),"-",substr(indNames(cats),1,8),sep="")
#   cat("Individual names were not unique and therefore adjusted.\n")
#   }
# 
# 
#   
# if (length(unique(popNames(cats)))!=length(popNames(cats))) 
#   {
#   popNames(cats) <- paste(1:length(popNames(cats)),"-",substr(popNames(cats),1,6),sep="")
#   cat("Subpopulation names were not unique and therefore adjusted.\n")
#   }
# 
#  
    #set directory where to save a file, defaults to tempdir (follow R policy)
  if (is.null(path.pgr)) 
  {
  path.pgr <- tempdir()

  }
  

  
#  setwd(path.pgr)
  
  #create a foldername folder if not existing...
  dirfiles <- list.dirs(path=path.pgr, recursive=FALSE)
  if (!(tolower (file.path(path.pgr,foldername))) %in% tolower(dirfiles)) {
    dir.create(file.path(path.pgr,foldername))
    cat("There is no ",foldername, " folder. I am trying to create it; \notherwise please create the folder manually. \n")
  }
  owd <-getwd()
  setwd(file.path(path.pgr, foldername))
  # conversion of lat longs to google map data (Mercator (dismo) wants to have long lat)
  

  coords=FALSE
  if (is.null(cats$other$latlong) & is.null(cats@other$xy)) coords=FALSE else {
  if (!is.null(cats@other$latlong)) cats@other$xy <- Mercator(cats@other$latlong[,c(2,1)])
  if (!is.null(cats@other$xy)) cats@other$latlong <- Mercator(cats@other$xy, inverse=TRUE)[,c(2,1)]  
  
  if ((nrow(cats@other$latlong) == length(indNames(cats))) & (nrow(cats@other$xy) == length(indNames(cats)) )) coords=TRUE
  


  } 
  # give cats a filename that can be seen in the snw chunks
  
  
  cats@other$filename<- fname
  cats@other$foldername<-foldername
  cats@other$path <- path.pgr
  #determine the type of map
  if ((mk.map==TRUE | mk.complete) & coords) 
  {
  cats@other$maptype=maptype
  cats@other$mapdotcolor =mapdotcolor
  cats@other$mapdotsize=mapdotsize
  cats@other$mapdotalpha=mapdotalpha
  cats@other$mapdottype=mapdottype
  cats@other$mapzoom=mapzoom
  
  }  
  
###################################
##### create a new environment to run knitr in it
 #pgr <- new.env(parent=.GlobalEnv)
 pgr <- new.env()
 assign("cats",cats,envir=pgr)
 #sumcats <- summary(cats, verbose=FALSE)
 #assign("sumcats",sumcats, envir=pgr )
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
if (is.null(path)) {stop("Could not find snw files in the PopGenReport library folder. Please check if the package is installed correctly (e.g.  installed.packages()[\"PopGenReport\",2]). \n")}
  #for testing:
  #path <- "d:\\bernd\\R\\popgenreport\\inst\\swchunks\\"
  #path<- "C:\\Aaron files\\popgenreport098\\PopGenReport_0.98\\PopGenReport\\swchunks\\"
  header.file <- readLines(paste(path,"header.snw",sep=""))
  required<- readLines(paste(path,"required.snw",sep=""))
  compl<-c(header.file,required) 
  
  cat("Compiling report...\n")
  if(mk.counts | mk.complete){
    cat("- General summary...\n")
    overview<-readLines(paste(path,"counts.snw",sep=""))
    compl<-c(compl,overview)
  }
  if (coords==FALSE) cat(" - No valid coordinates were provided. \n   Be aware you need to provide a coordinate (or NA) for each individual\n   and the coordinate heading in slot @other has to be 'latlong' or 'xy'.\n   Some of the analyses require coordinates and will be skipped!\n") 
  if ((mk.map==TRUE | mk.complete) & coords){
    cat("- Map of individuals...\n")  
    mapping<-  readLines(paste(path,"map.snw",sep=""))
    compl<-c(compl,mapping)
  }
  
  if (mk.locihz | mk.complete){
    cat("- Statistics on population heterogeneity ...\n")  
    popheterozygosity <- readLines(paste(path,"locihz.snw",sep=""))
    compl<-c(compl,popheterozygosity)
  }
  
  if (mk.allele.dist | mk.complete){
    cat("- Allelic distances ...\n")  
    numloci<-length(cats@loc.n.all)
    alleledistn <- readLines(paste(path,"allele.dist.snw",sep=""))
    compl<-c(compl,alleledistn)
}
if (mk.fst| mk.complete){
     cat("- Pairwise Fst ...\n")  
  popfst<-readLines(paste(path,"fst.snw",sep=""))
  compl<-c(compl,popfst)
}
if (mk.null.all | mk.complete){
     cat("- Checking for null alleles ...\n")  
     null.stat<-readLines(paste(path,"null.all.snw",sep=""))
  compl<-c(compl,null.stat)
}
if (mk.allel.rich | mk.complete){
     cat("- Allelic richness ...\n")  
     all.stat<-readLines(paste(path,"allel.rich.snw",sep=""))
  compl<-c(compl,all.stat)
}
if (mk.differ.stats | mk.complete){
     cat("- Pairwise differentiations ...\n")  
     diff.stat<-readLines(paste(path,"differ.stats.snw",sep=""))
  compl<-c(compl,diff.stat)
}

if (mk.hwe | mk.complete){
  cat("- Test for Hardy-Weinberg-Equilibrium ...\n") 
  cat("  !! You may get warnings when running HWE tests, if the test is based\n")
  cat(" on an entry in the chi-square table which is less than five.!! \n")
  popHWEll<-readLines(paste(path,"hwe.snw",sep=""))
  compl<-c(compl,popHWEll)
}

if ((mk.gd.kosman==TRUE | mk.complete) & coords){
  cat("- Kosman & Leonard 2005 genetic distances...\n")
  kosman<-readLines(paste(path,"gd.kosman.snw",sep=""))
  compl<-c(compl,kosman)
}

if ((mk.gd.smouse==TRUE | mk.complete) & coords){
  cat("- Smouse & Peakall 1999 genetic distances...\n")
  smouse<-readLines(paste(path,"gd.smouse.snw",sep=""))
  compl<-c(compl,smouse)
}
if ((mk.spautocor==TRUE | mk.complete) & coords){
  cat("- Spatial autocorrelation following Smouse & Peakall 1999 ...\n")
  spa<-readLines(paste(path,"spautocor.snw",sep=""))
  compl<-c(compl,spa)
}
if (mk.pcoa==TRUE | mk.complete){
  cat("- Principal coordinate analysis following Jombart et al. 2009...\n")
  pca<-readLines(paste(path,"pcoa.snw",sep=""))
  compl<-c(compl,pca)
}

if (mk.custom==TRUE){
  cat("- Run customised snw file, custom.snw ...\n")
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
cat(paste("Analysing data ...\n", sep=""))
#Sweave(paste(fname,".rnw",sep=""), output=paste(fname,".tex",sep=""), quiet=FALSE, driver=mydriver)
flush.console()
knit(input=rnwfile, output=texfile, quiet=TRUE, envir=pgr)

if (mk.pdf==TRUE)
{
cat(paste("Creating pdf from: ",rnwfile," ...\n",sep=""))
knit2pdf(texfile, texfile)
cat(paste("Finished.\nCheck ",fname,".pdf for results.\n", sep=""))
}

if (mk.Rcode) {
  cat(paste("Creating R code from: ",rnwfile,"...\n"), sep="")
  rfile <-paste(fname,".R",sep="")
  purl(input=rnwfile, output=rfile)
#  Stangle(paste(fname,".rnw",sep=""))
}
    

cat(paste("All files are available in the folder: \n",file.path(path.pgr, foldername),"\n",sep=""))

#reset working directory to previous
setwd(owd)
return(pgr$allresults)
}
