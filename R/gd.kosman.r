#' Individual genetic distance calculation based on Kosman & Leonhard 2005
#' 
#' Calculates pairwise genetic distances between all individuals using the
#' individual genetic distance measure of Kosman and Leonard (2005). This
#' function is similiar to the dist.codom in the package mmod. The two
#' functions differ in their treatment of individuals with missing data.
#' dist.codom omits individuals from the calculation of pairwise individual
#' genetic distances while this function includes individuals with missing
#' data. This is done by simply calculating the mean individual pairwise
#' genetic distance over all loci for which there are values. Note that
#' depending on your computers capabilities, you may run into memory errors
#' when this function is used on datasets with large numbers of individuals
#' (>1000). Additionally, the time for this function to run can be lengthy.
#' Using a PC with a 3.5 GHz processor to calculate pairwise genetic distances
#' between individuals with 36 diploid loci it took 0.3 seconds for 100
#' individuals, 5 seconds for 500 individuals, 21 seconds for 1000 individuals,
#' 84 seconds for 2000 individuals, and 194 seconds for 3000 individuals.
#' 
#' 
#' @param population this is the \code{\link{genind}} object the analysis will
#' be based on.
#' @return Returns a list with two distance matrices. The first (geneticdist)
#' contains pairwise individual genetic distances for each individual within a
#' population, the second (loci_used) provides the number of loci that were
#' used to calculate pairwise individual genetic distances between a pair of
#' individuals.
#' @author Aaron Adamack, aaron.adamack@@canberra.edu.au
#' @seealso \code{\link{popgenreport}}
#' @references Kosman E., Leonard K.J. 2005. Similarity coefficients for
#' molecular markers in studies of genetic relationships between individuals
#' for haploid, diploid, and polyploidy species. Molecular Ecology 14:415-424.
#' @examples
#' 
#' \dontrun{
#' data(bilby)
#' popgenreport(bilby, mk.gd.kosman = TRUE, mk.pdf=FALSE)
#' }
#' #to get a pdf output you need to have a running Latex version installed on your system.
#' #popgenreport(bilby, mk.gd.kosman = TRUE, mk.pdf=TRUE)
#' @export
gd.kosman <- function(population){
  # check to see if the function was passed a genind object
  if (!is(population,"genind")) {
    stop("You did not provide a valid genind object! Script stopped!")
  }
  
  # getting the number of loci
  n <- length(locNames(population)) 
  
  # getting the ploidy of an animal....
  uniqueploidy<-unique(population@ploidy)
  if(length(uniqueploidy)==1){
    ploidy<-uniqueploidy
  } else if(uniqueploidy>1){
    message("Your data set has multiple ploidies, please separate loci by ploidy")
    stop("Script stopped!")

  } else if(uniqueploidy<=0){
    stop("Your dataset has an invalid ploidy (ploidy<=0). Script stopped!")

  }
  
  # this calculates the manhattan distance between each individual and adjusts for ploidy..
  matrices <- lapply(seploc(population), function(l) as.matrix(dist(l@tab, "manhattan")/(2*ploidy)))
  
  # if the value is missing, mark it with a 1, if it is real, mark it 0
  missing <- lapply(matrices, function(m) ifelse(is.na(m), 1, 0))
  
  # This is going to replace missing numbers with 0. 
  replaced <- lapply(matrices, function(m) ifelse(is.na(m),0,m))
  
  
  loci.used<-(n-Reduce("+", missing))
  colnames(loci.used) <- indNames(population)
  rownames(loci.used) <- indNames(population)
  
  # This sums the values across lists and then divides by the number of loci compared less loci with missing numbers  
  d.fast<-(Reduce("+", replaced)/loci.used)
  colnames(d.fast) <- indNames(population)
  rownames(d.fast) <- indNames(population)
  
  # clean up matrices for export
  d.fast[upper.tri(d.fast, diag = FALSE)] <- NA
  loci.used[upper.tri(loci.used, diag = FALSE)] <- NA
  kosman.out <- list(geneticdist = d.fast, loci_used = loci.used)
  return(kosman.out)
}
