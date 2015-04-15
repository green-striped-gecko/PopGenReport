gd.kosman <- function(population){
  # check to see if the function was passed a genind object
  if (class(population) != "genind") {
    message("You did not provide a valid genind object! Script stopped!")
    return
  }
  
  # getting the number of loci
  n <- length(population@loc.names) 
  
  # getting the ploidy of an animal....
  ploidy <- population@ploidy 
  
  # this calculates the manhattan distance between each individual and adjusts for ploidy..
  matrices <- lapply(seploc(population), function(l) as.matrix(dist(l@tab, "manhattan")/ploidy))
  
  # if the value is missing, mark it with a 1, if it is real, mark it 0
  missing <- lapply(matrices, function(m) ifelse(is.na(m), 1, 0))
  
  # This is going to replace missing numbers with 0. 
  replaced <- lapply(matrices, function(m) ifelse(is.na(m),0,m))
  
  
  loci.used<-(n-Reduce("+", missing))
  colnames(loci.used) <- population@ind.names
  rownames(loci.used) <- population@ind.names
  
  # This sums the values across lists and then divides by the number of loci compared less loci with missing numbers  
  d.fast<-(Reduce("+", replaced)/loci.used)
  colnames(d.fast) <- population@ind.names
  rownames(d.fast) <- population@ind.names
  
  # clean up matrices for export
  d.fast[upper.tri(d.fast, diag = FALSE)] <- NA
  loci.used[upper.tri(loci.used, diag = FALSE)] <- NA
  kosman.out <- list(geneticdist = d.fast, loci_used = loci.used)
  return(kosman.out)
}
