#' Bilby data set
#' 
#' This is a synthetic sample data set (in \code{\link{genind}} format) of
#' microsatellite data.
#' 
#' @name bilby
#' @format genlight object
#' @docType data
#' @author Bernd Gruber \email{(bernd.gruber@@canberra.edu.au}
#' @keywords datasets
#' data(bilby)
#' bilby
"bilby"



#' A simulated genind data set with spatial coordinates
#' 
#' This data set is used to demonstrate the use of the
#' \code{\link{landgenreport}} and \code{\link{genleastcost}} functions. It is
#' a simple spatial \code{genind} object with 100 individuals in 10 populations
#' and 20 loci with up to 20 alleles per loci.
#' 
#' 
#' @name landgen
#' @docType data
#' @author Bernd Gruber \email{(bernd.gruber@@canberra.edu.au}
#' @seealso \code{\link{landgenreport}}, \code{\link{genleastcost}}
#' @keywords datasets
#' @examples
#' data(landgen)
#' summary(landgen)
"landgen"


#' A genlight object created via the read.genetable functions [possum data set from Sarre et al. 2015]
#'
#' @name possums
#' @format genind object
#' @docType data
#' @author Bernd Gruber \email{(bernd.gruber@@canberra.edu.au}
#' @keywords datasets
#' @references 
#' Sarre, S.D., Aitken, N., Adamack, A.T., Macdonald, A.J., Gruber, B. & Cowan, P. (2014). Creating new evolutionary pathways through bioinvasion: The population genetics of brushtail possums in New Zealand. Molecular Ecology, 23, 3419-3433.
#' @examples
#' possums
"possums"


