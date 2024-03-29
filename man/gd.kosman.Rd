% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gd.kosman.r
\name{gd.kosman}
\alias{gd.kosman}
\title{Individual genetic distance calculation based on Kosman & Leonhard 2005}
\usage{
gd.kosman(population)
}
\arguments{
\item{population}{this is the \code{\link{genind}} object the analysis will
be based on.}
}
\value{
Returns a list with two distance matrices. The first (geneticdist)
contains pairwise individual genetic distances for each individual within a
population, the second (loci_used) provides the number of loci that were
used to calculate pairwise individual genetic distances between a pair of
individuals.
}
\description{
Calculates pairwise genetic distances between all individuals using the
individual genetic distance measure of Kosman and Leonard (2005). This
function is similiar to the dist.codom in the package mmod. The two
functions differ in their treatment of individuals with missing data.
dist.codom omits individuals from the calculation of pairwise individual
genetic distances while this function includes individuals with missing
data. This is done by simply calculating the mean individual pairwise
genetic distance over all loci for which there are values. Note that
depending on your computers capabilities, you may run into memory errors
when this function is used on datasets with large numbers of individuals
(>1000). Additionally, the time for this function to run can be lengthy.
Using a PC with a 3.5 GHz processor to calculate pairwise genetic distances
between individuals with 36 diploid loci it took 0.3 seconds for 100
individuals, 5 seconds for 500 individuals, 21 seconds for 1000 individuals,
84 seconds for 2000 individuals, and 194 seconds for 3000 individuals.
}
\examples{

\dontrun{
data(bilby)
popgenreport(bilby, mk.gd.kosman = TRUE, mk.pdf=FALSE)
}
#to get a pdf output you need to have a running Latex version installed on your system.
#popgenreport(bilby, mk.gd.kosman = TRUE, mk.pdf=TRUE)
}
\references{
Kosman E., Leonard K.J. 2005. Similarity coefficients for
molecular markers in studies of genetic relationships between individuals
for haploid, diploid, and polyploidy species. Molecular Ecology 14:415-424.
}
\seealso{
\code{\link{popgenreport}}
}
\author{
Aaron Adamack, aaron.adamack@canberra.edu.au
}
