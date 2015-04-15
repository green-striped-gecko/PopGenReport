\name{allele.dist}
\alias{allele.dist}
\title{
Counts and visualises allele frequencies across loci and subpopulations
}
\description{
Counts the number of observations for each combination of allele variant and subpopulation for each locus. Calculates relative allele proportions for each subpopulations and then produces a heatmap using that data. 
}
\usage{
allele.dist(population, mk.figures = TRUE)
}
\arguments{
  \item{population}{this is the \code{\link{genind}} object the analysis will be based on }
  \item{mk.figures}{ if set to FALSE no figures are plotted. Default is TRUE.
}
}
\value{
Produces heatmaps of the relative allele frequencies for each subpopulation at each locus and returns a list containing the counts (count) for each combination of allele and subpopulation and the relative frequencies of alleles by subpopulation (frequency) for each locus. The color bars on the heatmaps shows the relative frequency of an allele within a subpopulation for a locus while the histogram gives an overall count for the number of combinations of allele and subpopulation with a relative frequency. 
}
\author{
Aaron Adamack, aaron.adamack@canberra.edu.au
}
\seealso{
\code{\link{popgenreport}}
}
\examples{
# data(bilby)
 #here we use only the first 50 individuals to speep up the example
# popgenreport(bilby, mk.allele.dist=TRUE, mk.pdf=FALSE)
 
#to get a pdf output you need to have a running Latex version installed on your system.
#popgenreport(bilby, mk.allele.dist=TRUE, mk.pdf=TRUE)
}
