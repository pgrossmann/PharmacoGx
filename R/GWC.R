########################
## Benjamin Haibe-Kains
## All rights Reserved
## October 23, 2013
########################

#################################################
## Normalize and import all data from Connectivity Map
##
## inputs:	
##  x1: vector of effect sizes (e.g., fold change or t statitsics) for the first experiment
##  p1: vector of p-values for each corresponding effect size for the first experiment
##  x2: effect size (e.g., fold change or t statitsics) for the second experiment
##  p2: vector of p-values for each corresponding effect size for the second experiment
##  ...: additional parameters to be passed to corWeighted
##
## outputs, a list of items:
##			- 
##
## http://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient#Calculating_a_weighted_correlation
## http://www.mathworks.com/matlabcentral/fileexchange/20846-weighted-correlation-matrix
## F. Pozzi, T. Di Matteo, T. Aste, "Exponential smoothing weighted correlations", The European Physical Journal B, Vol. 85, No 6, 2012. DOI: 10.1140/epjb/e2012-20697-x
#################################################

`GWC` <- 
function (x1, p1, x2, p2, ...) {
  
  ## scaled weights
  p1 <- -log10(p1)
  p1 <- p1 / sum(p1, na.rm=TRUE)
  p2 <- -log10(p2)
  p2 <- p2 / sum(p2, na.rm=TRUE)
  w <- p1 + p2
  ## compute genome-wide connectivity score
  res <- corWeighted(x=x1, y=x2, w=w, ...)
	return(res)
}
