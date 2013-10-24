########################
## Benjamin Haibe-Kains
## All rights Reserved
## October 23, 2013
########################

#################################################
## Extract the type of GeneChip from Affymetrix CEL files
##
## inputs:	
##      - 
##
## outputs, a list of 2 items:
##			- 
##
## Notes:	duration is not taken into account as only 4 perturbations lasted 12h, the other 6096 lasted 6h
#################################################

## courtesy of Matthew McCall
`celfileChip` <-
function (filename) {
  	require(affyio)
  	h <- affyio::read.celfile.header(filename, info="full")
  	return(as.character(h$cdfName))
  }

## End
