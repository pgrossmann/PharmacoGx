########################
## Benjamin Haibe-Kains
## All rights Reserved
## October 23, 2013
########################

#################################################
## Extract the date and hour of hybridization from Affymetrix CEL files
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
`celfileDateHour` <-
function (filename) {
	h <- affyio::read.celfile.header(filename, info="full")
	#ddate <- grep("/", strsplit(h$DatHeader, " ")[[1]], value=TRUE)
	#ddate <- strsplit(ddate, split="/")[[1]]
	#CC <- ifelse(substr(ddate[3],1,1)=="9", "19", "20")
	if(length(h$ScanDate) > 0) {
	    h$ScanDate <- gsub(pattern="T", replacement=" ", x=h$ScanDate)
	    ddate <- strsplit(h$ScanDate, " ")[[1]]
    } else { ddate <- rep(NA, 2)}
    names(ddate) <- c("day", "hour")
	return(ddate)
}

## End
