########################
## Benjamin Haibe-Kains
## All rights Reserved
## October 23, 2013
########################

#################################################
## Get fRMA normalized CMAP data from InSilicoDB
##
## 
#################################################

`getCMAP` <- 
function (standardization=c("none", "quantile"), gene=TRUE, nthread=1, verbose=FALSE) {

  # require(inSilicoDb2)

  InSilicoLogin(login="bhaibeka@gmail.com", password="747779bec8a754b91076d6cc1f700831")
  # inSilicoDb2::getCurationInfo(dataset="ISDB12026")
  if(verbose) { message("Downloading the Connectivity Map dataset from InSilicoDB") }
  esets <- inSilicoDb2::getDatasets(dataset="ISDB12026", norm="ORIGINAL", curation="24391", features="PROBE")
  InSilicoLogout()
  
  ## standardized of the expression data between the two microarray platforms
  
  ## gene centric expression
  for(i in 1:length(esets)) {
    esets[[i]] <- MetaGx::probeGeneMapping(platform=c("MISC", "GPL8300", "GPL96", "GPL97", "GPL570", "GPL1352"), ...)
  }
  
  ## merge the 2 microarray platforms
  
  
  
  
}

## End
