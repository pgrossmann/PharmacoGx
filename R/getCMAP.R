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
function (std=c("combat", "quantile", "none"), gene=TRUE, nthread=1, verbose=FALSE) {

  # require(inSilicoDb2)
  # require(sva)

  std <- match.arg(std)
  
  InSilicoLogin(login="bhaibeka@gmail.com", password="747779bec8a754b91076d6cc1f700831")
  # inSilicoDb2::getCurationInfo(dataset="ISDB12026")
  if(verbose) { message("Downloading the Connectivity Map dataset from InSilicoDB") }
  platfs <- inSilicoDb2::getPlatforms(dataset="ISDB12026")
  esets <- inSilicoDb2::getDatasets(dataset="ISDB12026", norm="FRMA", curation="24391", features="PROBE")
  InSilicoLogout()
  
  ## merge esets
  eset <- MetaGx::platformMerging(esets=esets)
  pheno <- Biobase::pData(eset)
  
  ## standardized of the expression data between the two microarray platforms
  switch (std, 
    "combat" = {
      batch <- as.factor(pheno[ ,"chiptype"])
      mod <- model.matrix(~ as.factor(xptype), data=pheno)
      exprs(eset) <- t(sva::ComBat(dat=exprs(eset), batch=batch, mod=mod))
    },
    "quantile" = {
      ## compute the median profile of HT_HG-U133A samples
      myx <- rownames(pheno)[!is.na(pheno[ ,"chiptype"]) & pheno[ ,"chiptype"] == "HT_HG-U133A"]
      tt <- apply(exprs(eset)[ , myx, drop=FALSE], 2, sort, method="quick")
      qnormvec.hthgu133acmap <- apply(tt, 1, median, na.rm=TRUE)
      ## normalize the data
      exprs(eset) <- t(MetaGx::normQuant(A=t(exprs(eset)), ties=TRUE, normvector=qnormvec.hthgu133acmap))
    }
  )
  
  ## gene centric expression
  eset <- MetaGx::probeGeneMapping(eset=eset, platform="GPL96", metghod="jetset")
  
  return (eset)
}

## End
