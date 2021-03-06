########################
## Benjamin Haibe-Kains
## All rights Reserved
## October 23, 2013
########################

########################################################
## Get fRMA or RMA normalized CMAP data from InSilicoDB
##
## 
########################################################

`getCMAP` <- 
function (std=c("combat", "quantile", "none"), gene=TRUE, normalization=c("rma", "frma"), verbose=FALSE) {

  std <- match.arg(std)
  normalization <- match.arg(normalization)
  esets <- list()
  
  if (normalization == "frma") {
    inSilicoDb::InSilicoLogin(login="bhaibeka@gmail.com", password="747779bec8a754b91076d6cc1f700831")
    # inSilicoDb::getCurationInfo(dataset="ISDB12026")
    if (verbose) { message("Downloading the Connectivity Map dataset from InSilicoDB") }
    platfs <- inSilicoDb::getPlatforms(dataset="ISDB12026")
    esets <- inSilicoDb::getDatasets(dataset="ISDB12026", norm="FRMA", curation="24805", features="PROBE")
    inSilicoDb::InSilicoLogout()
  }
  if (normalization == "rma"){
    datavector <- RMAnormalizationCMAP(gene=gene)
    data.cmap <- ExpressionSet(t(datavector$data.cmap))
    pData(data.cmap) <- datavector$sampleinfo.cmap
    fData(data.cmap) <- datavector$annot.cmap
    esets <- data.cmap
#     for(i in 1:length(unique(pData(data.cmap)[,"chiptype"]))) {
#       chipData <- pData(data.cmap)[pData(data.cmap)[,"chiptype"] == unique(pData(data.cmap)[,"chiptype"])[i],]
#       expression <- exprs(data.cmap)[,rownames(chipData)]
#       eset <- ExpressionSet(expression)
#       fData(eset) <- fData(data.cmap)
#       pData(eset) <- chipData
#       esets[[i]] <- eset  
#     }
  }
  ## merge esets
  if (verbose) { message("Merging CMAP1 and CMAP2") }
  eset <- MetaGx::platformMerging(esets=esets)
  Biobase::pData(eset)[Biobase::pData(eset) == "" | Biobase::pData(eset) == "NA"] <- NA
  ## check column format
  ############ NEED TO FIX THE DOTS ########################
  Biobase::pData(eset)[ , "concentration_M"] <- as.numeric(Biobase::pData(eset)[ , "concentration_M"])
  Biobase::pData(eset)[ , "duration_h"] <- as.numeric(Biobase::pData(eset)[ , "duration_h"])
  
  ## standardized of the expression data between the two microarray platforms
  if (std != "none" && verbose) { message("Standardizing the data between CMAP1 and CMAP2") }
  switch (std, 
    "combat" = {
      batch <- as.factor(Biobase::pData(eset)[ ,"chiptype"])
      mod <- model.matrix(~ as.factor(xptype), data=Biobase::pData(eset))
      exprs(eset) <- sva::ComBat(dat=exprs(eset), batch=batch, mod=mod)
    },
    "quantile" = {
      ## compute the median profile of HT_HG-U133A samples
      myx <- rownames(Biobase::pData(eset))[!is.na(Biobase::pData(eset)[ ,"chiptype"]) & Biobase::pData(eset)[ ,"chiptype"] == "HT_HG-U133A"]
      tt <- apply(Biobase::exprs(eset)[ , myx, drop=FALSE], 2, sort, method="quick")
      qnormvec.hthgu133acmap <- apply(tt, 1, median, na.rm=TRUE)
      ## normalize the data
      exprs(eset) <- t(MetaGx::normQuant(A=t(exprs(eset)), ties=TRUE, normvector=qnormvec.hthgu133acmap))
    }
  )
    
  ## gene centric expression
  if (gene && normalization=="frma") {
    if (verbose) { message("Gene centric data") }
    eset <- MetaGx::probeGeneMapping(eset=eset, platform="GPL96", method="jetset")
  }
  if (gene && normalization=="rma"){
    colnames(fData(eset))[which(colnames(fData(eset))=="jetset.EntrezID")] <- "ENTREZID"
    colnames(fData(eset))[which(colnames(fData(eset))=="jetset.symbol")] <- "SYMBOL"
  }

  ## get cell line information
  # celline <- sort(unique(Biobase::pData(eset)[ , "cell_id"]))
  # celline <- read.csv(file.path(system.file("extdata", package="PharmacoGx"), "cell_line_matching_all.csv"), stringsAsFactors=FALSE)
  # rownames(tissue.type) <- tissue.type[ , 1]
  # celline <- cbind("tissue.type"=tissue.type[match(celline, tissue.type[ , "cell_id"]), "tissue.type"], celline)
  celline <- cbind("cell_id"=c("HL60", "MCF7", "PC3", "SKMEL5", "ssMCF7"),
    "cell_id_standard"=c("HL-60", "MCF7", "PC-3", "SK-MEL-5", "ssMCF7"),
    "tissue.type"=c("haematopoietic_and_lymphoid_tissue", "breast", "prostate", "skin", "breast"))
  rownames(celline) <- as.character(celline[ , "cell_id"])
  
  ## get drug information
  tt <- Biobase::pData(eset)[ , c("drug_id", "cmap_name")]
  tt <- tt[!is.na(tt[ , "drug_id"]) & !duplicated(tt[ , "drug_id"]), , drop=FALSE]
  tt <- tt[!is.na(tt[ , "drug_id"]) & tt[ , "drug_id"] != "NA", , drop=FALSE]
  rownames(tt) <- as.character(tt[ , "drug_id"])
  druginfo <- tt
  ## ChemBank identifiers (CBID) provided by Justin Lamb
  tt <- read.csv(file.path(system.file("extdata", package="PharmacoGx"), "cmap_CBID_Lamb.csv"), stringsAsFactors=FALSE)
  ## WARNING: the cmap id do not correspond to the original annotations
  # rownames(tt) <- paste("drug_cmap", tt[ , "cmap_name_id"], sep=".")
  myx <- match(druginfo[ , "cmap_name"], tt[ , "cmap_name"])
  druginfo <- data.frame(druginfo, "ChemBank.id"=tt[myx, "CBID"], stringsAsFactors=FALSE)
  ## gather drug information from the CMAP website
  # http://www.broadinstitute.org/cmap/cmap_instances_02.xls
  tt <- gdata::read.xls(file.path(system.file("extdata", package="PharmacoGx"), "cmap_instances_02.xls"), stringsAsFactors=FALSE)
  tt <- tt[!is.na(tt[ , "cmap_name"]) & !duplicated(tt[ , "cmap_name"]), c("cmap_name", "vendor", "catalog_number", "catalog_name")]
  myx <- match(druginfo[ , "cmap_name"], tt[ , "cmap_name"])
  druginfo <- data.frame(druginfo, tt[myx, -1, drop=FALSE], stringsAsFactors=FALSE)
  
  ## get concentrations for each drugs
  tt <- table(Biobase::pData(eset)[ , "drug_id"], Biobase::pData(eset)[ , "concentration_M"])
  tt <- apply(tt, 1, function (x) {
    return (paste(names(x)[!is.na(x) & x > 0], collapse=", "))
  })
  drugconc <- cbind("drug_id"=rownames(druginfo), "concentration"=tt[rownames(druginfo)])
  
  ## get duration
  tt <- table(Biobase::pData(eset)[ , "drug_id"], Biobase::pData(eset)[ , "duration_h"])
  tt <- apply(tt, 1, function (x) {
    return (paste(names(x)[!is.na(x) & x > 0], collapse=", "))
  })
  drugduration <- cbind("drug_id"=rownames(druginfo), "duration"=tt[rownames(druginfo)])  
  
  return (list("expression"=eset, "cell.line"=celline, "drug"=list("info"=druginfo, "concentration"=drugconc, "duration"=drugduration)))
}

## End
