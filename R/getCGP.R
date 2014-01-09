########################
## Benjamin Haibe-Kains
## All rights Reserved
## October 23, 2013
########################

#################################################
## Get fRMA normalized CGP data from InSilicoDB
##
## 
#################################################

`getCGP` <- 
function (gene=TRUE, tmpdir="tmp", delete.tmpdir=FALSE, cosmic.annotation=FALSE, cosmic.version="v67_241013", replicates=c("last", "first", "all", "mean", "median"), verbose=FALSE) {

  replicates <- match.arg(replicates)
  
  ## create directories for temporary files
  if(!file.exists(tmpdir)) { dir.create(tmpdir, showWarnings=FALSE, recursive=TRUE) }
  
  badchars <- "[\xb5]|[]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
  
  if (verbose) { message("Downloading the genomic data of the Cancer Genome Project from InSilicoDB") }
  InSilicoLogin(login="bhaibeka@gmail.com", password="747779bec8a754b91076d6cc1f700831")
  # inSilicoDb2::getCurationInfo(dataset="ISDB12210")
  platf <- inSilicoDb2::getPlatforms(dataset="ISDB12210")
  eset <- inSilicoDb2::getDatasets(dataset="ISDB12210", norm="FRMA", curation="24802", features="PROBE")
  InSilicoLogout()
  
  ## only one platform, may be subject to change
  platf <- platf[[1]]
  eset <- eset[[1]]
  
  colnames(Biobase::pData(eset)) <- gsub(badchars, "_", colnames(Biobase::pData(eset)))
  
  ## gene centric expression
  if (gene) {
    if (verbose) { message("Gene centric data") }
    eset <- MetaGx::probeGeneMapping(eset=eset, platform="GPL96", method="jetset")
  }
  
  ## replicated experiments
  pheno <- Biobase::pData(eset)
  switch(replicates,
    "first" = {
      if (verbose) { message("First experiment of each replicate is kept") }
      iix <- order(pheno[ , "file_day"], pheno[ , "file_hour"], decreasing=FALSE, na.last=TRUE)
      ix <- rownames(pheno)[iix][!duplicated(pheno[iix, "cell_id"])]
      Biobase::exprs(eset) <- Biobase::exprs(eset)[ , ix, drop=FALSE]
      Biobase::pData(eset) <- Biobase::pData(eset)[ix, , drop=FALSE]
    },
    "last" = {
      if (verbose) { message("Last experiment of each replicate is kept") }
      iix <- order(pheno[ , "file_day"], pheno[ , "file_hour"], decreasing=TRUE, na.last=TRUE)
      ix <- rownames(pheno)[iix][!duplicated(pheno[iix, "cell_id"])]
      Biobase::exprs(eset) <- Biobase::exprs(eset)[ , ix, drop=FALSE]
      Biobase::pData(eset) <- Biobase::pData(eset)[ix, , drop=FALSE]
    },
    "mean" = {
      if (verbose) { message("Mean of replicates is computed") }
      stop(sprintf("Method replicates %s not implemented yet", replicates))
    },
    "median" = {
      if (verbose) { message("Median of replicates is computed") }
      stop(sprintf("Method replicates %s not implemented yet", replicates))
      ix <- pheno[duplicated(pheno[ , "cell_id"]), "cell_id"]
      nn <- NULL
      ex <- matrix(NA, nrow=nrow(Biobase::fData(eset)), ncol=sum(!duplicated(pheno[ , "cell_id"]), na.rm=TRUE), dimnames=list(rownames(Biobase::fData(eset)), nn))
      for (i in 1:length(ix)) {
        iix <- !is.na(pheno[ , "cell_id"]) & pheno[ , "cell_id"] == ix[i]
        xx <- apply(Biobase::exprs(eset)[ , iix], 1, get(replicates), na.rm=TRUE)
      }
    },
    "all" = {
      if (verbose) { message("All replicates are kept in the expression set") }
    }  
  )
  
  ## replace sample names by cell line names
  nn <- genefu::rename.duplicate(Biobase::pData(eset)[ , "cell_id"], sep="_")$new.x
  rownames(Biobase::pData(eset)) <- colnames(Biobase::exprs(eset)) <- nn
  
  ## get drug sensivity data
  if (verbose) { message("Download and format drug sensitivity data") }
  tmpfiles <- NULL

  ## download drug sensitivity (release 2)
  myfn <- file.path(tmpdir, "cgp_drug_sensitivity.csv")
  if(!file.exists(myfn)) {
    if (verbose) { message("Download drug sensitivity measurements") }
    dwl.status <- download.file(url="ftp://ftp.sanger.ac.uk/pub4/cancerrxgene/releases/release-2.0/gdsc_manova_input_w2.csv", destfile=myfn, method="curl")
    if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
    tmpfiles <- c(tmpfiles, myfn)
  }

  ## download drug concentration (release 2)
  myfn <- file.path(tmpdir, "cgp_drug_concentration.csv")
  if(!file.exists(myfn)) {
    if (verbose) { message("Download screening drug concentrations") }
    dwl.status <- download.file(url="ftp://ftp.sanger.ac.uk/pub4/cancerrxgene/current_release/gdsc_compounds_conc_w2.csv", destfile=myfn, method="curl")
    if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
    tmpfiles <- c(tmpfiles, myfn)
  }

  ## download cell line annotations and COSMIC IDs
  myfn <- file.path(tmpdir, "celline_annotations.RData")
  if(!file.exists(myfn)) {
    if (verbose) { message("Download cell lines annotations") }
    ## annotations from GDSC (Genomics of Drug Sensitivity in Cancer)
    dwl.status <- download.file(url="ftp://ftp.sanger.ac.uk/pub4/cancerrxgene/current_release/gdsc_cell_lines_w2.csv", destfile=file.path(tmpdir, "cgp_celline_collection.csv"), method="curl")
    tmpfiles <- c(tmpfiles, file.path(tmpdir, "cgp_celline_collection.csv"))
    if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
    celline.gdsc <- read.csv(file=file.path(tmpdir, "cgp_celline_collection.csv"), stringsAsFactors=FALSE)
    celline.gdsc[celline.gdsc == "" | celline.gdsc == " " | celline.gdsc == "  "] <- NA
    celline.gdsc <- celline.gdsc[!is.na(celline.gdsc[ , "CELL_LINE_NAME"]), , drop=FALSE]
    dupln <- unique(celline.gdsc[ , "CELL_LINE_NAME"][duplicated(celline.gdsc[ , "CELL_LINE_NAME"])])
    celline.gdsc <- celline.gdsc[!duplicated(celline.gdsc[ , "CELL_LINE_NAME"]), , drop=FALSE]
    celline.gdsc[ , "COSMIC_ID"] <- as.character(celline.gdsc[ , "COSMIC_ID"])
    rownames(celline.gdsc) <- celline.gdsc[ , "CELL_LINE_NAME"]
    ## annotations from COSMIC
    if (cosmic.annotation) {
      dwl.status <- download.file(url=sprintf("ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/data_export/CosmicCellLineProject_%s.tsv.gz", cosmic.version), destfile=file.path(tmpdir, sprintf("CosmicCellLineProject_%s.tsv.gz", cosmic.version)), method="curl")
      if(dwl.status != 0) { stop("Download failed, please rerun the pipeline! It may be that there is a new version of the file CosmicCellLineProject, please look at ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/data_export/ and update the script accordingly ...") }
      ## untar
      res <- R.utils::gunzip(filename=file.path(tmpdir, sprintf("CosmicCellLineProject_%s.tsv.gz", cosmic.version)), overwrite=TRUE)
      tmpfiles <- c(tmpfiles, file.path(tmpdir, sprintf("CosmicCellLineProject_%s.tsv", cosmic.version)))
      celline.cosmic <- read.csv(file=file.path(tmpdir, sprintf("CosmicCellLineProject_%s.tsv", cosmic.version)), sep="\t", stringsAsFactors=FALSE)
      celline.cosmic[celline.cosmic == "" | celline.cosmic == " " | celline.cosmic == "  "] <- NA
      ## remove cell line with no name
      celline.cosmic <- celline.cosmic[!is.na(celline.cosmic[ , "Sample.name"]), , drop=FALSE]
      ## merge the gene targets
      dupln <- unique(celline.cosmic[ , "Sample.name"][duplicated(celline.cosmic[ , "Sample.name"])])
      tt <- celline.cosmic
      iix.rm <- NULL
      for(i in 1:length(dupln)) {
        duplix <- celline.cosmic[ ,"Sample.name"] == dupln[i]
        iix <- sort((which(duplix)), decreasing=FALSE)[1]
        iix.rm <- c(iix.rm, setdiff(which(duplix), iix))
        tt[iix, "Gene.name"] <- paste(celline.cosmic[duplix, "Gene.name"], collapse="///")
        tt[iix, "UniProt.ID"] <- paste(celline.cosmic[duplix, "UniProt.ID"], collapse="///")
        tt[iix, "Zygosity"] <- paste(celline.cosmic[duplix, "Zygosity"], collapse="///")
        tt[iix, "CDS_MUT_SYNTAX"] <- paste(celline.cosmic[duplix, "CDS_MUT_SYNTAX"], collapse="///")
        tt[iix, "AA_MUT_SYNTAX"] <- paste(celline.cosmic[duplix, "AA_MUT_SYNTAX"], collapse="///")
        tt[iix, "NCBI36.genome.position"] <- paste(celline.cosmic[duplix, "NCBI36.genome.position"], collapse="///")
        tt[iix, "GRCh37.genome.position"] <- paste(celline.cosmic[duplix, "GRCh37.genome.position"], collapse="///")
      }
      tt <- tt[-iix.rm, , drop=FALSE]
      rownames(tt) <- tt[ , "Sample.name"]
      celline.cosmic <- tt
    } else {
      celline.cosmic <- cbind("Sample.name"=celline.gdsc[ , "CELL_LINE_NAME"], "ID_sample"=celline.gdsc[ , "COSMIC_ID"])
      rownames(celline.cosmic) <- celline.cosmic[ , "Sample.name"]
    }
    ## merge GDSC and COSMIC annotations through COSMIC_ID
    iix <- which(!is.na(celline.gdsc[ , "COSMIC_ID"]) & !is.element(celline.gdsc[ , "COSMIC_ID"], celline.cosmic[ , "ID_sample"]))
    if (length(iix) == 0) {
      tt <- data.frame(matrix(NA, nrow=nrow(celline.cosmic) + length(iix), ncol=ncol(celline.cosmic), dimnames=list(c(rownames(celline.cosmic), rownames(celline.gdsc)[iix]), colnames(celline.cosmic))))
      tt[rownames(celline.cosmic), ] <- celline.cosmic
      tt[rownames(celline.gdsc)[iix], "Sample.name"] <- celline.gdsc[iix, "CELL_LINE_NAME"]
      tt[rownames(celline.gdsc)[iix], "ID_sample"] <- celline.gdsc[iix, "COSMIC_ID"]
      celline <- tt
      colnames(celline)[match(c("Sample.name", "ID_sample"), colnames(celline))] <- c("CELL_LINE_NAME", "COSMIC_ID")
    } else {
      if (cosmic.annotation) {
        celline <- celline.cosmic
        colnames(celline)[match(c("Sample.name", "ID_sample"), colnames(celline))] <- c("CELL_LINE_NAME", "COSMIC_ID")
      } else {
        celline <- celline.gdsc
      }
    }
    save(list=c("celline.cosmic", "celline.gdsc", "celline"), compress=TRUE, file=myfn)
  } else {
    load(myfn)
  }

  ## download drug information
  myfn <- file.path(tmpdir, "nature11005-s2.zip")
  if(!file.exists(myfn)) {
    if (verbose) { message("Download drug information") }
    dwl.status <- download.file(url="http://www.nature.com/nature/journal/v483/n7391/extref/nature11005-s2.zip", destfile=myfn, method="curl")
    if(dwl.status != 0) { stop("Download failed, please rerun the script!") }
  }
  ff <- as.character(unzip(zipfile=file.path(tmpdir, "nature11005-s2.zip"), list=TRUE)[1, 1])
  unzip(zipfile=file.path(tmpdir, "nature11005-s2.zip"), exdir=file.path(tmpdir))
  tmpfiles <- c(tmpfiles, file.path(tmpdir, "Supplementary_data_final_Apr6.xls"))
  
  ## phenotype for the drugs
  if (verbose) { message("Read drug sensitivity measurements") }
  myfn <- file.path(tmpdir, "cgp_drug_sensitivity.RData")
  if(!file.exists(myfn)) {
    drugpheno <- read.csv(file.path(tmpdir, "cgp_drug_sensitivity.csv"), stringsAsFactors=FALSE)
    drugpheno[drugpheno == "" | drugpheno == " "] <- NA
    save(list="drugpheno", compress=TRUE, file=myfn)
    tmpfiles <- c(tmpfiles, myfn)
  } else { load(myfn) }
  ## format column names
  coln2 <- unlist(drugpheno[1, ,drop=TRUE])
  coln2[coln2 == ""] <- NA
  drugpheno <- drugpheno[-1, ,drop=FALSE]
  coln <- colnames(drugpheno)
  coln2[is.na(coln2)] <- coln[is.na(coln2)]
  coln2 <- genefu::rename.duplicate(x=coln2, sep="_dupl")$new.x
  myx <- sapply(sapply(strsplit(coln2, "_"), function(x) { return(x[[1]]) }), Hmisc::all.is.numeric)
  coln2[myx] <- paste("drugid", gsub(pattern=badchars, replacement="_", x=toupper(coln2[myx])), sep="_")
  colnames(drugpheno) <- coln2
  ## drug identifiers and names
  dn <- toupper(gsub(badchars, "", sapply(strsplit(coln, "_"), function(x) { return(x[[1]]) })))
  ## manual curation for drug names starting with a figure
  dn[!is.na(dn) & dn == "X17AAG"] <- "17AAG"
  dn[!is.na(dn) & dn == "X681640"] <- "681640"
  did <- sapply(strsplit(coln2, "_"), function(x) { if(x[[1]] == "drugid") { return(x[[2]]) } else { return(NA) } })
  drugnid <- cbind("drug.name"=dn, "drug_id"=did)[!is.na(did) & !duplicated(did), ]
  rownames(drugnid) <- paste("drugid", drugnid[ , "drug_id"], sep="_")

  ## cell line identifiers
  dupln <- duplicated(drugpheno[ ,"Cell.Line"])
  if(sum(dupln) > 1) { warning("some cell lines are duplicated, only the first instance is kept") }
  drugpheno <- drugpheno[!dupln, , drop=FALSE]
  if(any(!is.element(drugpheno[ ,"Cell.Line"], celline[ , "CELL_LINE_NAME"]))) { warning("Some cell line with drug sensitivity data have no annotations") }
  celln <- as.character(drugpheno[ ,"Cell.Line"])
  drugpheno <- data.frame("cell_id"=celln, drugpheno, stringsAsFactors=FALSE)
  rownames(drugpheno) <- celln

  ## get mutational data, i.e., protein coding variants
  ## Genetic mutation data for cancer genes. Includes MSI status (1 = unstable and 0 = stable) and gene-fusions. A binary code 'x::y' description is used for each gene where 'x' identifies a coding variant and 'y' indicates copy number information from SNP6.0 data. For gene fusions, cell lines are identified as fusion not-detected (0) or the identified fusion is given. The following abbreviations are used: not analysed (na), not detected or wild-type (wt), no copy number information (nci).
  ## we assume that AKT2 and WT1 are the first and last genes in the file
  rangeg <- which(colnames(drugpheno) == "AKT2"):which(colnames(drugpheno) == "WT1")
  mutation <- as.matrix(drugpheno[ , rangeg, drop=FALSE])
  mutation <- apply(X=mutation, MARGIN=c(1, 2), FUN=function(x) {
    x <- unlist(strsplit(x, split="::"))
    if(length(x) == 2) {
      if(!is.na(x[[1]]) && (x[[1]] == "na" || x[[1]] == "p.?" || x[[1]] == "p.0?")) {
        x <- NA
      } else {
        x <- x[[1]]
      }
    } else { x <- NA }
    return(x)
  })
  
  ## url for cell line collection
  celline <- data.frame("cell_id"=as.character(celline[ , "CELL_LINE_NAME"]), celline, stringsAsFactors=FALSE)
  ## add url based on COSMIC IDs
  uurl <- paste("http://cancer.sanger.ac.uk/cell_lines/sample/overview?id=", celline[ , "COSMIC_ID"], sep="")
  uurl[is.na(celline[ , "COSMIC_ID"])] <- NA
  celline <- data.frame("cell_id"=celline[ , "cell_id"], "link"=uurl, celline[ , !is.element(colnames(celline), "cell_id")], stringsAsFactors=FALSE)

  ## make sure that cell_id are not factors
  celline[, "cell_id"] <- as.character(celline[, "cell_id"])
  Biobase::pData(eset)[, "cell_id"] <- as.character(Biobase::pData(eset)[, "cell_id"])
  drugpheno[, "cell_id"] <- as.character(drugpheno[, "cell_id"])
  
  ## union of all cell line with data
  cellnall <- sort(unique(c(row.names(Biobase::pData(eset)), rownames(mutation), as.character(drugpheno[ ,"cell_id"]))))
  
  ## update sampleinfo
  dd <- data.frame(matrix(NA, ncol=ncol(Biobase::pData(eset)), nrow=length(cellnall), dimnames=list(cellnall, colnames(Biobase::pData(eset)))), check.names=FALSE)
  dd[rownames(Biobase::pData(eset)), colnames(Biobase::pData(eset))] <- Biobase::pData(eset)
  Biobase::pData(eset) <- dd
  
  ## update gene expressions
  dd <- matrix(NA, nrow=nrow(Biobase::exprs(eset)), ncol=length(cellnall), dimnames=list(rownames(Biobase::exprs(eset)), cellnall))
  dd[rownames(Biobase::exprs(eset)), colnames(Biobase::exprs(eset))] <- Biobase::exprs(eset)
  Biobase::exprs(eset) <- dd
  
  ## update drugpheno
  dd <- data.frame(matrix(NA, ncol=ncol(drugpheno), nrow=length(cellnall)))
  rownames(dd) <- cellnall
  colnames(dd) <- colnames(drugpheno)
  newlev <- sapply(drugpheno, levels)
  newlev$cell_id <- cellnall
  dd <- genefu::setcolclass.df(df=dd, colclass=sapply(drugpheno, class), factor.levels=newlev)
  dd[rownames(drugpheno),colnames(drugpheno)] <- drugpheno
  dd[ ,"cell_id"] <- cellnall
  drugpheno <- dd

  ## update mutation
  dd <- matrix(NA, ncol=ncol(mutation), nrow=length(cellnall), dimnames=list(cellnall, colnames(mutation)))
  dd[rownames(mutation), colnames(mutation)] <- mutation
  mutation <- t(dd)
  
  ## update celline
  dd <- data.frame(matrix(NA, ncol=ncol(celline), nrow=length(cellnall), dimnames=list(cellnall, colnames(celline))), check.names=FALSE)
  iix <- intersect(rownames(celline), cellnall)
  dd[iix, colnames(celline)] <- celline[iix, , drop=FALSE]
  celline <- dd
  celline[ , "cell_id"] <- celline[ , "CELL_LINE_NAME"] <- rownames(celline)
  ## annotate cell lines with curated tissue type
  tissue.type <- read.csv(file.path(system.file("extdata", package="PharmacoGx"), "cell_line_collection_all.csv"), stringsAsFactors=FALSE)
  rownames(tissue.type) <- tissue.type[ , 1]
  celline <- cbind("tissue.type"=tissue.type[match(celline[ , "cell_id"], tissue.type[ , "cell_id"]), "tissue.type"], celline)

  ## drug information
  if (verbose) { message("Read drug information") }
  myfn2 <- file.path(tmpdir, "nature_supplinfo_druginfo_cgp.RData")
  if(!file.exists(myfn2)) {
    druginfo <- gdata::read.xls(xls=file.path(tmpdir, "Supplementary_data_final_Apr6.xlsx"), sheet=4)
    druginfo[druginfo == "" | druginfo == " "] <- NA
    save(list="druginfo", compress=TRUE, file=myfn2)
  } else { load(myfn2) }
  druginfo <- data.frame("drug_id"=gsub(pattern =badchars, replacement="", x=toupper(druginfo[ ,"Drug.ID"])), druginfo, stringsAsFactors=FALSE)
  rownames(druginfo) <- druginfo[ ,"drug_id"] <- paste("drugid", as.character(druginfo[ ,"drug_id"]), sep="_")

  ## drug concentration
  if (verbose) { message("Read drug concentration") }
  drugconc <- read.csv(file.path(tmpdir, "cgp_drug_concentration.csv"), stringsAsFactors=FALSE)
  drugconc[drugconc == "" | drugconc == " "] <- NA
  drugconc <- data.frame("drug.name"=toupper(gsub(badchars, "", drugconc[ ,"Compound.Name"])), drugconc, stringsAsFactors=FALSE)
  if(all(!is.element(drugconc[ , "drug.name"], drugnid[ , "drug.name"]))) { stop("Screening concentration for drugs ithout identifiers!") }
  rownames(drugconc) <- rownames(drugnid)[match(drugconc[ , "drug.name"], drugnid[ , "drug.name"])]
  drugconc <- data.frame("drug_id"=rownames(drugconc), drugconc, stringsAsFactors=FALSE)

  ## combine all drugs
  dix <- sort(unique(c(rownames(druginfo), rownames(drugconc), paste("drugid", sapply(strsplit(colnames(drugpheno)[grep("^drugid_", colnames(drugpheno))], "_"), function(x) { return(x[[2]]) }), sep="_"))))
  ## update druginfo
  druginfo2 <- data.frame(matrix(NA, nrow=length(dix), ncol=ncol(druginfo), dimnames=list(dix, colnames(druginfo))))
  newlev <- sapply(druginfo, levels)
  newlev$drug_id <- sapply(strsplit(dix, split="_"), function(x) { return(x[2]) })
  druginfo2 <- genefu::setcolclass.df(df=druginfo2, colclass=sapply(druginfo, class), factor.levels=newlev)
  druginfo2[match(rownames(druginfo), dix), colnames(druginfo)] <- druginfo
  druginfo2[ , "drug_id"] <- newlev$drug_id
  druginfo <- druginfo2
  ## update drugconc
  drugconc2 <- data.frame(matrix(NA, nrow=length(dix), ncol=ncol(drugconc), dimnames=list(dix, colnames(drugconc))))
  newlev <- sapply(drugconc, levels)
  newlev$drug_id <- sapply(strsplit(dix, split="_"), function(x) { return(x[2]) })
  drugconc2 <- genefu::setcolclass.df(df=drugconc2, colclass=sapply(drugconc, class), factor.levels=newlev)
  drugconc2[match(rownames(drugconc), dix), colnames(drugconc)] <- drugconc
  drugconc2[ , "drug_id"] <- newlev$drug_id
  drugconc <- drugconc2

  ## report concentrations per cell line and per drug
  drugconc2 <- data.frame(matrix(NA, nrow=nrow(drugconc) * length(cellnall), ncol=6, dimnames=list(paste(rep(cellnall, times=nrow(drugconc)), rep(rownames(drugconc), each=length(cellnall)), sep="..."), c("cell_id", "drug_id", "drug_name", "nbr_conc_tested", "min_Dose_uM", "max_Dose_uM"))))
  drugconc2[ , "cell_id"] <- rep(cellnall, times=nrow(drugconc))
  drugconc2[ , "drug_id"] <- rep(rownames(drugconc), each=length(cellnall))
  drugconc2[ , "drug_name"] <- rep(as.character(drugconc[ ,"drug.name"]), each=length(cellnall))
  ## as mentioned in the supplementary information of Garnett et al., a single cell line is used on each plate and treated with 28 different drugs over a 9-pt, 256-fold concentration range
  drugconc2[ , "nbr_conc_tested"] <- 9
  drugconc2[ , "min_Dose_uM"] <- rep(drugconc[ , "Min.Concentration.micromolar."], each=length(cellnall))
  drugconc2[ , "max_Dose_uM"] <- rep(drugconc[ , "Max.Concentration.micromolar."], each=length(cellnall))
  drugconc <- drugconc2
  drugconc[, "cell_id"] <- as.character(drugconc[, "cell_id"])

  ## IC50 in micro molar
  if (verbose) { message("Extracting IC50 values") }

  myx <- sapply(strsplit(colnames(drugpheno), "_"), function(x) { return(all(x[c(length(x)-1, length(x))] == c("IC", "50"))) })
  ic50 <- drugpheno[ ,myx,drop=FALSE]
  nn <- dimnames(ic50)
  nn[[2]] <- gsub("_IC_50", "", nn[[2]])
  ic50 <- apply(ic50, 2, as.numeric)
  dimnames(ic50) <- nn
  ic50 <- exp(ic50)

  ## sensitivity calling using waterfall plot
  ic50.call <- NULL
  for(i in 1:ncol(ic50)) {
    ic50.call <- cbind(ic50.call, callingWaterfall(x=ic50[ ,i], type="ic50", intermediate.fold=c(4, 1.2, 1.2), cor.min.linear=0.95, plot=FALSE, name=sprintf("%s (CGP)", colnames(ic50)[i])))
  }
  dimnames(ic50.call) <- dimnames(ic50)

  ## activity area
  if (verbose) { message("Extracting AUC values") }

  ## continuous values
  myx <- sapply(strsplit(colnames(drugpheno), "_"), function(x) { return(all(x[c(length(x))] == c("AUC"))) })
  auc <- drugpheno[ , myx, drop=FALSE]
  nn <- dimnames(auc)
  nn[[2]] <- gsub("_AUC", "", nn[[2]])
  auc <- apply(auc, 2, as.numeric)
  ## AUC for sensitivity
  auc <- 1 - auc
  dimnames(auc) <- nn
  
  ## sensitivity calling using waterfall plot
  auc.call <- NULL
  for(i in 1:ncol(auc)) {
    auc.call <- cbind(auc.call, callingWaterfall(x=auc[ ,i], type="actarea", intermediate.fold=c(4, 1.2, 1.2), cor.min.linear=0.95, plot=FALSE, name=sprintf("%s (CGP)", colnames(auc)[i])))
  }
  dimnames(auc.call) <- dimnames(auc)
   
  if (delete.tmpdir){
    ## delete temporary files
    sapply(tmpfiles, function (x) { unlink(x, recursive=TRUE) })
  }

  return (list("expression"=eset, "mutation"=mutation, "cell.line"=celline, "drug"=list("info"=druginfo, "concentration"=drugconc, "duration"=NULL, "ic50"=ic50, "ic50.call"=ic50.call, "auc"=auc, "auc.call"=auc.call)))
}

## End
