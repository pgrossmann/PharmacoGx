########################
## Benjamin Haibe-Kains
## All rights Reserved
## October 23, 2013
########################

#	Normalize and import all data from Connectivity Map

`normalizationCGP` <- 
function (datadir=file.path("data", "CGP"), tmpdir="tmp", cosmic.version="v67_241013", nthread=1) {
  
  ## number of cores for parallel processing
  nbcore <- nthread
  availcore <- parallel::detectCores()
  if (missing(nbcore) || nbcore > availcore) { nbcore <- availcore }
  options("mc.cores"=nbcore)
  
  ## create directories for temporary files and raw data
  if(!file.exists(tmpdir)) { dir.create(tmpdir, showWarnings=FALSE, recursive=TRUE) }
  rawdir <- file.path(datadir, "raw_ge")
  if(!file.exists(rawdir)) { dir.create(rawdir, showWarnings=FALSE, recursive=TRUE) }

  ########################
  ## download data
  ########################
  ftpdir <- "ftp://ftp.ebi.ac.uk//pub/databases/microarray/data/experiment/MTAB/E-MTAB-783/"
  myfn <- file.path(rawdir, "celfile_timestamp.RData")
  if(!file.exists(myfn)) {
    message("\nDownload genomic data\n")
  
    ## download and compress CEL files
    celfile.timestamp <- celfn <- NULL
    i <- 1
    while(i <= 9) {
      ## assuming there are only 9 zip archives (need to check if the update version has more)
     dwl.status <- download.file(url=sprintf("%s/E-MTAB-783.raw.%i.zip", ftpdir, i), destfile=file.path(tmpdir, sprintf("E-MTAB-783.raw.%i.zip", i)))
     if(dwl.status != 0) {
       message("\t-> download failed, let's try again ...")
       file.remove(file.path(tmpdir, sprintf("E-MTAB-783.raw.%i.zip", i)))
       i <- i - 1
      } else {
         ## unzip archive
         fff <- unzip(zipfile=file.path(tmpdir, sprintf("E-MTAB-783.raw.%i.zip", i)), list=TRUE)
         celfile.timestamp <- c(celfile.timestamp, as.character(fff[ ,"Date"]))
         celfn <- c(celfn, as.character(fff[ ,"Name"]))
         res <- unzip(zipfile=file.path(tmpdir, sprintf("E-MTAB-783.raw.%i.zip", i)), exdir=rawdir)
         ## compress each CEL file individually using gzip
         sapply(file.path(rawdir, as.character(fff[ ,"Name"])), R.utils::gzip, overwrite=TRUE)
         i <- i + 1
       }
    }
    celfile.timestamp <- t(sapply(strsplit(celfile.timestamp, split=" "), function(x) { return(x) }))
    dimnames(celfile.timestamp) <- list(celfn, c("file.day", "file.hour"))
   
    # unlink(file.path(tmpdir), recursive=TRUE)
    write.csv(celfile.timestamp, file=file.path(rawdir, "celfile_timestamp.csv"))
    save(list=c("celfile.timestamp"), compress=TRUE, file=myfn)
  }

  ## download sample information
  myfn <- file.path(tmpdir, "E-MTAB-783.sdrf.txt")
  if(!file.exists(myfn)) {
    message("\nDownload sample information\n")
    dwl.status <- download.file(url=sprintf("%s/E-MTAB-783.sdrf.txt", ftpdir), destfile=myfn)
    if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
    file.copy(from=file.path(tmpdir, "E-MTAB-783.sdrf.txt"), to=file.path(rawdir, "E-MTAB-783.sdrf.txt"))
  }
  
  ## download drug sensitivity (release 2)
  myfn <- file.path(tmpdir, "gdsc_manova_input_w2.csv")
  if(!file.exists(myfn)) {
    message("\nDownload drug sensitivity measurements\n")
    dwl.status <- download.file(url="ftp://ftp.sanger.ac.uk/pub4/cancerrxgene/releases/release-2.0/gdsc_manova_input_w2.csv", destfile=myfn)
    if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
    file.copy(from=file.path(tmpdir, "gdsc_manova_input_w2.csv"), to=file.path(rawdir, "cgp_drug_sensitivity.csv"))
  }

  ## download drug concentration (release 2)
  myfn <- file.path(tmpdir, "gdsc_compounds_conc_w2.csv")
  if(!file.exists(myfn)) {
    message("\nDownload screening drug concentrations\n")
    dwl.status <- download.file(url="ftp://ftp.sanger.ac.uk/pub4/cancerrxgene/current_release/gdsc_compounds_conc_w2.csv", destfile=myfn)
    if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
    file.copy(from=file.path(tmpdir, "gdsc_compounds_conc_w2.csv"), to=file.path(rawdir, "cgp_drug_concentration.csv"))
  }
  
  ## download cell line annotations and COSMIC IDs
  ## annotations from COSMIC cell line project
  dda <- cosmic.version
  myfn <- file.path(tmpdir, "celline_annotations.RData")
  if(!file.exists(myfn)) {
    message("\nDownload cell lines annotations\n")
    dwl.status <- download.file(url=sprintf("ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/data_export/CosmicCellLineProject_%s.tsv.gz", dda), destfile=file.path(tmpdir, sprintf("CosmicCellLineProject_%s.tsv.gz", dda)))
    if(dwl.status != 0) { stop("Download failed, please rerun the pipeline! It may be that there is a new version of the file CosmicCellLineProject, please look at ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/data_export/ and update the script accordingly ...") }
    ## untar
    res <- R.utils::gunzip(filename=file.path(tmpdir, sprintf("CosmicCellLineProject_%s.tsv.gz", dda)), overwrite=TRUE)
    file.copy(from=file.path(tmpdir, sprintf("CosmicCellLineProject_%s.tsv", dda)), to=file.path(rawdir, "cosmic_celline_collection.csv"))
    celline.cosmic <- read.csv(file=file.path(rawdir, "cosmic_celline_collection.csv"), sep="\t", stringsAsFactors=FALSE)
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
    ## annotations from GDSC (Genomics of Drug Sensitivity in Cancer)
    dwl.status <- download.file(url="ftp://ftp.sanger.ac.uk/pub4/cancerrxgene/current_release/gdsc_cell_lines_w2.csv", destfile=file.path(tmpdir, "gdsc_cell_lines_w2.csv"))
    file.copy(from=file.path(tmpdir, "gdsc_cell_lines_w2.csv"), to=file.path(rawdir, "cgp_celline_collection.csv"))
    if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
    celline.gdsc <- read.csv(file=file.path(rawdir, "cgp_celline_collection.csv"), stringsAsFactors=FALSE)
    celline.gdsc[celline.gdsc == "" | celline.gdsc == " " | celline.gdsc == "  "] <- NA
    celline.gdsc <- celline.gdsc[!is.na(celline.gdsc[ , "CELL_LINE_NAME"]), , drop=FALSE]
    dupln <- unique(celline.gdsc[ , "CELL_LINE_NAME"][duplicated(celline.gdsc[ , "CELL_LINE_NAME"])])
    celline.gdsc <- celline.gdsc[!duplicated(celline.gdsc[ , "CELL_LINE_NAME"]), , drop=FALSE]
    rownames(celline.gdsc) <- celline.gdsc[ , "CELL_LINE_NAME"]
    ## merge GDSC and COSMIC annotations through COSMIC_ID
    iix <- which(!is.na(celline.gdsc[ , "COSMIC_ID"]) & !is.element(celline.gdsc[ , "COSMIC_ID"], celline.cosmic[ , "ID_sample"]))
    tt <- data.frame(matrix(NA, nrow=nrow(celline.cosmic) + length(iix), ncol=ncol(celline.cosmic), dimnames=list(c(rownames(celline.cosmic), rownames(celline.gdsc)[iix]), colnames(celline.cosmic))))
    tt[rownames(celline.cosmic), ] <- celline.cosmic
    tt[rownames(celline.gdsc)[iix], "Sample.name"] <- celline.gdsc[iix, "CELL_LINE_NAME"]
    tt[rownames(celline.gdsc)[iix], "ID_sample"] <- celline.gdsc[iix, "COSMIC_ID"]
    celline.cgp <- tt
    save(list=c("celline.cosmic", "celline.gdsc", "celline.cgp"), compress=TRUE, file=myfn)
  }

  ## download drug information
  myfn <- file.path(tmpdir, "nature11005-s2.zip")
  if(!file.exists(myfn)) {
    message("\nDownload drug information\n")
    dwl.status <- download.file(url="http://www.nature.com/nature/journal/v483/n7391/extref/nature11005-s2.zip", destfile=myfn)
    if(dwl.status != 0) { stop("Download failed, please rerun the script!") }
  }
  ff <- as.character(unzip(zipfile=file.path(tmpdir, "nature11005-s2.zip"), list=TRUE)[1, 1])
  unzip(zipfile=file.path(tmpdir, "nature11005-s2.zip"), exdir=file.path(tmpdir))
  file.copy(from=file.path(tmpdir, ff), to=file.path(rawdir, "nature_supplementary_information.xls"))

  ########################
  ## normalize and format data
  ########################
  myfn <- file.path(tmpdir, "cgp_frma.RData")
  if(!file.exists(myfn)) {

    data(hthgu133afrmavecs)
    data(hthgu133acdf)

    ## CEL file names
    celfn <- list.celfiles(rawdir, full.names=TRUE)
    celfns <- list.celfiles(rawdir, full.names=FALSE)
    ## experiments' names
    names(celfn) <- names(celfns) <- gsub(".CEL.gz", "", celfns)
    ## chip type and date
    chipt <- sapply(celfn, celfileChip)
    chipd <- t(sapply(celfn, celfileDateHour))
    ## reorder CEL files by hybridization time or timestamp
    myx <- NULL
    if(any(!complete.cases(chipd))) {
      ## all hybridization dates are not available
      load(file.path(rawdir, "celfile_timestamp.RData"))
      if(!all(is.element(celfns, paste(rownames(celfile.timestamp), "gz", sep=".")))) { stop("Timestamp is not available for all CEL files!") }
        celfile.timestamp <- celfile.timestamp[match(celfns, paste(rownames(celfile.timestamp), "gz", sep=".")), , drop=FALSE] 
        myx <- order(celfile.timestamp[ ,"file.day"], celfile.timestamp[ ,"file.hour"], decreasing=FALSE)
    } else {
        myx <- order(chipd[ ,"day"], chipd[ ,"hour"], decreasing=FALSE)
    }
    celfn <- celfn[myx]
    celfns <- celfns[myx]
    chipt <- chipt[myx]
    chipd <- chipd[myx, , drop=FALSE]
    celfile.timestamp <- celfile.timestamp[myx, , drop=FALSE]

    ## read info about drugs and experiments

    ## phenotype for the drugs
    message("\nRead drug sensitivity measurements\n")
    myfn2 <- file.path(tmpdir, "cgp_drug_sensitivity.RData")
    if(!file.exists(myfn2)) {
      drugpheno <- read.csv(file.path(rawdir, "cgp_drug_sensitivity.csv"), stringsAsFactors=FALSE)
      drugpheno[drugpheno == "" | drugpheno == " "] <- NA
      save(list="drugpheno", compress=TRUE, file=myfn2)
    } else { load(myfn2) }
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
    drugnid <- cbind("drug.name"=dn, "drug.id"=did)[!is.na(did) & !duplicated(did), ]
    rownames(drugnid) <- paste("drugid", drugnid[ , "drug.id"], sep="_")

    ## cell line identifiers
    dupln <- duplicated(drugpheno[ ,"Cell.Line"])
    if(sum(dupln) > 1) { warning("some cell lines are duplicated, only the first instance is kept") }
    drugpheno <- drugpheno[!dupln, , drop=FALSE]
    if(any(!is.element(drugpheno[ ,"Cell.Line"], celline.cgp[ , "Sample.name"]))) { stop("Some cell line names are not included in the COSMIC database") }
    celln <- drugpheno[ ,"Cell.Line"]
    drugpheno <- data.frame("cellid"=celln, drugpheno)
    rownames(drugpheno) <- celln
  
    ## protein coding variants
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

    ## info about each experiment
    message("\nRead sample information\n")
    sampleinfo <- read.csv(file.path(rawdir, "E-MTAB-783.sdrf.txt"), sep="\t", stringsAsFactors=FALSE)
    sampleinfo[sampleinfo == "" | sampleinfo == " "] <- NA
    ## curate cell line names
    sampleinfo[sampleinfo[ , "Source.Name"] == "MZ2-MEL.", "Source.Name"] <- "MZ2-MEL"
    iix <- which(!duplicated(sampleinfo[ , "Source.Name"]) & !is.element(sampleinfo[ , "Source.Name"], celline.cgp[ , "Sample.name"]))
    if(length(iix) > 0) {
      ## enrich the list of cell lines
      tt <- matrix(NA, nrow=length(iix), ncol=ncol(celline.cgp), dimnames=list(sampleinfo[iix, "Source.Name"], colnames(celline.cgp)))
      tt[ , "Sample.name"] <- sampleinfo[iix, "Source.Name"]
      celline.cgp <- rbind(celline.cgp, tt)
    }
    fn <- gsub(patter="[.]CEL", replacement="", x=sampleinfo[ ,"Array.Data.File"])
    if(any(!is.element(fn[!is.na(fn)], names(celfns)))) { stop("some CEL files are missing for the CGP project") }
    rownames(sampleinfo) <- fn
    sampleinfo <- sampleinfo[names(celfn), , drop=FALSE]
    sampleinfo <- data.frame("samplename"=names(celfns), "filename"=celfns, "chiptype"=chipt, "hybridization.date"=chipd[ ,"day"], "hybridization.hour"=chipd[ ,"hour"], "file.day"=celfile.timestamp[ ,"file.day"], "file.hour"=celfile.timestamp[ ,"file.hour"], "batch"=NA, "cellid"=sampleinfo[ , "Source.Name"], sampleinfo)
    sampleinfo2 <- sampleinfo
    
    ## remove duplcated cell line hybridization
    sampleinfo <- sampleinfo[!duplicated(sampleinfo[ ,"cellid"]), , drop=FALSE]
    rownames(sampleinfo) <- sampleinfo[ ,"cellid"]

    ## update of cgp cell line collection
    celline.cgp <- data.frame("cellid"=as.character(celline.cgp[ , "Sample.name"]), celline.cgp)
    celline.cgp[ , "cellid"] <- as.character(celline.cgp[ , "cellid"])
    ## add url based on COSMIC IDs
    uurl <- paste("http://cancer.sanger.ac.uk/cell_lines/sample/overview?id=", celline.cgp[ , "ID_sample"], sep="")
    uurl[is.na(celline.cgp[ , "ID_sample"])] <- NA
    celline.cgp <- data.frame("cellid"=celline.cgp[ , "cellid"], "link"=uurl, celline.cgp[ , !is.element(colnames(celline.cgp), "cellid")])

    ## drugpheno
    cellnall <- sort(unique(c(as.character(sampleinfo[ ,"cellid"]), as.character(drugpheno[ ,"cellid"]))))
    dd <- data.frame(matrix(NA, ncol=ncol(drugpheno), nrow=length(cellnall), dimnames=list(cellnall, colnames(drugpheno))))
    newlev <- sapply(drugpheno, levels)
    newlev$cellid <- cellnall
    dd <- setcolclass.df(df=dd, colclass=sapply(drugpheno, class), factor.levels=newlev)
    dd[rownames(drugpheno),colnames(drugpheno)] <- drugpheno
    dd[ ,"cellid"] <- cellnall
    drugpheno <- dd
  
    ## mutation
    dd <- matrix(NA, ncol=ncol(mutation), nrow=length(cellnall), dimnames=list(cellnall, colnames(mutation)))
    dd[rownames(mutation), colnames(mutation)] <- mutation
    rownames(dd) <- cellnall
    mutation <- dd

    ## drug information
    message("\nRead drug information\n")
    myfn2 <- file.path(tmpdir, "nature_supplinfo_druginfo_cgp.RData")
    if(!file.exists(myfn2)) {
      druginfo <- gdata::read.xls(xls=file.path(rawdir, "nature_supplementary_information.xls"), sheet=4)
      druginfo[druginfo == "" | druginfo == " "] <- NA
      save(list="druginfo", compress=TRUE, file=myfn2)
    } else { load(myfn2) }
    druginfo <- data.frame("drugid"=gsub(pattern =badchars, replacement="", x=toupper(druginfo[ ,"Drug.ID"])), druginfo)
    rownames(druginfo) <- paste("drugid", as.character(druginfo[ ,"drugid"]), sep="_")

    ## drug concentration
    message("\nRead drug concentration\n")
    drugconc <- read.csv(file.path(rawdir, "cgp_drug_concentration.csv"), stringsAsFactors=FALSE)
    drugconc[drugconc == "" | drugconc == " "] <- NA
    drugconc <- data.frame("drug.name"=toupper(gsub(badchars, "", drugconc[ ,"Compound.Name"])), drugconc)
    if(all(!is.element(drugconc[ , "drug.name"], drugnid[ , "drug.name"]))) { stop("Screening concentration for drugs ithout identifiers!") }
    rownames(drugconc) <- rownames(drugnid)[match(drugconc[ , "drug.name"], drugnid[ , "drug.name"])]
    drugconc <- data.frame("drugid"=rownames(drugconc), drugconc)

    ## combine all drugs
    dix <- sort(unique(c(rownames(druginfo), rownames(drugconc), paste("drugid", sapply(strsplit(colnames(drugpheno)[grep("^drugid_", colnames(drugpheno))], "_"), function(x) { return(x[[2]]) }), sep="_"))))
    ## update druginfo
    druginfo2 <- data.frame(matrix(NA, nrow=length(dix), ncol=ncol(druginfo), dimnames=list(dix, colnames(druginfo))))
    newlev <- sapply(druginfo, levels)
    newlev$drugid <- sapply(strsplit(dix, split="_"), function(x) { return(x[2]) })
    druginfo2 <- setcolclass.df(df=druginfo2, colclass=sapply(druginfo, class), factor.levels=newlev)
    druginfo2[match(rownames(druginfo), dix), colnames(druginfo)] <- druginfo
    druginfo2[ , "drugid"] <- newlev$drugid
    druginfo <- druginfo2
    ## update drugconc
    drugconc2 <- data.frame(matrix(NA, nrow=length(dix), ncol=ncol(drugconc), dimnames=list(dix, colnames(drugconc))))
    newlev <- sapply(drugconc, levels)
    newlev$drugid <- sapply(strsplit(dix, split="_"), function(x) { return(x[2]) })
    drugconc2 <- setcolclass.df(df=drugconc2, colclass=sapply(drugconc, class), factor.levels=newlev)
    drugconc2[match(rownames(drugconc), dix), colnames(drugconc)] <- drugconc
    drugconc2[ , "drugid"] <- newlev$drugid
    drugconc <- drugconc2

    ## report concentrations per cell line and per drug
    drugconc2 <- data.frame(matrix(NA, nrow=nrow(drugconc) * length(cellnall), ncol=6, dimnames=list(paste(rep(cellnall, times=nrow(drugconc)), rep(rownames(drugconc), each=length(cellnall)), sep="..."), c("cellid", "drugid", "drug.name", "nbr.conc.tested", "min.Dose.uM", "max.Dose.uM"))))
    drugconc2[ , "cellid"] <- rep(cellnall, times=nrow(drugconc))
    drugconc2[ , "drugid"] <- rep(rownames(drugconc), each=length(cellnall))
    drugconc2[ , "drug.name"] <- rep(as.character(drugconc[ ,"drug.name"]), each=length(cellnall))
    ## as mentioned in the supplementary information of Garnett et al., a single cell line is used on each plate and treated with 28 different drugs over a 9-pt, 256-fold concentration range
    drugconc2[ , "nbr.conc.tested"] <- 9
    drugconc2[ , "min.Dose.uM"] <- rep(drugconc[ , "Min.Concentration.micromolar."], each=length(cellnall))
    drugconc2[ , "max.Dose.uM"] <- rep(drugconc[ , "Max.Concentration.micromolar."], each=length(cellnall))
    drugconc <- drugconc2

    ## normalization
    message("\nNormalize gene expression data\n")
    # rr <- just.rma(filenames=celfn)
    ## frma normalization using parallel
    splitix <- parallel::splitIndices(nx=length(celfn), ncl=nbcore)
    splitix <- splitix[sapply(splitix, length) > 0]
    res <- parallel::mclapply(splitix, function(x, celfn) {
      ## fRMA
      tt <- celfn[x]
      names(tt) <- NULL
      abatch <- affy::read.affybatch(filenames=tt)
      rr <- frma::frma(object=abatch, summarize="robust_weighted_average", verbose=TRUE, input.vecs=hthgu133afrmavecs)
      rr <- exprs(rr)
    }, celfn=celfn)
    datat <- t(do.call(cbind, res))

    ## build annotation matrix
    message("Build annotation matrix")
    mart.db <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org", dataset="hsapiens_gene_ensembl")
    ## select the best probe for a single gene
    js <- jetset::jscores(chip="hgu133a", probeset=colnames(datat))
    js <- js[colnames(datat), , drop=FALSE]
    ## identify the best probeset for each entrez gene id
    geneid1 <- as.character(js[ ,"EntrezID"])
    names(geneid1) <- rownames(js)
    geneid2 <- sort(unique(geneid1))
    names(geneid2) <- paste("geneid", geneid2, sep=".")
    gix1 <- !is.na(geneid1)
    gix2 <- !is.na(geneid2)
    geneid.common <- intersect(geneid1[gix1], geneid2[gix2])
    ## probes corresponding to common gene ids
    gg <- names(geneid1)[is.element(geneid1, geneid.common)]
    gid <- geneid1[is.element(geneid1, geneid.common)]
    ## duplicated gene ids
    gid.dupl <- unique(gid[duplicated(gid)])
    gg.dupl <- names(geneid1)[is.element(geneid1, gid.dupl)]
    ## unique gene ids
    gid.uniq <- gid[!is.element(gid, gid.dupl)]
    gg.uniq <- names(geneid1)[is.element(geneid1, gid.uniq)]
    ## which are the best probe for each gene
    js <- data.frame(js, "best"=FALSE)
    js[gg.uniq, "best"] <- TRUE
    ## data for duplicated gene ids
    if(length(gid.dupl) > 0) {	
    	## use jetset oevrall score to select the best probeset
    	myscore <- js[gg.dupl,"overall"]
    	myscore <- cbind("probe"=gg.dupl, "gid"=geneid1[gg.dupl], "score"=myscore)
    	myscore <- myscore[order(as.numeric(myscore[ , "score"]), decreasing=TRUE, na.last=TRUE), , drop=FALSE]
    	myscore <- myscore[!duplicated(myscore[ , "gid"]), , drop=FALSE]
    	js[myscore[ ,"probe"], "best"] <- TRUE
    }
    ## more annotations from biomart
    ugid <- sort(unique(js[ ,"EntrezID"]))
    ss <- "entrezgene"
    gene.an <- biomaRt::getBM(attributes=c(ss, "ensembl_gene_id", "hgnc_symbol", "unigene", "description", "chromosome_name", "start_position", "end_position", "strand", "band"), filters="entrezgene", values=ugid, mart=mart.db)
    gene.an[gene.an == "" | gene.an == " "] <- NA
    gene.an <- gene.an[!is.na(gene.an[ , ss]) & !duplicated(gene.an[ , ss]), , drop=FALSE]
    gene.an <- gene.an[is.element(gene.an[ ,ss], ugid), ,drop=FALSE]
    annot <- data.frame(matrix(NA, nrow=ncol(datat), ncol=ncol(gene.an)+1, dimnames=list(colnames(datat), c("probe", colnames(gene.an)))))
    annot[match(gene.an[ , ss], js[ ,"EntrezID"]), colnames(gene.an)] <- gene.an
    annot[ ,"probe"] <- colnames(datat)
    colnames(js)[colnames(js) != "best"] <- paste("jetset", colnames(js)[colnames(js) != "best"], sep=".")
    annot <- data.frame(annot, "EntrezGene.ID"=js[ ,"jetset.EntrezID"], js)
  
    ## save the full dataset, with duplicates
    sampleinfo.full.cgp <- sampleinfo2
    data.full.cgp <- datat
    rownames(data.full.cgp) <- gsub(".CEL.gz", "", rownames(data.full.cgp))
    data.full.cgp <- data.full.cgp[match(rownames(sampleinfo.full.cgp), rownames(data.full.cgp)), , drop=FALSE]
    annot.full.cgp <- annot
    save(list=c("data.full.cgp", "annot.full.cgp", "sampleinfo.full.cgp"), compress=TRUE, file=file.path(tmpdir, "cgp_full_frma.RData"))
    write.table(sampleinfo.full.cgp, row.names=FALSE, sep=";", file=file.path(tmpdir, "sampleinfo_full_cgp.csv"))
    
    ## match the experiment labels
    myx <- rownames(sampleinfo)[match(rownames(datat), as.character(sampleinfo[ ,"filename"]))]
    datat <- datat[!is.na(myx), , drop=FALSE]
    myx <- myx[!is.na(myx)]
    rownames(datat) <- myx

    ## keep only experiments for which we have all the info
    myx <- fold(intersect, rownames(drugpheno), rownames(sampleinfo), rownames(datat))
    data.ge.cgp <- datat[myx, ,drop=FALSE]
    annot.ge.cgp <- annot
    sampleinfo.ge.cgp <- sampleinfo[myx, , drop=FALSE]
    drugpheno.ge.cgp <- drugpheno[myx, , drop=FALSE]
    druginfo.ge.cgp <- druginfo
    drugconc.ge.cgp <- drugconc[is.element(drugconc[ ,"cellid"], myx), , drop=FALSE]
    mutation.cgp <- mutation[myx, , drop=FALSE]

    ## make sure that cellid are not factors
    celline.cgp[, "cellid"] <- as.character(celline.cgp[, "cellid"])
    sampleinfo.ge.cgp[, "cellid"] <- as.character(sampleinfo.ge.cgp[, "cellid"])
    drugpheno.ge.cgp[, "cellid"] <- as.character(drugpheno.ge.cgp[, "cellid"])
    drugconc.ge.cgp[, "cellid"] <- as.character(drugconc.ge.cgp[, "cellid"])

    message("Save data")
    write.csv(celline.cgp, file=file.path(tmpdir, "cell_line_collection_cgp.csv"), row.names=FALSE)
    write.csv(annot.ge.cgp, file=file.path(tmpdir, "annot_ge_cgp.csv"), row.names=FALSE)
    write.csv(t(data.ge.cgp), file=file.path(tmpdir, "data_ge_cgp.csv"))
    write.csv(drugpheno.ge.cgp, file=file.path(tmpdir, "drugpheno_ge_cgp.csv"), row.names=FALSE)
    write.csv(druginfo.ge.cgp, file=file.path(tmpdir, "druginfo_ge_cgp.csv"), row.names=FALSE)
    write.csv(drugconc.ge.cgp, file=file.path(tmpdir, "drugconc_ge_cgp.csv"), row.names=FALSE)
    write.csv(sampleinfo.ge.cgp, file=file.path(tmpdir, "sampleinfo_ge_cgp.csv"), row.names=FALSE)
    write.csv(mutation.cgp, file=file.path(tmpdir, "mutation_cgp.csv"), row.names=TRUE)
    save(list=c("data.ge.cgp", "annot.ge.cgp", "sampleinfo.ge.cgp", "mutation.cgp", "drugpheno.ge.cgp", "druginfo.ge.cgp", "drugconc.ge.cgp", "celline.cgp"), compress=TRUE, file=myfn)
  }
  
}



## End


