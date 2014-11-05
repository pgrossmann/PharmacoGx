######################################################################################################
###
### Function to download, normalize, and curate TGGATES data sets
###
### Author: Patrick Grossmann, patrick@jimmy.harvard.edu, pgrossmann@github.com/, Smisery1@github.com
###
### Mar 24, 2014
###
#######################################################################################################


#' Downloads toxicogenomic data compromising drug information and gene expression response data 
#' from the TG-GATEs databse.
#' 
#' Gene expression is normalized by RMA algorithm, and saved to an Biobase expression set.
#' Should the gene expression data hosted by Array Express be already bundled in an expression set ("*eSet_orig.rda"),
#' expression data is not re-downloaded. Curation is always re-done with current sdrf file anyway.
#' SDRF file MUST be provided or re-downloaded within a total dataset!
#' 
#' @param identifier [character] accession number in Array Express
#' @param outdidr [character] directory to work in (download, etc.) and output results
#' @param unzip [logical] if downloaded zip files are required to be unzipped
#' @param sourcedir [character] path to local sourcedir, if data has been downloaded already, 
#' or NULL otherwise
#' @param verbose [logical] print additional processing information
#' @param resultPref [character] prefix or result file. will be appended by ".rda"
#' @param accessions may contain more than one tissue type from the same sample (e.g. 799);
#' in this case tissue will determine which one to return to get unique sample names.
#' @param return [ExpressionSet] expression set object, this object is named <resultPref>
normalize.TGGATES <- function(identifier, outdir=sprintf("normalizeTGGATES_%s",identifier), 
                              unzip=TRUE, sourcedir=NULL, verbose=FALSE, 
                              resultPref=sprintf("%s.eSet",identifier), tissue=c("liver","kidney")) {
  #-----------------------------------------------------------------
  ### setup code
  #-----------------------------------------------------------------
  
  require(ArrayExpress) # for getAE and Biobase 
  require(rat2302.db)
  require(hgu133plus2.db)
  require(BufferedMatrixMethods) # rma
  require(affyio) # for BufferedMatrix.justRMA() (bug in BufferedMatrix.justRMA ?)
  
  ### global variables ###
  resultPref <- gsub("-","",resultPref)
  tissue <- match.arg(tissue)
  
  ### aux. functions ###
  
  pVerbose <- function(message) { if (verbose) print(message) }

  #' function NOT used atm!
  #' @param samples [data.frame/table] subtable of duplica
  #' @param return row
  duplicate.samples <- function(samples) {
    stop("Don't use this function now")
    samples[1, , drop=F]
  }
  
  #' DEAL with duplicate samples.
  #'
  #' @param phenotable [data.frame/table] full phenotype table (sdrf file)
  #' @param return [data.frame/table] subset of phenotype table without duplicates
  duplicate.pheno <- function(phenotable) {
    dups <- duplicated(phenotable[,"Source.Name"])
    pVerbose(sprintf("%s duplicated entries are removed",sum(dups)))
    phenotable[!dups, , drop=F]
  }
  
  #' TODO: insource from below
  #' DEAL with duplicate genes.
  #'
  duplicate.genes <- function(geneexpressions) {
    stop("Don't use this function now")
  }
  
  #' TODO: insource from below
  #' PERFORM mapping
  #'
  map.genes <- function(probelist,identifier) {
    stop("Don't use this function now")
  }
  
  #' Datasets may have multiple tissues (liver, kidney).
  #' This function grabs source names for only a single tissue
  #' and if duplicated source names the first one is taken. 
  #' Duplicated source names can occure BEFORE applying this function
  #' as multiple samples can be extracted for an individual rat 
  #' (e.g. both liver and kidney)
  take.tissue <- function(phenotable, tissue) {
    tt <- phenotable$Factor.Value.OrganismPart.
    pheno.sub <- phenotable[tt==tissue, ,drop=F]
    if (nrow(pheno.sub) == 0)
      stop(sprintf("Data does not contain any %s data"), tissue)
    pheno.sub
  }
  
  #-----------------------------------------------------------------
  ### setup data
  #-----------------------------------------------------------------
  message("setting up data")
  
  # initially assume data needs to be downloaded
  local <- FALSE
  
  if (!is.null(sourcedir)) {
    if (!file.exists(sourcedir)) stop("sourcedir does not exist!")
    pVerbose(sprintf("will be reading data from sourcedir %s", sourcedir))
    local <- TRUE # if data has been downloaded already
  } else {
    pVerbose(sprintf("will download data to outdir %s", outdir))
    sourcedir <- outdir
      # outdir
    if (!file.exists(sourcedir)) dir.create(sourcedir)
  }
  
  # don't use outdir from here on #
  
  ## download data if not done yet ##
  pVerbose("processing raw array data")
  getAE(identifier, type="full", extract=unzip, sourcedir=sourcedir, local=local, path=sourcedir)
  
  ## get important file names ##
  pVerbose("listing meta data")
  adfNames <- list.files(sourcedir, pattern="adf.txt")
  sdrfNames <- list.files(sourcedir, pattern="sdrf.txt")
  idfNames <- list.files(sourcedir, pattern="idf.txt")
  celNames <- list.files(sourcedir, pattern="CEL")
  
  # print file names if verbose
  pVerbose("meta data file names:")
  sapply(c(adfNames, sdrfNames, idfNames), function(x) pVerbose(x))
  
  ## check if meta file names make sense
  stopifnot(sapply(c(adfNames, sdrfNames, idfNames), length) == 1)
  
  adf <- file.path(sourcedir, adfNames)
  sdrf <- file.path(sourcedir, sdrfNames)
  idf <- file.path(sourcedir, idfNames)
  #cel <- file.path(prependCEL, celNames)
  
  ###################################################################
  ### read data and normalization (rma)
  ### The magic starts here ! ###
  ###################################################################
  message("read in phenotype data and expression data")
 
  ## read in phenotype data ##
  pVerbose("reading in phenotype data")
  pheno <- read.csv(sdrf, sep="\t", stringsAsFactors=F)
  
  ## read in feature data ##
  #pVerbose("reading in feature data")
  #featu <- read.table(adf, sep="\t", stringsAsFactors=F)
  
  ## read in experiment data ##
#   pVerbose("reading in experiment data")
  #experiment <- read.table(idf, sep="\t", stringsAsFactors=F)

  ## read in, normalize, and save initial expression data ##
  pVerbose("reading in expression data")
  eSet.orig.fn <- file.path(sourcedir, sprintf("%s_eSet_orig.rda",identifier))
  if (!file.exists(eSet.orig.fn)) {
    message("creating normalized expression set")
    ## this is kind of messy! but cannot be avoided as justRMA prepends getwd()....
    wd <- getwd() # use this to set and unset wd
    tryCatch({
        setwd(sourcedir) # go to sourcedir because justRMA prepends getwd()
        #eSet.orig <- affy::justRMA(filenames=celNames)
        eSet.orig <- BufferedMatrixMethods::BufferedMatrix.justRMA(filenames=celNames)
      },
      error=function(e) {
        setwd(wd)
        print(e)
        stop("something wrong with a CEL file!!")
      }
    )
    setwd(wd)
    save(eSet.orig, file=eSet.orig.fn)
    message(sprintf("eSet successfully stored to %s", eSet.orig.fn))
  } else load(eSet.orig.fn)
  
  #-----------------------------------------------------------------
  ### curation, depending on organism (rat, human, etc)
  #-----------------------------------------------------------------
  message("curating data")  
  
  # work on expression data seperately; consider to rm(eSet.orig}
  geneex <- exprs(eSet.orig)
  
  ## curate phenotype data ##
  
  ## if accession contains more than one tissue type
  containsOneTissue <- ifelse(!is.null(pheno$Factor.Value.OrganismPart.), FALSE, TRUE)
  if (!containsOneTissue) {
    pVerbose("getting data for only one tissue")
    pheno <- take.tissue(pheno, tissue=tissue)
  } else pheno <- pheno
  pVerbose("remove phenotype duplicates")
  pheno.unique <- duplicate.pheno(pheno)
  rownames(pheno.unique) <- pheno.unique[,"Source.Name"] # to also check if really unique
#   print(sprintf("%s duplicate sample(s) names have been removed",nrow(pheno)-nrow(pheno.unique)))
  # keep unique expressions
  geneex <- geneex[,pheno.unique[,"Array.Data.File"]]
  # assign sample phenotype labels to expression data instead of CEL file names
  colnames(geneex) <- 
  pheno.unique[match(colnames(geneex),pheno.unique[,"Array.Data.File"]),"Source.Name"]
  
  ## map probes to gene IDs ##
  
  pVerbose("chose chip database")
  db <-
  switch(identifier,
    "E-MTAB-797"=rat2302.db,
    "E-MTAB-798"=hgu133plus2.db,
    "E-MTAB-799"=rat2302.db,
    "E-MTAB-800"=rat2302.db)
  
  pVerbose("map probes to entrez gene ids")
  gids <- AnnotationDbi::select(db,rownames(geneex),col=c("ENTREZID","PROBEID"))
  
  ## remove probes not matchable to gene id
  gids.fullmatched <- gids[!is.na(gids[,"ENTREZID"]), , drop=F]
  print(sprintf("%s probes did not match an entrez gene and were removed",nrow(gids)-nrow(gids.fullmatched)))
  
  ## select unique genes
  pVerbose("select unique genes")
  probes2gids <- gids.fullmatched[,"ENTREZID"]
  names(probes2gids) <- gids.fullmatched[,"PROBEID"]
  pVerbose("geneid.map(probes2gids,t(geneex),probes2gids) performed")
  geneMapResults <- genefu::geneid.map(probes2gids,t(geneex)[,names(probes2gids)],probes2gids) # performs unique mapping
  gids.fullmatched.unique <- geneMapResults$geneid1[!is.na(geneMapResults$geneid1)] # named vector (names are probe IDs) of gene IDs
  print(sprintf("Of %s probes, %s were unique entrez genes were selected", length(probes2gids), length(gids.fullmatched.unique)))
  
  ## assign gene ids to expression data and subset accordingly ##
  pVerbose("subsetting expression data") # by sample IDs
  geneex <- geneex[names(gids.fullmatched.unique),]
  idNames <- paste0("geneid.", gids.fullmatched.unique) # use gene IDs instead of probe IDs
  rownames(geneex) <- idNames
  pVerbose("creating expression set with full annotation")
  varmetadata <- data.frame(labelDescription=colnames(pheno.unique), row.names=colnames(pheno.unique))
  pd <- new("AnnotatedDataFrame", data=pheno.unique, varMetadata=varmetadata)
  eSet.processed <- new("ExpressionSet", exprs=geneex, phenoData=pd, annotation=annotation(eSet.orig))
  ## save expression set to rdata
  if (!containsOneTissue) {
    resultPref <- sprintf("%s_%s",resultPref,tissue)
  }
  fullEset.fn <- file.path(sourcedir,sprintf("%s.rda",resultPref))
  eSetName <- resultPref
  assign(eSetName, eSet.processed)
  save(list=c(eSetName), file=fullEset.fn, compress=TRUE)
  message(sprintf("sucessfully saved full processed %s data to %s", identifier, fullEset.fn))
  get(eSetName) # return eSet
}


# normalize.TGGATES("E-MTAB-797", sourcedir="tmp", unzip=F,verbose=T)