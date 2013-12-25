
#' @param identifier [character] accession number in Array Express
#' @param tmp [character] directory to work in (download, etc.)
#' @param unzip [logical] if downloaded zip files are required to be unzipped
#' @param sourcedir [character] path to local sourcedir, if data has been downloaded already, 
#' or NULL otherwise
normalize.TGGATES <- function(identifier, tmpdir="tmp", unzip=TRUE, sourcedir=NULL) {
  require(ArrayExpress) # for getAE and Biobase  

  #-----------------------------------------------------------------
  ### setup
  #-----------------------------------------------------------------
  
  # initially assume data needs to be downloaded
  local <- FALSE
  
  if (!is.null(sourcedir)) {
    if (!file.exists(sourcedir)) stop("sourcedir does not exist!")
    local <- TRUE # if data has been downloaded already
  } else {
    sourcedir <- tmpdir
  }
  
  ## download data if not done yet ##
  getAE(identifier, type="full", extract=unzip, sourcedir=sourcedir, local=local, path=sourcedir)
  
  ## get important file names ##
  adfNames <- list.files(sourcedir, pattern="adf.txt")
  sdrfNames <- list.files(sourcedir, pattern="sdrf.txt")
  idfNames <- list.files(sourcedir, pattern="idf.txt")
  celNames <- list.files(sourcedir, pattern="CEL")
  
  stopifnot(sapply(c(adfNames, sdrfNames, idfNames), length) == 1)
  
  adf <- file.path(sourcedir, adfNames)
  sdrf <- file.path(sourcedir, sdrfNames)
  idf <- file.path(sourcedir, idfNames)
  #cel <- file.path(sourcedir, celNames)
  
  #-----------------------------------------------------------------
  ### read data and normalization (rma)
  #-----------------------------------------------------------------
 
  ## read in phenotype data ##
  pheno <- read.csv(sdrf, sep="\t", stringsAsFactors=F)
  
  ## read in feature data ##
  #featu <- read.table(adf, sep="\t", stringsAsFactors=F)
  
  ## read in experiment data ##
  #experiment <- read.table(idf, sep="\t", stringsAsFactors=F)

  ## read in, normalize, and save initial expression data ##
  eSet.orig.fn <- file.path(sourcedir, "eSet.orig.rda")
  if (!file.exists(eSet.orig.fn)) {
  print("hi")
  browser()
    eSet.orig <- affy::justRMA(filenames=celNames)
    save(eSet.orig, file=eSet.orig.fn)
  } else load(eSet.orig.fn)
  
  #-----------------------------------------------------------------
  ### curation
  #-----------------------------------------------------------------
  
  ## assign correct phenotype labels ##
  
  #celnames2samplenames <- 
  #colnames(exprs) <- 
  
  
}