
#' @param identifier [character] accession number in Array Express
#' @param tmp [character] directory to work in (download, etc.) and output results
#' @param unzip [logical] if downloaded zip files are required to be unzipped
#' @param sourcedir [character] path to local sourcedir, if data has been downloaded already, 
#' or NULL otherwise
#' @param verbose [logical] print additional processing information
normalize.TGGATES <- function(identifier, tmpdir="tmp", unzip=TRUE, sourcedir=NULL, verbose=FALSE) {
  #-----------------------------------------------------------------
  ### setup code
  #-----------------------------------------------------------------
  
  require(ArrayExpress) # for getAE and Biobase 
  require(rat2302.db)
  
  ### aux. functions ###
  
  pVerbose <- function(message) { if (verbose) print(message) }

  #-----------------------------------------------------------------
  ### setup data
  #-----------------------------------------------------------------
  print("setting up data")
  
  # initially assume data needs to be downloaded
  local <- FALSE
  
  # tmpdir
  if (!file.exists(tmpdir)) dir.create(tmpdir)
  
  if (!is.null(sourcedir)) {
    if (!file.exists(sourcedir)) stop("sourcedir does not exist!")
    pVerbose("will be reading data from sourcedir")
    local <- TRUE # if data has been downloaded already
  } else {
    pVerbose("will download data to tmpdir")
    sourcedir <- tmpdir
  }
  
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
  
  stopifnot(sapply(c(adfNames, sdrfNames, idfNames), length) == 1)
  
  adf <- file.path(sourcedir, adfNames)
  sdrf <- file.path(sourcedir, sdrfNames)
  idf <- file.path(sourcedir, idfNames)
  #cel <- file.path(sourcedir, celNames)
  
  #-----------------------------------------------------------------
  ### read data and normalization (rma)
  #-----------------------------------------------------------------
  print("normalizing expression data")
 
  ## read in phenotype data ##
  pVerbose("reading in phenotype data")
  pheno <- read.csv(sdrf, sep="\t", stringsAsFactors=F)
  
  ## read in feature data ##
  #pVerbose("reading in feature data")
  #featu <- read.table(adf, sep="\t", stringsAsFactors=F)
  
  ## read in experiment data ##
  pVerbose("reading in experiment data")
  #experiment <- read.table(idf, sep="\t", stringsAsFactors=F)

  ## read in, normalize, and save initial expression data ##
  pVerbose("reading in expression data")
  eSet.orig.fn <- file.path(tmpdir, "eSet.orig.rda")
  if (!file.exists(eSet.orig.fn)) {
    message("creating normalized expression set")
    eSet.orig <- affy::justRMA(filenames=celNames)
    save(eSet.orig, file=eSet.orig.fn)
    message(sprintf("eSet successfully stored to %s", eSet.orig.fn))
  } else load(eSet.orig.fn)
  
  #-----------------------------------------------------------------
  ### curation, depending on organism (rat, human, etc)
  #-----------------------------------------------------------------
  print("curating data")  
  
  ## map probes to gene IDs ##
  
  db <-
  switch(identifier,
    "E-MTAB-797"=rat2302.db)
  
  ## assign correct phenotype labels
  
  #celnames2samplenames <- 
  #colnames(exprs) <- 
  
  
}


normalize.TGGATES("E-MTAB-797", verbose=T)