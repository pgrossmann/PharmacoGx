########################
## Benjamin Haibe-Kains
## All rights Reserved
## October 23, 2013
########################

#################################################
## Normalize and import all data from Connectivity Map
##
## inputs:	
##      - data: gene expression data matrix
##			- drug: single or vector of drug(s) of interest; if a vector of drugs is provided, they will be considered as being the same drug and will be jointly analyszed
##			- drug.id: drug used in each experiment
##			- drug.concentration: drug concentration used in each experiment
##			- cell: cell line for each experiment
##			- xptype: type of experiment (perturbation or control)
##      - batch: 
##			- single: Shoudl the statitsics be computed for each cell line individually?
##      - nthread: number of parallel threads (bound to the maximum number of cores available)
##
## outputs, a list of 2 items:
##			- meta: a data.frame including the statistics (standardidzed mean difference Hedges' g, its standard error, sample size for class 1 and class 2, p-value, p-value for heterogeneity between cell lines and probe annotations)
##			- list of data.frame with similar results for each cell line separately if any
##
## Notes:	duration is not taken into account as only 4 perturbations lasted 12h, the other 6096 lasted 6h
#################################################

`rankgenesCMAP` <- 
function (data, drug, drug.id, drug.concentration, cell, xptype, batch, single.cell=FALSE, nthread=1, verbose=FALSE) {
  if (nthread != 1) {
    require(parallel)
    availcore <- parallel::detectCores()
    if (missing(nthread) || nthread < 1 || nthread > availcore) {
      nthread <- availcore
    }
  }
	if (any(c(length(drug.id), length(drug.concentration), length(cell), length(xptype), length(batch)) != nrow(data))) {
    stop("length of drug.id, drug.concentration, cell, xptype and batch should be equal to the number of rows of data!")
  }
	names(drug.id) <- names(drug.concentration) <- names(cell) <- names(batch) <- rownames(data)
	if (!all(complete.cases(cell, xptype, batch))) {
    stop("cell, batch, and xptype should not contain missing values!")
  }
  ## is the drug in CMAP?
  drugix <- drug.id %in% drug
  if (sum(drugix) == 0) {
    warning(sprintf("Drug(s) %s not in CMAP", paste(drug, collapse=", ")))
    return(list("all.cell"=NULL, "single.cell"=NULL))
  }
	## select xps with controls or with the drug(s) of interest
	iix <- is.na(drug.id) | drugix
	data <- data[iix, ,drop=FALSE]
	drug.id <- drug.id[iix]
	drug.concentration <- drug.concentration[iix]
	cell <- cell[iix]
	xptype <- xptype[iix]
  batch <- batch[iix]
	
	res.cell <- NULL
  
  ## build input matrix
  inpumat <- NULL
	## for each batch/vehicle of perturbations+controls (test within each batch/vehicle to avoid batch effect)
  ubatch <- sort(unique(batch[!is.na(xptype) & xptype == "perturbation"]))
  names(ubatch) <- paste("batch", ubatch, sep="")
  
  for (bb in 1:length(ubatch)) {
		## identify the perturbations and corresponding control experiments
    xpix <- rownames(data)[complete.cases(batch, xptype) & batch == ubatch[bb] & xptype == "perturbation"]
    ctrlix <- rownames(data)[complete.cases(batch, xptype) & batch == ubatch[bb] & xptype == "control"]
    
		if (all(!is.na(c(xpix, ctrlix))) && length(xpix) > 0 && length(ctrlix) > 0) {
      if (!all(is.element(ctrlix, rownames(data)))) {
        stop("data for some control experiments are missing!")
      }
			if (verbose) {
        cat(sprintf("cell %s: batch %i/%i -> %i vs %i\n", ucell[cc], bb, length(ubatch), length(xpix), length(ctrlix)))
      }
      ## transformation of drug concentrations values
      # drug.concentration[!is.na(drug.concentration) & drug.concentration != 0] * 10^6
      conc <- drug.concentration * 10^6
      inpumat <- rbind(inpumat, data.frame("treated"=c(rep(1, length(xpix)), rep(0, length(ctrlix))), "cell"=c(cell[xpix], cell[ctrlix]), "batch"=paste("batch", c(batch[xpix], batch[ctrlix]), sep=""), "concentration"=c(conc[xpix], conc[ctrlix])))
		}
	}
  inpumat[ , "cell"] <- factor(inpumat[ , "cell"], ordered=FALSE)
  inpumat[ , "batch"] <- factor(inpumat[ , "batch"], ordered=FALSE)
  
  if (nrow(inpumat) < 3 || length(sort(unique(inpumat[ , "concentration"]))) < 2) {
    ## not enough experiments in CMAP
    warning(sprintf("Not enough data for drug(s) %s", paste(drug, collapse=", ")))
    return(list("all.cell"=NULL, "single.cell"=NULL))
  }
  
  ## linear model for cmap data
  cmap.lm <- function(x, concentration, cell, batch) {
    nc <- c("estimate", "se", "n", "tstat", "fstat", "pvalue")
    if (length(sort(unique(concentration))) < 2) {
      warning("No drug concentrations tested")
      tt <- rep(NA, length(nc))
      names(tt) <- nc
      return(tt)
    }
    ff0 <- sprintf("x ~ 1")
    ff <- sprintf("%s + concentration", ff0)
    if (length(sort(unique(cell))) > 1) { 
      ff0 <- sprintf("%s + cell", ff0)
      ff <- sprintf("%s + cell", ff)
    }
    if (length(sort(unique(batch))) > 1) {
      ff0 <- sprintf("%s + batch", ff0)
      ff <- sprintf("%s + batch", ff)
    }
    dd <- data.frame("x"=x, "concentration"=concentration, "cell"=cell, "batch"=batch)
    nn <- sum(complete.cases(dd))
    if(nn < 3) {
      tt <- c("estimate"=NA, "se"=NA, "n"=nn, "tsat"=NA, "fstat"=NA, "pvalue"=NA)
    } else {
      mm0 <- lm(formula=ff0, data=dd, model=FALSE, x=FALSE, y=FALSE, qr=TRUE)
      mm <- lm(formula=ff, data=dd, model=FALSE, x=FALSE, y=FALSE, qr=TRUE)
      mmc <- anova(mm0, mm)
      mm <- summary(mm)
      ## extract statistics
      tt <- c("estimate"=mm$coefficients["concentration", "Estimate"], "se"=mm$coefficients["concentration", "Std. Error"], "n"=nn, "tsat"=mm$coefficients["concentration", "t value"], "fstat"=mmc$F[2], "pvalue"=mmc$'Pr(>F)'[2])
    }
    names(tt) <- nc
    return(tt)
  }
     
  res <- NULL
  ucell <- sort(unique(as.character(inpumat[ , "cell"])))
  lcell <- list("all"=ucell)
  if(single.cell) {
    lcell <- c(lcell, as.list(ucell))
    names(lcell)[-1] <- ucell
  }
  for(ll in 1:length(lcell)) {
    ## select the cell line of interest
    inpumat2 <- inpumat[!is.na(inpumat[ , "cell"]) & is.element(inpumat[ , "cell"], lcell[[ll]]), , drop=FALSE]
    inpumat2 <- inpumat2[complete.cases(inpumat2), , drop=FALSE]
    if (nrow(inpumat2) < 3 || length(sort(unique(inpumat2[ , "concentration"]))) < 2) {
      ## not enough experiments in CMAP
      nc <- c("estimate", "se", "n", "tstat", "fstat", "pvalue")
      rest <- matrix(NA, nrow=nrow(data), ncol=length(nc), dimnames=list(rownames(data), nc))
    } else {
  		## test perturbation vs control
      if(nthread > 1) {
        ## parallel threads
        splitix <- parallel::splitIndices(nx=ncol(data), ncl=nthread)
        splitix <- splitix[sapply(splitix, length) > 0]
        mcres <- parallel::mclapply(splitix, function(x, data, inpumat) {
    			res <- t(apply(data[rownames(inpumat), x, drop=FALSE], 2, cmap.lm, concentration=inpumat[ , "concentration"], cell=inpumat[ , "cell"], batch=inpumat[ , "batch"]))
          return(res)
    		}, data=data, inpumat=inpumat2)
        rest <- do.call(rbind, mcres)
      } else {
        ##################
        ## DEBUG
        # browser()
        # genec <- intersect(rownames(angiogenes), colnames(data))
        # aa <- angiogenes[genec, "NetworkClassification"]
        # rest <- t(apply(data[rownames(inpumat2), genec, drop=FALSE], 2, cmap.lm, concentration=inpumat2[ , "concentration"], cell=inpumat2[ , "cell"], batch=inpumat2[ , "batch"]))
        # rr <- apply(data[rownames(inpumat2), genec, drop=FALSE], 2, function(x, concentration) { return(mean(x[concentration != 0], na.rm=TRUE) - mean(x[concentration == 0], na.rm=TRUE)) }, concentration=inpumat2[ , "concentration"])
        # boxplot(rr ~ aa, outline=FALSE, main="Simple fold change")
        # X11()
        # boxplot(rest[ , "estimate"] ~ aa, outline=FALSE, main="Coefficient in linear model")
        ##################
        rest <- t(apply(data[rownames(inpumat2), , drop=FALSE], 2, cmap.lm, concentration=inpumat2[ , "concentration"], cell=inpumat2[ , "cell"], batch=inpumat2[ , "batch"]))
      }
    }
    rest <- cbind(rest, "fdr"=p.adjust(rest[ , "pvalue"], method="fdr"))
    res <- c(res, list(rest))
  }
  names(res) <- names(lcell)
	return(res)
}

## End
