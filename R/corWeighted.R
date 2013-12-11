########################
## Benjamin Haibe-Kains
## All rights Reserved
## October 23, 2013
########################

#################################################
## Compute a weighted correlation coefficient
## inspired from package boot
##
## inputs:	
##      - 
##
## outputs, a list of items:
##			- 
##
#################################################

`corWeighted` <- 
function (x, y, w, method=c("pearson", "spearman"), alternative=c("two.sided", "greater", "less"), permutation.test=FALSE, nperm=1000, nthread=1, setseed, na.rm=FALSE) {
 
  if(missing(w)) { w <- rep(1, nrow(d))/nrow(d) }
  method <- match.arg(method)
  if(method == "spearman") { d <- apply(d, 2, rank) }
  alternative <- match.arg(alternative)
  
  wcor <- function (d, w) {
    s <- sum(w)
    m1 <- sum(d[, 1L] * w)/s
    m2 <- sum(d[, 2L] * w)/s
    res <- (sum(d[, 1L] * d[, 2L] * w)/s - m1 * m2)/sqrt((sum(d[, 1L]^2 * w)/s - m1^2) * (sum(d[, 2L]^2 * w)/s - m2^2))
    return(res)
  }
  
  ## remove missing values
  ccix <- complete.cases(x, y, w)
  if(!all(ccix) && !na.rm) { warning("Missing values are present") }
  if(sum(ccix) < 3) {
    res <- NA
    if(permutation.test) { res <- list("rho"=NA, "p"=NA) }
    return(res)
  }
  x <- x[ccix]
  y <- y[ccix]
  w <- w[ccix]
  
  res <- wcor(d=cbind(x, y), w=w)
  
  p <- NULL
  if (permutation.test) {
    require(parallel)
    if (!missing(setseed)) { set.seed(setseed) }
    splitix <- parallel::splitIndices(nx=nperm, ncl=nthread)
    if (!is.list(splitix)) { splitix <- list(splitix) }
    splitix <- splitix[sapply(splitix, length) > 0]
    mcres <- parallel::mclapply(splitix, function(x, xx, yy, ww) {
      pres <- sapply(x, function(x, xx, yy, www) {
        ## permute the data and the weights
        d2 <- cbind(xx[[sample(1:length(xx))]], yy[sample(1:length(yy))])
        w2 <- ww[sample(1:length(ww))]
        return(wcor(d=d2, w=w2))
      }, xx=xx, yy=yy, ww=ww)
      return(pres)
    }, xx=x, yy=y, ww=w)
    perms <- do.call(c, mcres)
    
    switch (alternative,
      "two.sided" = { 
        if (res < 0) { p <- sum(perms <= res) } else { p <- sum(perms >= res) }
        if (p == 0) { p <- 1 / (nperm + 1) } else { p <- p / nperm }
        p <- p * 2
      },
      "greater" = {
        p <- sum(perms >= res) 
        if (p == 0) { p <- 1 / (nperm + 1) } else { p <- p / nperm }
      },
      "less" = {
        p <- sum(perms <= res) 
        if (p == 0) { p <- 1 / (nperm + 1) } else { p <- p / nperm }
      })
    res <- list("rho"=res, "p"=p)
  }
  return(res)
}

## End
