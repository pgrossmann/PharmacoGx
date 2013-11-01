########################
## Benjamin Haibe-Kains
## All rights Reserved
## October 23, 2013
########################

geneDrugAssocs <- 
function (x, y, z, method=c("lm", "cindex")) {
  ## function to compute the geen-drug associations
  ## Input:
  ##  x: gene expressions
  ##  y: drug sensitivity data
  ##  z: tissu type
  ##
  ## Output:
  ##  a vector of statistics
  
  method <- match.arg(method)
  switch(method, 
    "lm" = {
      ## standardized coefficient using linear regression controlled for tissue type
      xx <- as.numeric(x)
      names(xx) <- names(x)
      ## avoid underestimation of sd due to truncated or placeholders values
      xxsd <- sd(xx[!duplicated(xx)], na.rm=TRUE)
      x <- xx / xxsd
      if(!is.factor(y)) {
        yy <- as.numeric(y)
        names(yy) <- names(y)
        ## avoid underestimation of sd due to truncated or placeholders values
        yysd <- sd(yy[!duplicated(yy)], na.rm=TRUE)
        y <- yy / yysd
      } else {
        if(sum(table(y) > 0) != 2) { stop("y must have 2 levels!") }
      }
      ff <- sprintf("y ~ x")
      nn <- sum(complete.cases(x, y))
      if(sum(table(z) > 0) > 1) {
        ff <- sprintf("%s + z", ff)
        nn <- sum(complete.cases(x, y, z))
      }
      if(nn >= 3) {
        cc <- summary(glm(formula(ff), family=ifelse(is.factor(y), "binomial", "gaussian")))
        if(is.element("x", rownames(cc$coefficients))) {
          res <- c("estimate"=cc$coefficients["x", 1], "se"=cc$coefficients["x", 2], "n"=nn, "stat"=cc$coefficients["x", 3], "pvalue"=cc$coefficients["x", 4])
          } else { res <- c("estimate"=NA, "se"=NA, "n"=0, "stat"=NA, "pvalue"=NA) }
      } else { res <- c("estimate"=NA, "se"=NA, "n"=nn, "stat"=NA, "pvalue"=NA) }
    },
    "cindex" = {
      ## Somers' D index stratified by tissue type
      requite(survcomp)
      cc <- survcomp::concordance.index(x=-x, cl=y, strat=z, outx=TRUE, method="noether", na.rm=TRUE)
      ss <- (cc$c.index - 0.5) / cc$se
      res <- c("estimate"=(cc$c.index - 0.5) * 2, "se"=cc$se * 2, "n"=cc$n, "stat"=ss, "pvalue"=cc$p.value)
    })
    return(res)
}

## End
