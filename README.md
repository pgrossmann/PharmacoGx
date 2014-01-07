PharmacoGx
==========

R package to perform large-scale pharmacogenomic analysis of datasets stored in InSilicoDB

Dependencies:

Install the R/Bioconductior dependencies:

pp <- c("Biobase", "BiocGenerics", "org.Hs.eg.db", "survival", "survcomp", "genefu", "mRMRe", "WriteXLS")
source("http://bioconductor.org/biocLite.R")
myrepos <- biocinstallRepos()
rr <- biocLite(pkgs=pp, dependencies=TRUE, type="source", destdir=".")

Download and install the standalone packages:

- jetset: http://goo.gl/mI5SW9
- InSilicoDb2: http://goo.gl/FWHt8W

Then run the following commands in your R session

pp <- c("inSilicoDb2_2.0.0.tar.gz", "jetset_1.6.0.tar.gz")
install.packages(pkgs=pp, repos=NULL, type="source")