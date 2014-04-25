PharmacoGx
==========

R package to perform large-scale pharmacogenomic analysis of datasets stored in InSilicoDB

Dependencies:

Install the R/Bioconductior dependencies:

pp <- c("InSilicoDb", "Biobase", "BiocGenerics", "org.Hs.eg.db", "survival", "survcomp", "genefu", "mRMRe", "WriteXLS")
source("http://bioconductor.org/biocLite.R")
myrepos <- biocinstallRepos()
rr <- biocLite(pkgs=pp, dependencies=TRUE, type="source", destdir=".")