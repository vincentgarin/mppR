## Test environments
* win-builder(devel and release)
* ubuntu 16.04, R 3.4.4
* ubuntu 14.04 (on travis-ci), R 3.5.0

## R CMD check results
There were no ERRORs or WARNINGs.

There were 2 NOTEs

* checking CRAN incoming feasibility ... NOTE
  
  New submission

  Possibly mis-spelled words in DESCRIPTION:
    QTLs (13:69)

  Quantitative trait loci (QTLs) standard abreviaion. E.g. r/qtl package.
  
* checking package dependencies ... NOTE
  
  Package suggested but not available for checking: 'clusthaplo'

  the function parent_cluster.mppData call the CRAN archived package 'clusthaplo'
  https://cran.r-project.org/src/contrib/Archive/clusthaplo/clusthaplo_1.2.tar.gz
  
  As an alternative, the user can provide a matrix with the parent clustering
  values calculated with another method of his/her choice.

## Vignette
The vignette is built on Maintainer machine.

## Other
* I read the CRAN Repository Policy Version $Revision: 3874 $ 
