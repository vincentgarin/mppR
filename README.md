mppR: Multi-Parent Population QTL Analysis
====


## Overview

mppR is an R package to perform QTL analysis of experimental multi-parent populations. The population must be composed of crosses between a set of at least three parents (e.g. factorial design, 'diallel', or nested association mapping). The functions cover data processing, QTL detection, and results visualization.

## Installation

mppR has different branches: "master", "mppR_clusthaplo", and "mppR_CRAN". The
"master" branch allows to perform MPP mixed model QTL detection calling the asreml-R package and function parent_cluster.mppData that call the archived R package clusthaplo for parent clustering. The "mppR_clusthaplo" branch contains the function parent_cluster.mppData that call the archived R package clusthaplo for parent clustering. The branch "mppR_CRAN" do not contain the mixed models and the call to clusthaplo.

```
devtools::install_github("vincentgarin/mppR", ref = "mppR_clusthaplo")

```

## Usage

[vignette of the package](inst/doc/mppR_vignette.pdf)

# Travis

[![Travis-CI Build Status](https://travis-ci.org/vincentgarin/mppR.svg?branch=master,mppR_CRAN)](https://travis-ci.org/vincentgarin/mppR)