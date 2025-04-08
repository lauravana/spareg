# spareg

**spareg** is an R package for building ensembles of predictive generalized linear models tailored for high-dimensional data. It leverages a powerful combination of **variable screening** and **random projection** techniques to efficiently tackle large-scale modeling problems.

[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![CRAN Status](https://www.r-pkg.org/badges/version/spareg)](https://CRAN.R-project.org/package=spareg)
[![R-CMD-check](https://github.com/lauravana/spareg/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/lauravana/spareg/actions)

---

## Features

- Efficient modeling with **large predictor sets**
- Modular S3 architecture for screening and projection techniques
- Competitive predictive performance with **low computational cost**
- Extensible framework for custom screening and projection methods
- Simple, user-friendly API for rapid experimentation

---

## Installation

You can install the development version of `spareg` from GitHub:

```r
# install.packages("devtools")
devtools::install_github("laura/spareg")
```

The CRAN version can be installed with 

```r
install.packages("spareg")
```
