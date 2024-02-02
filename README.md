# geno2r

<!-- badges: start -->

[![R-CMD-check](https://github.com/oumarkme/geno2r/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/oumarkme/geno2r/actions/workflows/R-CMD-check.yaml) [![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)

<!-- badges: end -->

The goal of geno2r is to easily import genotype data from various formats into R.

[HTSlib](https://github.com/samtools/htslib/) is wrapped and imported by [Rsamtools](https://bioconductor.org/packages/release/bioc/html/Rsamtools.html), to enhance the VCF reading speed. In order to efficiently read a VCF file, it should be compressed with [bgzip](https://www.htslib.org/doc/bgzip.html) and indexed with [tabix](https://www.htslib.org/doc/tabix.html).

## Installation

You can install the most recent version of geno2r from GitHub with:

``` r
#install.packages("remotes")
remotes::install_github("oumarkme/geno2r", force=TRUE)
```

## Example

This is a basic example which shows you how to read VCF data to R with [`read_vcf()`](reference/read_vcf.html):

``` r
library(geno2r)
read_vcf(file = "human.vcf.gz", range="22:49379357")   # or
read_vcf(file = "human.vcf.gz", range="22:49379357-49379357")   # or
read_vcf(file = "human.vcf.gz", range=NULL)
```

## Example data

You can find an example VCF file on [Github](https://github.com/oumarkme/geno2r/tree/main/vcf_data_example).
