
# geno2r

<!-- badges: start -->
<!-- badges: end -->

The goal of geno2r is to easily import genotype data from various formats into R. 
- [RSamtools](https://bioconductor.org/packages/release/bioc/html/Rsamtools.html) package imports wrapped [htslib](https://github.com/samtools/htslib/) to enhance the VCF reading speed. In order to efficiently read a VCF file, it should be compressed with [bgzip](https://www.htslib.org/doc/bgzip.html) and indexed with [tabix](https://www.htslib.org/doc/tabix.html).
- The [genio](https://cran.r-project.org/web/packages/genio/index.html) package offers functions that enable the reading of [PLINK binary format data](https://www.cog-genomics.org/plink/).



## Installation

You can install the development version of geno2r from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("oumarkme/geno2r")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(geno2r)
read_vcf(file = "example.vcf.gz", range="chr1:1000-1500")
```

