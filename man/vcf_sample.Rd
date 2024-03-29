% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcfSample.R
\name{vcf_sample}
\alias{vcf_sample}
\title{Load sample names from VCF file}
\usage{
vcf_sample(file)
}
\arguments{
\item{file}{The VCF file.}
}
\value{
Sample names available in the VCF file.
}
\description{
This function extract sample names from the VCF file. It is mandatory to compress the VCF file using bgzip and create an index using tabix.
}
