
#' Read VCF to R
#'
#' @description
#' This function reads the VCF file into R and converts it to a data.table format. While genomic range is given, [HTSlib](https://github.com/samtools/htslib), a wrapped C library, reads directly from a region instead of loading the full VCF file to memory. It is mandatory to compress the VCF file using bgzip and create an index using tabix.
#'
#'
#' @param file The VCF file
#' @param range Genome range to read. e.g. chr1:150-2000. If you want to read all variants, leave the field blank or enter "NULL" in the input.
#' @import data.table
#' @import Rsamtools
#' @import stringr
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @return VCF file in data.table format.
#'
#' @export
read_vcf = function(file, range = NULL){

  # check if the vcf and index files exist
  if(! file.exists(file)){
    stop("VCF file not found.")
  }
  if(! file.exists(paste0(file, ".tbi"))){
    stop(paste("Index file not found.", paste0(file, ".tbi"), "should be placed in the same directory as the VCF file."))
  }

  # Headers of the VCF file
  vcfHeader = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", as.character(Rsamtools::scanBcfHeader(file)[[file]]$Sample))

  # Read whole VCF file
  if(is.null(range)){
    vcf = data.table::fread(file)
  }

  # Read particle
  if(! is.null(range)){

    # get range
    range = stringr::str_split(range, ":|-")[[1]]
    chr = range[1]
    start = as.numeric(range[2])
    end = ifelse(is.na(as.numeric(range[3])), as.numeric(range[2]), as.numeric(range[3]))

    # load file
    vcf = unlist(Rsamtools::scanTabix(file, param=GenomicRanges::GRanges(seqnames = "1", IRanges::IRanges(start, width = end-start+1))))
    if(length(vcf) == 1){
      vcf = data.table::as.data.table(t(unlist(stringr::str_split(vcf, "\t"))))
    }else{
      vcf = data.table::fread(text = vcf)
    }
  }

  # rename header
  colnames(vcf) = vcfHeader

  # return
  return(vcf)
}

