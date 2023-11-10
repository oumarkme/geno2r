
#' Read vcf to R
#'
#' @param file the VCF file
#' @param range Genome range to read. e.g. chr1:150-2000
#' @import data.table
#' @import Rsamtools
#' @import stringr
#' @import GenomicRanges
#' @import IRanges
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

