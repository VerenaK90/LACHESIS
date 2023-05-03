#' Extract the copy number information per SSNV
#'
#' This function takes as input a vcf file and assigns each mutation a copy number state based on the output by ACEseq.
#' @param vcf Mutation file in VCF representation as a list (as returned by read.vcf from package bedR).
#' @param cnv ACEseq file as data frame
#' @return a list containing the coverage ratio, B-allele frequency, genotype, total copy number and A allele count per mutation
#' Extract.copy.number.info.per.SSNV()
#' @export


Extract.copy.number.info.per.SSNV <- function(vcf, cnv){

  res <- t(apply(vcf$vcf, 1, function(x){
    x <- unlist(x)
    tmp <- cnv[cnv$X.chromosome==x[1],]
    tmp <- tmp[as.numeric(tmp$start) <= as.numeric(x[2]) & as.numeric(tmp$end) >= as.numeric(x[2]),]
    if(nrow(tmp)>0){
      coverage.ratio <- unlist(tmp[28])
      baf <- unlist(tmp[30])
      genotype <- unlist(tmp[36])
      tcn <- unlist(tmp[37])
      A <- unlist(tmp[34])
    }else{
      coverage.ratio <- NA
      baf <- NA
      genotype <- NA
      tcn <- NA
      A <- NA
    }
    c(coverage.ratio = coverage.ratio, baf = baf, genotype = genotype, tcn = tcn, A=A)
  }))

  res <- as.data.frame(res)
  colnames(res) <- c("coverage.ratio", "baf", "genotype", "tcn", "A")
  res$coverage.ratio<- as.numeric(as.character(res$coverage.ratio))
  res$baf <- as.numeric(res$baf)
  res$tcn <- as.numeric(res$tcn)
  res$A <- as.numeric(res$A)
  res

}
