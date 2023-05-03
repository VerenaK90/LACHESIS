#' Extracts purity and ploidy estimates from ACEseq file name
#'
#' @param aceseq.filename The name of the ACEseq file.
#' @return purity and ploidy.
#' Extract.purity.ploidy.from.ACEseq()
#' @export

Extract.purity.ploidy.from.ACEseq <- function(aceseq.filename){
  ploidy. <- sapply(aceseq.filename, function(x){strsplit(x, split="extra")[[1]][2]})
  purity. <- sapply(ploidy., function(x){strsplit(x, split="_")[[1]][2]})
  purity. <- as.numeric(sapply(purity., function(x){strsplit(x, split=".txt")[[1]][1]}))
  ploidy. <- as.numeric(sapply(ploidy., function(x){strsplit(x, split="_")[[1]][1]}))

  return(list(purity=purity., ploidy = ploidy.))
}
