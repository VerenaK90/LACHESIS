#' Convert mutation counts into mutation times
#'
#' The number of clonal mutations found in a clonal peak does not necessarily correspond to "mutation time" as mutations are lost from the non-amplified clonal peak by a chromosomal gain. This function corrects for such losses, rendering mutation times comparable.
#' @param clonal.mutations A vector with clonal mutation counts for different copy number and allele states. Corresponds to a column of the clonal.mutation.matrix as returned by Count.clonal.mutations().
#' @param max.CN maximal copy number to be considered. Defaults to 4.  Which information is to be retrieved? Possible values are "readcounts", "depth", "VAF", "AA_change". Defaults to "readcounts".
#' Mutation.time.converter()
#' @export

Mutation.time.converter <- function(clonal.mutations, max.CN=4){

  ## indicator vectors for the copy numbers and B-alleles that were taken into account
  copy.number.indicator <- c(unlist(sapply(1:max.CN, function(x){rep(x, each=x)})))
  B.allele.indicator <- c(unlist(sapply(1:max.CN, function(x){1:x})))

  ## Add the number of mutations, which do not contribute to the monosomic peak due to amplification, on the respective number of alleles

  clonal.mutations[which(B.allele.indicator==1)] <- clonal.mutations[which(B.allele.indicator==1)] +
    sapply(unique(copy.number.indicator), function(x){
      tmp <- which(copy.number.indicator==x & B.allele.indicator!=1)
      if(length(tmp)==0){
        return(0)
      }else{
        return(sum(clonal.mutations[tmp]*B.allele.indicator[tmp]))
      }
    })


  ## divide by the copy number to obtain mutations per copy
  clonal.mutations[which(B.allele.indicator==1)] <- clonal.mutations[which(B.allele.indicator==1)]/seq(1,max.CN)
  ## divide by 2 if A=B
  clonal.mutations[which(B.allele.indicator==(copy.number.indicator-B.allele.indicator) & B.allele.indicator!=1)] <- clonal.mutations[which(B.allele.indicator==(copy.number.indicator-B.allele.indicator) & B.allele.indicator!=1)]/2
  clonal.mutations
}
