#' Converts a user-specified bed-file with copy number information into a standardized format that can be used as input for downstream analysis
#'
#' @author {Verena Körber}
#' @description Convert a user-specified bed-file with copy number information into a standardized format. Perform various quality checks on the file input and return the clean and standardized data-frame. Written on May 31st, 2023.
#'
#' @param cn.info A data frame containing the copy number information. Requires columns for the chromosome number, start and end of the segment, the number of A- and B-alleles
#' @param chr.col name or number of the column containing the chromosome number
#' @param start.col name or number of the column containing the first position of the segment
#' @param end.col name or number of the column containing the last position of the segment
#' @param A.col name or number of the column containing the number of A alleles
#' @param B.col name or number of the column containing the number of B alleles
#' @param merge.tolerance the maximal distance below which adjacent segments with equal copy number are merged. Defaults to 10^5 bp
#' @return A standardized data frame with copy number information per segment.
#' readCNV()
#' @export

readCNV <- function(cn.info, chr.col, start.col, end.col, A.col, B.col, merge.tolerance = 10^5){

  if(!chr.col %in% colnames(cn.info)){
    stop("Error: chromosome information is not provided.")
  }else if(!start.col %in% colnames(cn.info)){
    stop("Error: start position is not provided.")
  }else if(!end.col %in% colnames(cn.info)){
    stop("Error: end position information is not provided.")
  }else if(!A.col %in% colnames(cn.info)){
    stop("Error: A allele counts are not provided.")
  }else if(!B.col %in% colnames(cn.info)){
    stop("Error: B allele counts are not provided.")
  }

  if(!(is.character(cn.info[,chr.col]) | is.numeric(cn.info[,chr.col]))){
    stop("Error: chromosome information must be string or numeric.")
  }else if(!is.numeric(cn.info[,start.col])){
    stop("Error: start position must be numeric.")
  }else if(!is.numeric(cn.info[,end.col])){
    stop("Error: end position must be numeric.")
  }

  ## in ACEseq output A and B alleles are characters - transform to numeric
  cn.info[,A.col] <- as.numeric(as.character(cn.info[,A.col]))
  cn.info[,B.col] <- as.numeric(as.character(cn.info[,B.col]))

  ## check chromosome format and amend if not 'chr1', 'chr2', etc.
  if(cn.info[1,chr.col]=="chr"){
    cn.info[,chr.col] <- gsub(pattern = "chr", replacement = "", x = cn.info[,chr.col])
  }

  ## subset on autosomes
  cn.info <- cn.info[as.numeric(cn.info[,chr.col]) %in% c(1:22),]

  ## order and re-name columns
  cn.info <- cn.info[,c(chr.col, start.col, end.col, A.col, B.col)]
  colnames(cn.info) <- c("Chr", "Start", "End", "A", "B")

  ## filter NA-segments
  cn.info <- cn.info[!is.na(cn.info$A) & !is.na(cn.info$B),]

  ## add total copy number
  cn.info$TCN <- cn.info$A + cn.info$B

  ## merge adjacent segments
  cn.info.merged <- cn.info[1,]
  i = 1
  for(j in 2:nrow(cn.info)){
    if(cn.info[j,"Chr"] != cn.info.merged[i,"Chr"]){
      cn.info.merged <- rbind(cn.info.merged, cn.info[j,])
      i <- i + 1
      next
    }
    if(cn.info[j,"Start"] - cn.info.merged[i,"End"] > merge.tolerance){
      cn.info.merged <- rbind(cn.info.merged, cn.info[j,])
      i <- i + 1
      next
    }
    if(cn.info[j,"A"] != cn.info.merged[i,"A"] | cn.info[j,"B"] != cn.info.merged[i,"B"]){
      cn.info.merged <- rbind(cn.info.merged, cn.info[j,])
      i <- i + 1
      next
    }
    cn.info.merged[i,"End"] <- cn.info[j,"End"]
  }

  cn.info <- cn.info.merged
  rm(cn.info.merged)

  return(cn.info)

}
