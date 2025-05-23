#' Converts a user-specified bed-file with copy number information into a standardized format that can be used as input for downstream analysis
#'
#' @author {Verena Körber}
#' @description Convert a user-specified bed-file with copy number information into a standardized format. Perform various quality checks on the file input and return the clean and standardized data-frame. If column identifiers for chromosomal positions and allele-specific copy number information are not provided, the function attempts to identify these columns based on standard nomenclature. If total copy number information is provided but allele-specific information is missing, the function assumes that the number of B alleles is the rounded off of half the total copy number.
#'
#' @param cn.info Path to the copy number information. Requires columns for the chromosome number, start and end of the segment, and either the total copy number or the number of A- and B-alleles
#' @param chr.col column index of chromosome number
#' @param start.col column index of first position of the segment
#' @param end.col column index of last position of the segment
#' @param A.col column index of the number of A alleles. If A and B are not provided, allele configuration are assumed as 1:1 for disomic, 2:1 for trisomic and 3:1 for tetrasomic regions.
#' @param B.col column index of the number of B alleles. If A and B are not provided, allele configuration are assumed as 1:1 for disomic, 2:1 for trisomic and 3:1 for tetrasomic regions.
#' @param tcn.col column index of the total copy number. Is computed to A + B if not provided.
#' @param max.cn maximum copy number to be included in the analysis. Defaults to 4.
#' @param merge.tolerance the maximum distance below which adjacent segments with equal copy number are merged. Defaults to 10^5 bp.
#' @param ignore.XY Ignore allosomes. Default TRUE
#' @param tumor.id Tumor ID, optional.
#' @examples
#' aceseq_cn = system.file("extdata", "ACESeq/NBE11_comb_pro_extra2.59_0.83.txt", package = "LACHESIS")
#' cn_data = readCNV(aceseq_cn)
#' ascat_cn = system.file("extdata", "ASCAT/S98.segments.txt", package = "LACHESIS")
#' cn_data = readCNV(ascat_cn)
#' @return A standardized data frame with copy number information per segment.
#' readCNV()
#' @import data.table
#' @importFrom graphics plot.new
#' @importFrom grDevices dev.off pdf
#' @importFrom stats cor density end plot.ecdf start
#' @export

readCNV <- function(cn.info = NULL, chr.col = NULL, start.col = NULL, end.col = NULL, A.col = NULL, B.col = NULL, tcn.col = NULL, merge.tolerance = 10^5, ignore.XY = TRUE, max.cn = 4, tumor.id = NULL){

  . <- Alt <- Chr <- Chromosome <- ECA_time_mean <- End_position <- ID <- MRCA_time_mean <- Ref <- Start <- Start_Position <- TCN <- chrom <- cnv.file <- end <- start <- t_alt_count <-t_depth <- t_ref_count <- t_vaf <- NULL

  ## Check input format
  if(is.null(cn.info)){
    stop("Error: missing cn.info! Please provide path to file with copy number information.")
  }

  ## Read cn.info and assume column index if not provided
  format.cnv <- .format_cnv_data(cn.info = cn.info, chr.col = chr.col, start.col = start.col, end.col = end.col, A.col = A.col, B.col = B.col, tcn.col = tcn.col)

  cn.info <- format.cnv$cn.info
  chr.col <- format.cnv$chr.col
  start.col <- format.cnv$start.col
  end.col <- format.cnv$end.col
  A.col <- format.cnv$A.col
  B.col <- format.cnv$B.col
  tcn.col <- format.cnv$tcn.col
  estimate.alleles <- format.cnv$estimate.alleles

  message("********** Read in ", nrow(cn.info), " segments with copy number information on ", length(unique(cn.info[[chr.col]])), " chromosomes.")

  ## round copy number to multiples of 0.2
  cn.info[[tcn.col]] <- as.numeric(as.character(cn.info[[tcn.col]]))
  ## for callers without subconal/clonal assignment, remove subclonals based on maximal deviation of 0.2
  cn.info <- cn.info[cn.info[[tcn.col]]%%1 <=  0.2 | cn.info[[tcn.col]]%%1 >= 0.8,,drop = FALSE]
  ## round the remaining copy numbers
  cn.info[[tcn.col]] <- round(cn.info[[tcn.col]])

  message("********** Removing ", sum(is.na(cn.info[[tcn.col]])), " segments without copy number information...")

  cn.info <- cn.info[!is.na(cn.info[[tcn.col]]),]

  if(nrow(cn.info)==0){
    stop("Error: no segments with copy number information provided.")
  }

  message("********** Removing ", sum(cn.info[[tcn.col]] > max.cn), " segments with copy number > ", max.cn, "...")

  cn.info <- cn.info[cn.info[[tcn.col]] <= max.cn & cn.info[[tcn.col]] > 0,]

  if(nrow(cn.info)==0){
    stop("Error: no segments with copy number information greater 0 and <= ", max.cn)
  }

  if(estimate.alleles){ # assume 1:1, 2:1, 2:2, ... configuration
    cn.info <- .estimate_alleles(cn.info, tcn.col)
  }

  ## in ACEseq output A and B alleles are characters - transform to numeric.
  ## moreover, round copy numbers to 1 digit (for e.g. input from PURPLE)
  cn.info[[A.col]] <- as.numeric(as.character(cn.info[[A.col]]))
  cn.info[[B.col]] <- as.numeric(as.character(cn.info[[B.col]]))

  ## for callers without subconal/clonal assignment, remove subclonals based on maximal deviation of 0.2 for the total CN (A+B), and based on a maximal deviation of 0.1 for each allele
  cn.info <- cn.info[cn.info[[A.col]]%%1 <=  0.1 | cn.info[[A.col]]%%1 >= 0.9,,drop = FALSE]
  cn.info <- cn.info[cn.info[[B.col]]%%1 <=  0.1 | cn.info[[B.col]]%%1 >= 0.9,,drop = FALSE]
  ## round the remaining copy numbers
  cn.info[[A.col]] <- round(cn.info[[A.col]])
  cn.info[[B.col]] <- round(cn.info[[B.col]])


  ## check chromosome format and amend if not 'chr1', 'chr2', etc.
  if (grepl("chr", cn.info[[chr.col]][1])){
    message("********** Change chromosome names to 1, 2, 3, ...")
    cn.info[[chr.col]] <- gsub(pattern = "chr", replacement = "", x = cn.info[[chr.col]])
  }

  ## subset on autosomes
  if(ignore.XY){
    message("********** Removing ", sum(cn.info[[chr.col]] %in% c("X", "Y")), " segments on allosomes...")
    cn.info <- cn.info[!cn.info[[chr.col]] %in% c("X", "Y"),]
    if(nrow(cn.info)==0){
      stop("No copy number information on autosomes!")
    }
  }

  ## order and re-name columns
  cn.info <- cn.info[, .(Chr = cn.info[[chr.col]],
                         Start = cn.info[[start.col]],
                         End = cn.info[[end.col]],
                         A = cn.info[[A.col]],
                         B = cn.info[[B.col]])]


  ## filter NA-segments
  message("********** Removing ", sum(is.na(cn.info[["A"]]) | is.na(cn.info[["B"]])) , " segments without B allele information...")
  cn.info <- cn.info[!is.na(cn.info[["A"]]) & !is.na(cn.info[["B"]])]

  ## add total copy number
  cn.info$TCN <- cn.info$A + cn.info$B

  ## merge adjacent segments
  message("********** Merging adjacent segments with equal copy number...")
  cn.info <- .merge_adjacent_segs(cn.info, merge.tolerance)
  message("********** Retaining ", nrow(cn.info), " segments with copy number information on ", length(unique(cn.info$Chr)), " chromosomes.")

  if(any(cn.info$End < cn.info$Start)){
    warning("Removing ", sum(cn.info$End < cn.info$Start), " segments with start > end")
    cn.info <- cn.info[cn.info$Start <= cn.info$End,]
  }
  if(sum(cn.info$End - cn.info$Start) < 3*10^8){
    warning("Less than 10% of the genome with valid copy number information.")
  }

  attr(cn.info, "ID") <- tumor.id

  return(cn.info)

}

.merge_adjacent_segs <- function(cn.info, merge.tolerance){
  cn.info <- split(cn.info, cn.info$Chr) # split by chromosome
  cn.info <- lapply(cn.info, function(x){
    to.keep <- 1 + which(!(x$Start[-1] - x$End[-length(x$End)] < merge.tolerance & x$A[-1] == x$A[-length(x$A)] & x$B[-1] == x$B[-length(x$B)])) # merge segments if they have the same allele counts and if start and end are less than the merge.tolerance apart from each other
    end <- c(x$End[to.keep - 1], x$End[length(x$End)])
    x <- x[c(1, to.keep),]
    x$End <- end
    x
  })
  cn.info <- do.call(rbind.data.frame, cn.info)
  return(cn.info)
}

.format_cnv_data <- function(cn.info = NULL, chr.col = NULL, start.col = NULL, end.col = NULL, A.col = NULL, B.col = NULL, tcn.col = NULL){

  if (is.data.frame(cn.info) || is.data.table(cn.info)) {
    cn.info <- data.table::setDT(cn.info)
    print("datatable")
  } else if (is.character(cn.info)) {
    if (!file.exists(cn.info)) {
      stop(paste("Please provide a valid CNV file path.", cn.info, "does not exist."))
    } else {
      cn.info <- data.table::fread(file = cn.info, sep = "\t", header = TRUE)
    }
  } else {
    stop("Please provide a file path, data.frame, or data.table for copy number information.")
  }

  estimate.alleles <- FALSE # will be set to TRUE if allele info is not provided, see below

  if(is.null(chr.col) || is.na(chr.col)){
    chr.col <- colnames(cn.info)[grepl("chr", colnames(cn.info), ignore.case = TRUE)] # try to match with standard nomenclature
    chr.col <- ifelse(length(chr.col) > 0, chr.col, 1)
    warning("No chromosome identifier provided, assuming ", chr.col)
  }else if(is.character(chr.col)){
    chr.col <- match.arg(arg = chr.col, choices = colnames(cn.info), several.ok = FALSE)
  }else if(is.numeric(chr.col)){
    if(chr.col > ncol(cn.info)){
      stop("Error: 'arg' should be between 1 and ", ncol(cn.info))
    }
  }else{
    stop("Error: 'arg' should be string or numeric.")
  }

  if(is.null(start.col) || is.na(start.col)){
    start.col <- colnames(cn.info)[grepl("start", colnames(cn.info), ignore.case = TRUE) | grepl("pos", colnames(cn.info), ignore.case = TRUE)] # try to match with standard nomenclature
    start.col <- ifelse(length(start.col) > 0, start.col[1], 2)
    warning("No start position identifier provided, assuming ", start.col)
  }else if(is.character(start.col)){
    start.col <- match.arg(arg = start.col, choices = colnames(cn.info), several.ok = FALSE)
  }else if(is.numeric(start.col)){
    if(start.col > ncol(cn.info)){
      stop("Error: 'arg' should be between 1 and ", ncol(cn.info))
    }
  }else{
    stop("Error: 'arg' should be string or numeric.")
  }

  if(is.null(end.col) || is.na(end.col)){
    end.col <- colnames(cn.info)[grepl("end", colnames(cn.info), ignore.case = TRUE) | grepl("pos2", colnames(cn.info), ignore.case = TRUE)] # try to match with standard nomenclature
    end.col <- ifelse(length(end.col) > 0, end.col[1], 3)
    warning("No end position identifier provided, assuming ", end.col)
  }else if(is.character(end.col)){
    end.col <- match.arg(arg = end.col, choices = colnames(cn.info), several.ok = FALSE)
  }else if(is.numeric(end.col)){
    if(end.col > ncol(cn.info)){
      stop("Error: 'arg' should be between 1 and ", ncol(cn.info))
    }
  }else{
    stop("Error: 'arg' should be string or numeric.")
  }


  if(is.null(A.col) || is.na(A.col)){
    A.col <- colnames(cn.info)[grepl("major", colnames(cn.info), ignore.case = TRUE)] # try to match with standard nomenclature
    if(length(A.col)==0){
      if("A" %in% colnames(cn.info)){
        A.col <- "A"  # assume standard nomenclature
        warning("A allele identifier not provided, assuming 'A'")
      }else if (!is.null(tcn.col) & !is.null(B.col)){
        cn.info$A <- cn.info[[tcn.col]] - cn.info[[B.col]]
        message("********** A allele identifier not provided, computing A = TCN - B.")
      }else{
        estimate.alleles <- TRUE # Assume alleles as standard 1:1, 2:1 or 2:2 configuration
        A.col <- "A"
        warning("Allele information is not provided and will be assumed 1:1 in disomic regions, 2:1 in trisomic regions, 2:2 in tetrasomic regions, ... .")
      }
    }else{
      A.col <- A.col[1]
      warning("A allele identifier not provided, assuming ", A.col)
    }
  }else if(is.character(A.col)){
    A.col <- match.arg(arg = A.col, choices = colnames(cn.info), several.ok = FALSE)
  }else if(is.numeric(A.col)){
    if(A.col > ncol(cn.info)){
      stop("Error: 'arg' should be between 1 and ", ncol(cn.info))
    }
  }else{
    stop("Error: 'arg' should be string or numeric.")
  }

  if(is.null(B.col) || is.na(B.col)){
    B.col <- colnames(cn.info)[grepl("minor", colnames(cn.info), ignore.case = TRUE)] # try to match with standard nomenclature
    if(length(B.col)==0){
      if("B" %in% colnames(cn.info)){
        B.col <- "B" # assume standard nomenclature
        warning("B allele identifier not provided, assuming 'B'")
      }else if (!is.null(tcn.col) & !estimate.alleles){
        cn.info$B <- cn.info[[tcn.col]] - cn.info[[A.col]]
        message("********** B allele identifier not provided, computing B = TCN - A.")
      }else{
        B.col <- "B"
      }
    }else{
      B.col <- B.col[1]
      warning("B allele identifier not provided, assuming ", B.col)
    }
  }else if(is.character(B.col)){
    B.col <- match.arg(arg = B.col, choices = colnames(cn.info), several.ok = FALSE)
  }else if(is.numeric(B.col)){
    if(B.col > ncol(cn.info)){
      stop("Error: 'arg' should be between 1 and ", ncol(cn.info))
    }
  }else{
    stop("Error: 'arg' should be string or numeric.")
  }

  if(is.null(tcn.col) || is.na(tcn.col)){
    if(!is.null(A.col) & !is.null(B.col)){
      cn.info$TCN <- as.numeric(cn.info[[A.col]]) + as.numeric(cn.info[[B.col]])
      tcn.col <- "TCN"
      message("********** Total copy number computed as A + B.")
    }else{
      tcn.col <- colnames(cn.info)[grepl("\\btcn\\b", colnames(cn.info), ignore.case = TRUE) | grepl("\\bcnt\\b", colnames(cn.info), ignore.case = TRUE) |
                                     grepl("\\bcopynumber\\b", colnames(cn.info), ignore.case = TRUE) | grepl("\\bcopy number\\b", colnames(cn.info), ignore.case = TRUE)] # try to match with standard nomenclature
      if(length(tcn.col)==0){
        stop("Error: TCN identifier is not provided and could not be inferred!")
      }
      tcn.col <- tcn.col[1]
      warning("TCN identifier is not provided, assuming ", tcn.col)
    }
  }else if(is.character(tcn.col)){
    tcn.col <- match.arg(arg = tcn.col, choices = colnames(cn.info), several.ok = FALSE)
  }else if(is.numeric(tcn.col)){
    if(tcn.col > ncol(cn.info)){
      stop("Error: 'arg' should be between 1 and ", ncol(cn.info))
    }
  }else{
    stop("Error: 'arg' should be string or numeric.")
  }

  cn.info[[chr.col]] <- as.character(cn.info[[chr.col]])

  if(!is.character(cn.info[[chr.col]])){
    stop("Error: chromosome information must be string or numeric.")
  }else if(!(is.character(cn.info[[tcn.col]]) | is.numeric(cn.info[[tcn.col]]))){
    stop("Error: total copy number must be string or numeric.")
  }else if(!is.numeric(cn.info[[start.col]])){
    stop("Error: start position must be numeric.")
  }else if(!is.numeric(cn.info[[end.col]])){
    stop("Error: end position must be numeric.")
  }

  return(list(cn.info = cn.info, chr.col = chr.col, start.col = start.col, end.col = end.col, A.col = A.col, B.col = B.col, tcn.col = tcn.col, estimate.alleles = estimate.alleles))
}

.estimate_alleles <- function(cn.info, tcn.col){
  cn.info$A <- ceiling(cn.info[[tcn.col]]/2)
  cn.info$B <- floor(cn.info[[tcn.col]]/2)

  return(cn.info)
}
