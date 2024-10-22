#' Uses a file with log2 ratio copy number to calculate the total copy number (tcn)
#'
#' @author {Maximilia Eggle}
#' @description #'Uses a file with log2 ratio copy number to calculate the total copy number (tcn)
#'
#' @param cn.info Path to file with log2 ratio copy number
#' @param logR_col Name of column with log value


.estimate_tcn_from_log2ratios <- function(cn.info, logR.col = NULL) {

  ## Reading as data table
  cn.info <- data.table::fread(cn.info)

  ## Identifying column with copy number information
  if (is.null(logR.col)) {
    logR.col <- colnames(cn.info)[
      grepl("\\btcn\\b", colnames(cn.info), ignore.case = TRUE) |
        grepl("\\bcnt\\b", colnames(cn.info), ignore.case = TRUE) |
        grepl("\\bcopynumber\\b", colnames(cn.info), ignore.case = TRUE) |
        grepl("\\bcopy number\\b", colnames(cn.info), ignore.case = TRUE) |
        grepl("\\bseg.mean\\b", colnames(cn.info), ignore.case = TRUE) |
        grepl("\\bseg.cn\\b", colnames(cn.info), ignore.case = TRUE) |
        grepl("\\bcn\\b", colnames(cn.info), ignore.case = TRUE) |
        grepl("\\blog2.copy.number\\b", colnames(cn.info), ignore.case = TRUE)
    ]
  }

  if (length(logR.col) == 0) {
    stop("Error: TCN identifier could not be inferred, please provide logR.col.")
  }

  ## Calculating tumor cell content
  cn.info[, tcn := 2 * (2 ^ get(logR.col))]
  cn.info[, tcn := round(tcn)]

  ## Writing data in a txt-file
  data.table::fwrite(cn.info, file = "converted_log_cnv.txt", sep = "\t")

  return(cn.info)
}


