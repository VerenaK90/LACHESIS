#' Plot histogram of VAF distribution
#' @description
#' Plot frequency distribution of variant allele frequencies
#'
#' @param vaf output produced by \code{\link{readVCF}}
#' @param vafbreak Interval size. Default 0.05
#' @param t_sample Sample name for tumor. Used for plot title. Default NULL
#' @examples
#' strelka_vcf = system.file("extdata", "strelka2.somatic.snvs.vcf.gz", package = "NBevolution")
#' s_data = readVCF(vcf = strelka_vcf, vcf_source = "strelka", ignore_XY = FALSE)
#' plotVAFdistr(s_data)
#' @return data.table of frequency table
#' @export
plotVAFdistr <- function(vaf = NULL, vafbreak = 0.05, t_sample = NULL){

  if(is.null(vaf)){
    stop("Missing input vaf. Use `readVCF` to extract vaf.")
  }

  par(mar = c(4, 4, 2, 1))
  h <-hist(vaf$t_vaf, breaks = seq(0, 1, vafbreak), col = "#34495e", border = "#bdc3c7", xlab = NA, ylab = NA, axes = FALSE, main = NA)
  axis(side = 1, at = seq(0, 1, 0.1))
  axis(side = 2, at = pretty(h$counts), las = 2)
  text(x = h$mids, y = h$counts, labels = h$counts, pos = 3, xpd = TRUE)

  if(is.null(x = t_sample)){
    t_sample = attr(vaf, which = 't_sample')
  }
  title(main = t_sample)

  mtext(side = 1, text = paste0("VAF (interval: ", vafbreak, ")"), line = 2)
  mtext(side = 2, text = "Frequency", line = 2.5)

  data.table::data.table(breaks = h$breaks, counts = c(0, h$counts))

}
