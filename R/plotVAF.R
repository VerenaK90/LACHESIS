#' Plot histogram of VAF distribution
#' @description
#' Plot frequency distribution of variant allele frequencies
#'
#' @param vaf output produced by \code{\link{readVCF}}
#' @param vaf.interval Interval size. Default 0.05
#' @param t_sample Sample name for tumor. Used for plot title. Default NULL
#' @param vaf.show.counts Show counter per break on the histogram. Default FALSE
#' @param vaf.show.density Show additional inset plot of density. Default TRUE
#' @param vaf.col Color to be used to fill the bars, default "#34495e"
#' @param vaf.border Border color, default "#bdc3c7"
#' @param srtcounts Text angle if `vaf.show.counts` is TRUE. Default 45
#' @param output.file Optional, will save the plot.
#' @param ... further arguments and parameters passed to other
#' LACHESIS functions.
#' @return VAF histogram
#' @examples
#' strelka_vcf <- system.file("extdata", "strelka2.somatic.snvs.vcf.gz",
#'     package = "LACHESIS"
#' )
#' s_data <- readVCF(
#'     vcf = strelka_vcf, vcf.source = "strelka",
#'     ignore.XY = FALSE
#' )
#' plotVAFdistr(s_data)
#' @importFrom graphics axis box grid hist mtext par text title
#' @export
plotVAFdistr <- function(vaf = NULL, vaf.interval = 0.05, t_sample = NULL,
                         vaf.show.counts = FALSE, vaf.show.density = TRUE,
                         vaf.col = "#34495e", vaf.border = "#bdc3c7",
                         srtcounts = 45, output.file = NULL, ...) {
    if (is.null(vaf)) {
        stop("Missing input vaf. Use `readVCF` to extract vaf.")
    }

    if (!is.null(output.file)) {
        pdf(output.file, width = 4, height = 4)
    }

    par(mar = c(4, 4, 2, 1), fig = c(0, 1, 0, 1))
    h <- hist(vaf$t_vaf,
        breaks = seq(0, 1, vaf.interval), col = vaf.col,
        border = vaf.border, xlab = NA, ylab = NA, axes = FALSE, main = NA
    )
    axis(side = 1, at = seq(0, 1, 0.1))
    axis(side = 2, at = pretty(h$counts), las = 2)
    if (vaf.show.counts) {
        text(
            x = h$mids, y = h$counts, labels = h$counts, pos = 3, xpd = TRUE,
            srt = srtcounts, cex = 0.7
        )
    }


    if (is.null(x = t_sample)) {
        t_sample <- attr(vaf, which = "t.sample")
    }
    title(main = t_sample)

    mtext(
        side = 1, text = paste0("VAF (interval: ", vaf.interval, ")"),
        line = 2
    )
    mtext(side = 2, text = "Frequency", line = 2.5)

    if (vaf.show.density) {
        par(fig = c(0.65, 1, 0.55, 0.95), new = TRUE, mar = c(3, 2, 1, 1))
        plot(density(vaf$t_vaf),
            frame.plot = FALSE, main = NA, axes = FALSE,
            xlim = c(0, 1), col = vaf.col
        )
        grid(nx = 5)
        box(lty = 1, lwd = 0.1)
        mtext(text = "Density", side = 2, cex = 0.7)
        axis(
            side = 1, at = seq(0, 1, 0.2), cex.axis = 0.65, line = -1,
            tick = FALSE
        )
    }

    data.table::data.table(breaks = h$breaks, counts = c(0, h$counts))

    if (!is.null(output.file)) {
        dev.off()
    }
}
