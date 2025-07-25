#' Plot normalized mutation density at copy number gain and MRCA per segment
#' @description
#' Visualizes results from \code{\link{MRCA}}. Top plot, histograms of mean mutation densities; bottom plots, timeline of early tumor evolution, showing mutation densities (mean and 95% CI) of individual chromosomal gains and mutation densities at ECA and MRCA.
#' @param mrcaObj output generated from \code{\link{MRCA}}
#' @param samp.name sample name, optional
#' @param min.seg.size minimal segment size to plot
#' @param mut.col.zero optional, the bar color for densities of mutations present on single copies.
#' @param mut.col.multi  optional, the bar color for densities of mutations present on multiple copies.
#' @param mut.border optional, the line color
#' @param mut.show.density optional; if `TRUE`, the density distribution of mutation densities on single copies will be shown in the histogram of mutation densities on multiple copies.
#' @param mut.breaks optional; the number of bins in the histogram.
#' @param mut.show.realtime logical; if `TRUE`, displays weeks post-conception on the evolutionary timeline.
#' @param mut.snv.rate optional; rate of accumulated SNVs per day in a diploid genome (i.e. 3.2 SNVs/day in neuroblastoma)
#' @param output.file optional; will save the plot.
#' @param ... further arguments and parameters passed to other LACHESIS functions.
#' @examples
#' snvs <- system.file("extdata", "NBE15", "snvs_NBE15_somatic_snvs_conf_8_to_10.vcf", package = "LACHESIS")
#' s_data <- readVCF(vcf = snvs, vcf.source = "dkfz")
#' aceseq_cn <- system.file("extdata", "NBE15", "NBE15_comb_pro_extra2.51_1.txt", package = "LACHESIS")
#' c_data <- readCNV(aceseq_cn)
#' nb <- nbImport(cnv = c_data, snv = s_data, purity = 1, ploidy = 2.51)
#' cl_muts <- clonalMutationCounter(nb)
#' norm_muts <- normalizeCounts(cl_muts)
#' mrca <- MRCA(norm_muts)
#' plotMutationDensities(mrca)
#' @export
#' @importFrom graphics abline Axis box grid hist mtext par rect text title arrows legend points polygon

plotMutationDensities <- function(mrcaObj = NULL, samp.name = NULL, min.seg.size = 10^7, mut.col.zero = "#4FB12B", mut.col.multi = "#176A02", mut.border = NULL, mut.show.density = TRUE, mut.breaks = NULL, mut.show.realtime = FALSE, mut.snv.rate = 3.2, output.file = NULL, ...){

  Seglength <- . <- A <- B <- variable <- value <- lines <- density <- chrom <- TCN <- Seglength <- n_mut_A <- n_mut_B <- n_mut_total <- density_total_mean <- density_A_mean <- density_B_mean <- density_total_lower <- density_total_upper <- density_A_lower <- density_A_upper <- density_B_lower <- density_B_upper <- p_total_to_mrca <- p_A_to_mrca <- p_B_to_mrca <- p_adj_total_to_mrca <- p_adj_A_to_mrca <- p_adj_B_to_mrca <- MRCA_qual <- p_A_to_eca <- p_B_to_eca <- p_adj_A_to_eca <- p_adj_B_to_eca <- A_time <- B_time <- NULL
  if(is.null(mrcaObj)){
    stop("Missing input. Please provide the output generated by MRCA()")
  }
  if(!is.null(output.file)){
    pdf(output.file, width = 7, height = 6)
  }

  to.plot <- data.table::melt(mrcaObj, id.vars = c("chrom", "TCN", "A", "B", "Seglength"),
                              measure.vars = c("density_total_mean", "density_A_mean", "density_B_mean"))

  # subset on segments larger min.seg.size; for plotting CNVs restrict to non-normal copy numbers
  to.plot <- to.plot[Seglength > min.seg.size &
                       !(variable == "density_total_A" & A == 1) &
                       !(variable == "density_total_B" & B == 1),]

  if(is.null(mut.breaks)){
    mut.breaks = 20
  }

  lo_mat <- matrix(data = c(1, 2, 3, 3), nrow = 2, ncol=2, byrow = TRUE)
  graphics::layout(mat = lo_mat, widths = c(1, 1, 2), heights = c(1, 1, 1))

  par(mar = c(3, 4, 3, 1))

  hist(to.plot[variable == "density_total_mean",value], xlim = c(0, 1.05 * max(to.plot[,value])),
       breaks = mut.breaks, col = mut.col.zero, border = mut.border, main = NA,
       xlab = NA, ylab = NA)

  title(main = paste("Single-copy SNV densities"), cex.main = 1.2)
  mtext(text = "SNVs per Mb", side = 1, line = 2.5, cex = 0.8)
  mtext(text = "No. of genomic segments", side = 2, line = 1.8, cex = 0.8)

  if(nrow(to.plot[(variable == "density_A_mean" & A > 1) |
                  (variable == "density_B_mean" & B > 1 & A != B),]) > 0){
    #Get y and x axis limits
    temp_d = density(to.plot[variable == "density_total_mean",value])
    temp_hist <- hist(to.plot[(variable == "density_A_mean" & A > 1) |
                                (variable == "density_B_mean" & B > 1 & A != B), value],
                      breaks = mut.breaks, plot = FALSE)
    temp_hist_y <- max(temp_hist$counts, na.rm = TRUE)
    temp_hist_x <- max(temp_hist$breaks, na.rm = TRUE)

    par(mar = c(3, 4, 3, 1))
    hist(to.plot[(variable == "density_A_mean" & A > 1) |
                   (variable == "density_B_mean" & B > 1 & A != B),value],
         xlim = c(0, max(c(temp_hist_x, max(temp_d$x, na.rm = TRUE)))),
         ylim = c(0, max(c(temp_hist_y, max(temp_d$y, na.rm = TRUE)))),
         breaks = mut.breaks, col = mut.col.multi, border = mut.border, main = NA,
         xlab = NA, ylab = NA)

    # add density of MRCA
    if(mut.show.density == TRUE & nrow(to.plot[variable == "density_total_mean"])>1){
      lines(density(to.plot[variable == "density_total_mean",value]), lty = 2)
    }

    title(main = paste("Multi-copy SNV densities"), cex.main = 1.2)
    mtext(text = "SNVs per Mb", side = 1, line = 2.5, cex = 0.8)
    mtext(text = "No. of genomic segments", side = 2, line = 1.8, cex = 0.8)

  }else{
    plot(NA, NA, xlim = c(0, 1), ylim = c(0, 1), xlab = NA, ylab = NA, main = NA, axes = FALSE, frame.plot = FALSE)
  }


  # Timeline summary
  if(mut.show.realtime){
    par(mar = c(3, 1, 5, 1), xpd = FALSE)
  }
  else{
    par(mar = c(3, 1, 3, 1), xpd = FALSE)
  }
  x.min = 0
  x.max = max(c(mrcaObj$density_total_upper, mrcaObj$density_A_upper, mrcaObj$density_B_upper), na.rm = TRUE)*1.3
  y.min = 0
  y.max.a = nrow(mrcaObj[A>1,])
  y.max = max(c(1, nrow(mrcaObj[A>1,]) + nrow(mrcaObj[B>1 & B!=A])))
  plot(NA, NA, xlim = c(x.min, x.max), ylim = c(y.min, y.max), xlab = NA, ylab = NA, main = NA, axes = FALSE, frame.plot = FALSE)
  Axis(side = 1, cex = 0.7)
  mtext(text = "SNVs per Mb", side = 1, line = 2, cex = 0.7)

  if(mut.show.realtime){
    weeks_pc <- c(12, 27, 38, 64, 90, 116)
    snvs_per_mb <- (weeks_pc - 2) * 7 * mut.snv.rate / (3300 * 2) # Converting SNVs per day to SNVs per Mb starting from Gastrulation (-2 weeks), assuming haploid genome of 3300Mb
    realtime_labels <- c("12w", "27w", "38w", "6m", "12m", "18m")
    axis(side = 3, at = c(x.min, snvs_per_mb, x.max), labels = c("", realtime_labels, ""), cex.axis = 0.7)
    segments(x0 = x.min, y0 = par("usr")[4], x1 = x.max, y1 = par("usr")[4], xpd = NA)
    mtext("Estimated time (weeks post conception and months postnatal)", side = 3, line = 2, cex = 0.7)

    title(main = paste("Evolutionary timeline of chromosomal gains and losses"), cex.main = 1.2, line = 3.2)
  }
  else{
    title(main = paste("Evolutionary timeline of chromosomal gains and losses"), cex.main = 1.2)
  }

  # ECA:
  polygon(c(attr(mrcaObj, "ECA_time_lower"), rep(attr(mrcaObj, "ECA_time_upper"),2),attr(mrcaObj, "ECA_time_lower")), c(rep(y.min, 2), rep(y.max,2)),
          col = mut.col.multi, border = NA)
  abline(v = attr(mrcaObj, "ECA_time_mean"), lty = 2)

  #MRCA:
  polygon(c(attr(mrcaObj, "MRCA_time_lower"), rep(attr(mrcaObj, "MRCA_time_upper"),2),attr(mrcaObj, "MRCA_time_lower")), c(rep(y.min, 2), rep(y.max,2)),
          col = mut.col.zero, border = NA)
  abline(v = attr(mrcaObj, "MRCA_time_mean"), lty = 2)

  signs <- c("ECA" = 19, "MRCA" = 17, "ECA/MRCA" = 15, "not mapped to ECA or MRCA" = 1)
  # A alleles:
  if(nrow(mrcaObj[A>1,])>0){
    points(mrcaObj[A>1,density_A_mean], 1:mrcaObj[,sum(A>1)], col=1:mrcaObj[,sum(A>1)], pch=signs[mrcaObj[A>1,A_time]])
    arrows(x0=mrcaObj[A>1,density_A_lower], y0=1:mrcaObj[,sum(A>1)], x1=mrcaObj[A>1,density_A_upper], y1=1:mrcaObj[,sum(A>1)], code=3, angle=90, length=0, col=1:mrcaObj[,sum(A>1)], lwd=1)
    legend("topright",box.lwd = 0, pch=signs[names(signs) %in% mrcaObj$A_time | names(signs) %in% mrcaObj$B_time], legend = names(signs[names(signs) %in% mrcaObj$A_time | names(signs) %in% mrcaObj$B_time]), cex = 0.7)
  }
  # B alleles:
  if(nrow(mrcaObj[B>1 & B!=A,])>0){
    points(mrcaObj[B>1 & B!=A,density_B_mean], (y.max.a+1):(y.max.a+mrcaObj[,sum(B>1 & B!=A)]), col=(y.max.a+1):(y.max.a+mrcaObj[,sum(B>1 & B!=A)]), pch=signs[mrcaObj[B>1 & B!=A,B_time]])
    arrows(x0=mrcaObj[B>1 & B!=A,density_B_lower], y0=(y.max.a+1):(y.max.a+mrcaObj[,sum(B>1 & B!=A)]), x1=mrcaObj[B>1 & B!=A,density_B_upper], y1=(y.max.a+1):(y.max.a+mrcaObj[,sum(B>1 & B!=A)]), code=3, angle=90, length=0, col=(y.max.a+1):(y.max.a+mrcaObj[,sum(B>1 & B!=A)]), lwd=1)
    legend("bottomright", box.lwd = 0, lty=1, col= c(1:mrcaObj[,sum(A>1)],  (y.max.a+1):(y.max.a+mrcaObj[,sum(B>1 & B!=A)])),
           legend = c(paste(paste0("chr", mrcaObj[A>1,chrom]), mrcaObj[A>1,TCN], mrcaObj[A>1,A], sep = "_"),
                      paste(paste0("chr", mrcaObj[B>1 & B!=A,chrom]), mrcaObj[B>1 & B!=A,TCN], mrcaObj[B>1 & B!=A,B], sep = "_")
           ), cex = 0.7, ncol = 2)
  }else{
    legend("bottomright", box.lwd = 0, lty=1, col= c(1:mrcaObj[,sum(A>1)]),
           legend = paste(paste0("chr", mrcaObj[A>1,chrom]), mrcaObj[A>1,TCN], mrcaObj[A>1,A], sep = "_"), cex = 0.7, ncol = 2)

  }

  if(!is.null(output.file)){
    dev.off()
  }
}
