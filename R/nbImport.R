#' Combine CNVs and SNVs
#' @description
#' Merges CNVs and SNVs into a single data.table. Each variant is assigned to its corresponding copy number segment and status.
#' @param cnv CNV data from \code{\link{readCNV}}
#' @param snv SNV data from \code{\link{readVCF}}
#' @param purity tumor cell content
#' @param ploidy average copy number in the tumor sample
#' @param sig.assign Logical. If TRUE, each variant will be assigned to the most likely mutational signature
#' @param ID sample name
#' @param sig.file File path to the SigAssignment output file, typically named "Decomposed_MutationType_Probabilities.txt".
#' @param sig.select A character vector of specific signatures to include in the analysis (e.g., c("SBS1", "SBS5", "SBS40") to focus on clock-like mutational processes).
#' @param min.p Numeric. The minimum probability threshold from the SigAssignment output that a variant must meet to be considered as matching a specific signature.
#' @examples
#' snvs <- system.file("extdata", "NBE15", "snvs_NBE15_somatic_snvs_conf_8_to_10.vcf", package = "LACHESIS")
#' s_data <- readVCF(vcf = snvs, vcf.source = "dkfz")
#' aceseq_cn <- system.file("extdata", "NBE15", "NBE15_comb_pro_extra2.51_1.txt", package = "LACHESIS")
#' c_data <- readCNV(aceseq_cn)
#' nb <- nbImport(cnv = c_data, snv = s_data, purity = 1, ploidy = 2.51)
#' @seealso \code{\link{plotNB}}
#' @return a data.table
#' @importFrom RColorBrewer Set3
#' @export

nbImport <- function(cnv = NULL, snv = NULL, purity = NULL, ploidy = NULL, sig.assign = FALSE, ID = NULL, sig.file = NULL, sig.select = NULL, min.p = NULL){

  end <- start <- NULL

  if(any(is.null(cnv), is.null(snv))){
    stop("Missing snv and cnv inputs!")
  }
  if(any(is.null(purity), is.null(ploidy))){
    stop("Missing purity and ploidy inputs!")
  }

  colnames(cnv)[1:3] <- c("chrom", "start", "end")
  data.table::setDT(x = cnv, key = c("chrom", "start", "end"))
  colnames(snv)[1:2] <- c("chrom", "start")
  snv[,end := start]
  data.table::setDT(x = snv, key = c("chrom", "start", "end"))

  sv <- data.table::foverlaps(x = snv, y = cnv, type = "within")

  if(nrow(sv) == 0){
    stop("No overlapping SNVs found within CNV regions")
  }

  if(nrow(sv[is.na(start)])){
    warning("Removed ", nrow(sv[is.na(start)]), " variants with no copy number overlaps")
    sv <- sv[!is.na(start)]
  }

  if(sig.assign == TRUE){
    t.sample <- attributes(sv)$t.sample
    assign.result <- .assign_signatures(sv, sig.file, ID, sig.select, min.p)
    sv <- assign.result$sv
    sig.colors <- assign.result$sig.colors
    attr(sv, "t.sample") <- t.sample
  }

  #Make columns more intuitive
  colnames(sv)[which(colnames(sv) == "i.start")] <- "snv_start"
  colnames(sv)[which(colnames(sv) == "i.end")] <- "snv_end"
  colnames(sv)[which(colnames(sv) == "start")] <- "cn_start"
  colnames(sv)[which(colnames(sv) == "end")] <- "cn_end"
  attr(sv, "cnv") <- cnv
  attr(sv, "purity") <- as.numeric(purity)
  attr(sv, "ploidy") <- as.numeric(ploidy)
  attr(sv, "sig.colors") <- sig.colors
  sv
}

.assign_signatures <- function(sv = NULL, sig.file = NULL, ID = NULL, sig.select = NULL, min.p = NULL) {

  if (is.null(sv)) {
    stop("Missing 'sv' input data!")
  }

  if (is.null(sig.file)) {
    stop("Missing 'SigAssignment' input data!")
  }

  sig.data <- fread(sig.file)

  setnames(sig.data, c("Sample Names", "Chr", "Pos"), c("Sample", "chrom", "i.start"))

  sbs.cols <- grep("^SBS", names(sig.data), value = TRUE)

  sig.data <- sig.data[, {
    max.p.sig <- which.max(.SD)
    list(
      Signature = sbs.cols[max.p.sig],
      Probability = .SD[[max.p.sig]]
    )
  }, by = .(Sample, chrom, i.start, MutationType), .SDcols = sbs.cols]

  if (!is.null(min.p)) {
    sig.data <- sig.data[Probability >= min.p]
  }

  sv[, Sample := ID]

  sv <- merge(sv, sig.data, by = c("Sample", "chrom", "i.start"), all.x = TRUE)

  if (!is.null(sig.select)) {
    sv <- sv[Signature %in% sig.select]
    sig.number <- length(sig.select)
    sig.colors <- setNames(.get_sig_colors(sig.number), sig.select)
  } else {
    sig.options <- unique(sig.data$Signature)
    sig.number <- length(sig.options)
    sig.colors <- setNames(.get_sig_colors(sig.number), sig.options)
  }

  sv[, "Sample" := NULL]

  return(list(sv = sv, sig.colors = sig.colors))
}

.get_sig_colors <- function(n, palette = "Set3", max.colors = 12) {
  base.colors <- brewer.pal(min(max.colors, n), palette)
  if (n > max.colors) {
    colorRampPalette(base.colors)(n)
  } else {
    base.colors
  }
}

#' Plot VAF distribution per copy number
#' @description
#' Visualizes results from  \code{\link{nbImport}}. Top plot, measured copy numbers along the genome; bottom plots, VAF histograms of SNVs stratified by copy number and minor/major allele count.
#' @param nb output generated from \code{\link{nbImport}}
#' @param ref.build Reference genome. Default `hg19`. Can be `hg18`, `hg19` or `hg38`
#' @param min.cn maximum copy number to be included in the plotting. Defaults to 2.
#' @param max.cn maximum copy number to be included in the plotting. Defaults to 4.
#' @param nb.col.abline optional, the color code for the abline.
#' @param nb.col.cn.2 optional, the color code if tcn = 2.
#' @param nb.col.cn optional, the color code if other copy numbers.
#' @param nb.col.hist optional, the color code for histograms.
#' @param nb.border, optional, the line color.
#' @param nb.breaks optional; the number of bins in the histogram.
#' @param samp.name Sample name. Optional. Default NULL
#' @param output.file optional, will save the plot.
#' @param sig.show plot stratified VAF histogram with assigned mutational signatures.
#' @param sig.output.file optional, will save the stratified VAF histogram with mutational signatures.
#' @param ... further arguments and parameters passed to other LACHESIS functions.
#' @examples
#' snvs = system.file("extdata", "NBE15", "snvs_NBE15_somatic_snvs_conf_8_to_10.vcf", package = "LACHESIS")
#' s_data <- readVCF(vcf = snvs, vcf.source = "dkfz")
#' aceseq_cn <- system.file("extdata", "NBE15", "NBE15_comb_pro_extra2.51_1.txt", package = "LACHESIS")
#' c_data <- readCNV(aceseq_cn)
#' nb <- nbImport(cnv = c_data, snv = s_data, purity = 1, ploidy = 2.51)
#' plotNB(nb)
#' @export
#' @importFrom graphics abline axis box grid hist mtext par rect text title
#' @import ggplot2

plotNB <- function(nb = NULL, ref.build = "hg19", min.cn = 2, max.cn = 4, nb.col.abline = "gray70", nb.col.cn.2 = "#7f8c8d", nb.col.cn = "#16a085", nb.col.hist = "#34495e", nb.border = NA, nb.breaks = 100, samp.name = NULL, output.file = NULL, sig.show = FALSE, sig.output.file = NULL, ...){

  chrom <- start <- t_vaf <- NULL

  if(is.null(nb)){
    stop("Missing input. Please provide the output generated by nbImport")
  }

  if(max.cn <= min.cn){
    stop("max.cn must be larger than min.cn")
  }

  if(!is.null(output.file)){
    pdf(output.file, width = 7, height = 9)
  }

  sig.colors <- attr(nb, "sig.colors")
  purity = attr(nb, "purity")

  segs <- attr(nb, "cnv")
  segs <- segs[order(chrom, start)]
  colnames(segs)[1:3] <- c("Chromosome", "Start_Position", "End_Position")
  segs <- .transformSegments(segmentedData = segs, build = ref.build)

  contig_lens <- cumsum(.getContigLens(build = ref.build))

  #n_copies <- length(min.cn:max.cn)
  n_copy_combs <- nrow(unique(nb[TCN >= min.cn & TCN <= max.cn,TCN, B]))
  n_copy_combs <- n_copy_combs + n_copy_combs%%2
  lo_mat <- matrix(data = c(rep(1, n_copy_combs/2), 2:(n_copy_combs/2+1), (n_copy_combs/2 + 2):(n_copy_combs+1)), nrow = 3, byrow = TRUE)
  # lo_mat <- matrix(data = c(rep(1, n_copies), 2:(n_copies+1)), (n_copies + 2), nrow = 3, byrow = TRUE)
  graphics::layout(mat = lo_mat, heights = c(3, 2, 2))
  par(mar = c(3, 4, 4, 3))
  plot(NA, ylim = c(0, max.cn), xlim = c(0, max(contig_lens)), axes = FALSE, xlab = NA, ylab = NA)
  abline(h = 1:max.cn, v = contig_lens, lty = 2, col = nb.col.abline, lwd = 0.4)
  rect(xleft = segs$Start_Position_updated, ybottom = segs$TCN-0.1, xright = segs$End_Position_updated, ytop = segs$TCN+0.1, col = ifelse(segs$TCN == 2, nb.col.cn.2, nb.col.cn), border = nb.border, lty = 3)
  contig_lens_mid <- c(contig_lens[1] / 2,
                       (contig_lens[-length(contig_lens)] + contig_lens[-1]) / 2)
  axis(side = 1, at = c(0, contig_lens), labels = FALSE, pos = 0)
  axis(side = 1, at = contig_lens_mid, labels = c(1:22, "X", "Y"), tick = FALSE, line = -0.5, cex.axis = 0.9)
  axis(side = 2, at = 0:max.cn, labels = 0:max.cn, las = 2)
  mtext(text = "Total CN", side = 2, line = 2, cex = 0.9)
  mtext(text = "Chromosome", side = 1, line = 2, cex = 0.9)
  title(main = ifelse(is.null(samp.name), yes = attr(nb, "t.sample"), no = samp.name))

  nb <- nb[TCN >= min.cn & TCN <= max.cn,]
  nb$TCN <- factor(nb$TCN, levels = 1:max.cn)
  nb <- split(nb, nb$TCN)

  for(cn in seq_along(nb)){
    nb. <- split(nb[[cn]], nb[[cn]]$B)
    ploidy = names(nb)[cn]
    for(b in seq_along(nb.)){
      tcn <- nb.[[b]]
      B <- names(nb.)[b]
      if(nrow(tcn) == 0){
        par(mar = c(3, 4, 3, 1))
        plot(NA, xlim = c(0, 1), ylim = c(0, 50), xlab = NA, ylab = NA, col = nb.col.hist, main = NA, frame.plot = FALSE)
        title(main = paste0("Total CN:", ploidy), cex.main = 1.2)
        mtext(text = "No. of SNVs", side = 2, line = 2.5, cex = 0.7)
        mtext(text = "VAF", side = 1, line = 1.8, cex = 0.7)
      }else{
        par(mar = c(3, 4, 3, 1))
        hist(tcn[,t_vaf], breaks = nb.breaks, xlim = c(0, 1), xlab = NA, ylab = NA, border = nb.border, col = nb.col.hist, main = NA)
        title(main = paste0("CN:", as.numeric(ploidy), " (", as.numeric(ploidy) - as.numeric(B), ":", as.numeric(B), ")"), cex.main = 1.2)
        mtext(text = "No. of SNVs", side = 2, line = 2.5, cex = 0.7)
        mtext(text = "VAF", side = 1, line = 1.8, cex = 0.7)
        if(!is.null(purity)){
          abline(v = .expectedClVAF(CN = as.numeric(names(nb)[cn]), purity = purity), lty = 2)
        }
      }
    }
  }

  if(!is.null(output.file)){
    dev.off()
  }

  if (sig.show == TRUE) {
    pdf(file = sig.output.file, width = 8, height = 6)

    for (cn in seq_along(nb)) {
      nb. <- split(nb[[cn]], nb[[cn]]$B)
      ploidy <- names(nb)[cn]

      for (b in seq_along(nb.)) {
        tcn <- nb.[[b]]
        B <- names(nb.)[b]

        if (nrow(tcn) > 0) {
          p <- ggplot(tcn, aes(x = t_vaf, fill = Signature)) +
            geom_histogram(bins = 100, color = "black", position = "stack") +
            scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
            labs(x = "VAF", y = "No. of SNVs") +
            ggtitle(paste0("CN:", ploidy, " (", as.numeric(ploidy) - as.numeric(B), ":", B, ")")) +
            scale_fill_manual(values = sig.colors) +
            theme_minimal()
          if (!is.null(purity)) {
            expected_vafs <- .expectedClVAF(CN = as.numeric(names(nb)[cn]), purity = purity)
            p <- p + geom_vline(xintercept = expected_vafs, linetype = "dashed")
          }

          print(p)
        }
      }
    }

    dev.off()
  }
}


#contig lengths for hg19, hg38 and hg18
.getContigLens <- function(build = "hg19"){

  if(build == 'hg19'){
    chr.lens = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663,
                 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540,
                 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566,
                 155270560, 59373566)
  } else if(build == 'hg18'){
    chr.lens = c(247249719, 242951149, 199501827, 191273063, 180857866, 170899992,
                 158821424, 146274826, 140273252, 135374737, 134452384, 132349534,
                 114142980, 106368585, 100338915, 88827254, 78774742, 76117153,
                 63811651, 62435964, 46944323, 49691432, 154913754, 57772954)
  } else if(build == 'hg38'){ #hg38
    chr.lens = c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
                 159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
                 114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
                 58617616, 64444167, 46709983, 50818468, 156040895, 57227415)
  } else{
    stop('Available reference builds: hg18, hg19, hg38')
  }

  chr.lens

}

#--- Change segment sizes into linear scale
.transformSegments <- function(segmentedData, build = 'hg19'){

  Start_Position <- End_Position <- Chromosome <- NULL

  build.opts <- c('hg19', 'hg18', 'hg38')

  if(!build %in% build.opts){
    stop('Available reference builds: hg18, hg19, hg38')
  }

  #Get chr lens
  chr.lens <- .getContigLens(build = build)

  segmentedData[,Start_Position := as.numeric(as.character(Start_Position))]
  segmentedData[,End_Position := as.numeric(as.character(End_Position))]

  #Replace chr x and y with numeric value (23 and 24) for better ordering
  segmentedData$Chromosome <- gsub(pattern = 'chr', replacement = '', x = segmentedData$Chromosome, fixed = TRUE)
  segmentedData$Chromosome <- gsub(pattern = 'X', replacement = '23', x = segmentedData$Chromosome, fixed = TRUE)
  segmentedData$Chromosome <- gsub(pattern = 'Y', replacement = '24', x = segmentedData$Chromosome, fixed = TRUE)

  segmentedData$Chromosome <- factor(x = segmentedData$Chromosome, levels = 1:24, labels = 1:24)

  segmentedData <- segmentedData[order(Chromosome, Start_Position, decreasing = FALSE)]

  seg.spl <- split(segmentedData, segmentedData$Chromosome)

  seg.spl.transformed <- seg.spl[[1]]
  if(nrow(seg.spl.transformed) > 0){
    seg.spl.transformed$Start_Position_updated <- seg.spl.transformed$Start_Position
    seg.spl.transformed$End_Position_updated <- seg.spl.transformed$End_Position
  }

  chr.lens.sumsum <- cumsum(chr.lens)

  for(i in 2:length(seg.spl)){

    x.seg <- seg.spl[[i]]
    if(nrow(x.seg) > 0){
      x.seg$Start_Position_updated <- x.seg$Start_Position + chr.lens.sumsum[i-1]
      x.seg$End_Position_updated <- x.seg$End_Position + chr.lens.sumsum[i-1]
    }
    seg.spl.transformed <- rbind(seg.spl.transformed, x.seg, fill = TRUE)
  }

  return(seg.spl.transformed)
}

# expected clonal VAFs for copy number CN at a given purity on autosomes
.expectedClVAF <- function(CN, purity){
  (1:CN)*purity/(purity*CN + 2*(1-purity))
}

