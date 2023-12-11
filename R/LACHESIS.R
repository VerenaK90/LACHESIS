#' Run MRCA density estimation for a set of tumors
#' @description
#' Takes a set of SNV and CNV files as input and outputs per-tumor SNV densities. Input can either be a tab-delimited file containing
#' the sample specifications or vectors giving direct paths to the sample files. CNV file requires columns for the chromosome number, start and end of the segment, and either the total copy number or the number of A- and B-alleles
#' @param input.files a tab-delimited sample-specification file, it must contain the sample name, the path to the SNV file, path to CNV file, and optionally purity, ploidy, cnv.chr.col, cnv.start.col, cnv.end.col, cnv.A.col, cnv.B.col, cnv.tcn.col. A template for this spreadsheet can be downloaded from ...
#' @param ids vector of sample names, will be ignored if `input.files` is specified.
#' @param cnv.files vector of cnv files in same order as ids; should be in tab-delimited format, will be ignored if `input.files` is specified.
#' @param snv.files vector of snv files in same order as ids; should be in vcf format, will be ignored if `input.files` is specified.
#' @param vcf.source Tool used for generating VCF file. Can be `strelka` or `mutect` or `dkfz`
#' @param purity vector tumor cell content in same order as ids; will be ignored if `input.files` is specified.
#' @param ploidy average copy number in the tumor sample in same order as ids; will be ignored if `input.files` is specified.
#' @param cnv.chr.col column index of chromosome number in cnv.files
#' @param cnv.start.col column index of first position of the segment
#' @param cnv.end.col column index of last position of the segment
#' @param cnv.A.col column index of the number of A alleles. If A and B are not provided, allele configuration are assumed as 1:1 for disomic, 2:1 for trisomic and 3:1 for tetrasomic regions.
#' @param cnv.B.col column index of the number of B alleles. If A and B are not provided, allele configuration are assumed as 1:1 for disomic, 2:1 for trisomic and 3:1 for tetrasomic regions.
#' @param cnv.tcn.col column index of the total copy number. Is computed to A + B if not provided.
#' @param output.dir link to directory in which output is to be stored.
#' @param age, optional, the age at diagnosis
#' @param OS.time, optional, overall survival time
#' @param OS, optional, overall survival indicator variable
#' @param EFS.time, optional, event-free survival time
#' @param EFS, optional, event-free survival indicator variable
#' optional arguments
#' @param min.vaf Remove variants with vcf below threshold. Default 0.01
#' @param min.depth Minimum required depth for a variant to be considered. Default 30.
#' @param min.cn maximum copy number to be included in the plotting. Defaults to 2.
#' @param max.cn maximum copy number to be included in the analysis. Defaults to 4.
#' @param chromosomes vector of chromosomes to include in SNV density estimation
#' @param merge.tolerance the maximum distance below which adjacent segments with equal copy number are merged. Defaults to 10^5 bp.
#' @param ignore.XY Ignore allosomes. Default TRUE
#' @param ref_build Reference genome. Default `hg19`. Can be `hg18`, `hg19` or `hg38`
#' @examples
#' input.files = system.file("extdata", "Sample_template.txt")
#' lachesis <- LACHESIS(input.files = input.files)
#' @seealso \code{\link{MRCA}} \code{\link{clonalMutationCounter}} \code{\link{normalizeCounts}}
#' @return a data.table
#' @export

LACHESIS <- function(input.files = NULL, ids = NULL, cnv.files = NULL, snv.files = NULL, vcf.source = NULL,
                     purity = NULL, ploidy = NULL,
                     cnv.chr.col = NULL, cnv.start.col = NULL, cnv.end.col = NULL, cnv.A.col = NULL,
                     cnv.B.col = NULL, cnv.tcn.col = NULL, age = NULL,
                     OS.time = NULL, OS = NULL, EFS.time = NULL, EFS = NULL, output.dir = NULL,  ...){

  if(is.null(input.files) & is.null(cnv.files)){
    stop("Missing input file!")
  }else if(is.null(input.files)){
    if(any(is.null(cnv.files), is.null(snv.files))){
      stop("Missing snv and cnv inputs!")
    }else if(length(cnv.files) != length(snv.files)){
      stop("Please provide snv and cnv input for every sample!")
    }
  }

  # collect ECA and MRCA densities for each tumor
  cohort.densities <- data.table::data.table(Sample_ID = character(),
                                             MRCA_time_mean = numeric(),
                                             MRCA_time_lower = numeric(),
                                             MRCA_time_upper = numeric(),
                                             ECA_time_mean = numeric(),
                                             ECA_time_lower = numeric(),
                                             ECA_time_upper = numeric(),
                                             Ploidy = numeric(),
                                             Purity = numeric(),
                                             Age = numeric(),
                                             OS.time = numeric(),
                                             OS = numeric(),
                                             EFS.time = numeric(),
                                             EFS = numeric())

  if(!is.null(input.files)){

    sample.specs <- data.table::fread(input.files, sep = "\t", header = T, stringsAsFactors = F)

    if(any(is.na(sample.specs[,ID]))){
      warning("No sample name provided for samples ", sample.specs[,which(is.na(ID))], "; sample name was set to 1 - ", sample.specs[,sum(is.na(ID))])
      sample.specs[, ID := as.character(ID)][is.na(ID), ID := which(is.na(ID))]
    }

    if(any(is.na(sample.specs[,cnv.file]))){
      warning("No CNV file provided for sample(s) ", paste(sample.specs[,ID[which(is.na(cnv.file))]], collapse=", "), "; sample(s) will be excluded")
      sample.specs[!is.na(cnv.file), ]
      if(nrow(sample.specs)==0){
        stop("No files retained! Stopping analysis.")
      }
    }

    if(any(is.na(sample.specs[,snv.file]))){
      warning("No SNV file provided for sample(s) ", paste(sample.specs[,ID[which(is.na(snv.file))]], collapse=", "), "; sample(s) will be excluded")
      sample.specs[!is.na(snv.file), ]
      if(nrow(sample.specs)==0){
        stop("No files retained! Stopping analysis.")
      }
    }

    sample.specs.spl <- split(sample.specs, sample.specs$ID)

    for(i in 1:length(sample.specs.spl)){

      x <- sample.specs.spl[[i]]
      x[,which(sapply(x, is.na)):=NULL] # remove NA entries

      if(is.null(x$ID)){
        stop()
      }
      message("Computing SNV density for sample ", x$ID)

      if(!is.null(output.dir)){
        dir.create(paste(output.dir, x$ID, sep="/")) # create per-sample output directory
      }else{
        warning("No output directory specified. Per-sample output will be discarded.")
      }

      cnv <- readCNV(cn.info = x$cnv.file, chr.col = x$cnv.chr.col, start.col = x$cnv.start.col,
                     end.col = x$cnv.end.col, A.col = x$cnv.A.col, B.col = x$cnv.B.col,
                     tcn.col = x$cnv.tcn.col, tumor.id = x$ID, ...)

      snv <- readVCF(vcf = x$snv.file, vcf.source = x$vcf.source, t.sample = x$ID, ...)

      nb <- nbImport(cnv = cnv, snv = snv, purity = x$purity, ploidy = x$ploidy)
      nb.p1 <- plotNB(nb = nb, samp.name = x$ID, ...)

      if(!is.null(output.dir)){
        pdf(paste(output.dir, x$ID, "VAF_histogram.pdf", sep="/"), width = 8, height = 6)
        print(nb.p1)
        dev.off()
      }

      raw.counts <- clonalMutationCounter(nbObj = nb, ...)
      norm.counts <- normalizeCounts(countObj = raw.counts)
      mrca <- MRCA(normObj = norm.counts, ...)

      this.tumor.density <- data.table::data.table(Sample_ID = x$ID,
                                                   MRCA_time_mean = attributes(mrca)$MRCA_time_mean,
                                                   MRCA_time_lower = attributes(mrca)$MRCA_time_lower,
                                                   MRCA_time_upper = attributes(mrca)$MRCA_time_upper,
                                                   ECA_time_mean = attributes(mrca)$ECA_time_mean,
                                                   ECA_time_lower = attributes(mrca)$ECA_time_lower,
                                                   ECA_time_upper = attributes(mrca)$ECA_time_upper,
                                                   Ploidy = x$ploidy,
                                                   Purity = x$purity,
                                                   Age = x$Age,
                                                   OS.time = x$OS.time,
                                                   OS = x$OS,
                                                   EFS.time = x$EFS.time,
                                                   EFS = x$EFS)

      cohort.densities <- merge(cohort.densities, this.tumor.density, all=T)

      # output the result for this sample
      if(!is.null(output.dir)){
        write(attributes(mrca)[c("purity", "ploidy", "MRCA_time_mean", "MRCA_time_lower", "MRCA_time_upper", "ECA_time_mean", "ECA_time_lower", "ECA_time_upper")],
              file = paste(output.dir, x$ID, paste("MRCA_densities_", x$ID, ".txt", sep=""), sep="/"))
        write.table(mrca, file = paste(output.dir, x$ID, paste("SNV_timing_per_segment_", x$ID, ".txt", sep=""), sep="/"))
      }

      mrca.p2 <- plotMutationDensities(mrcaObj = mrca, samp.name = x$ID, ...)

      if(!is.null(output.dir)){
        pdf(paste(output.dir, x$ID, "SNV_densities.pdf", sep="/"), width = 8, height = 6)
        print(mrca.p2)
        dev.off()
      }
    }
    rm(sample.specs.spl)
  }else{

    if(any(is.na(ids))){
      warning("No sample name provided for samples ", which(is.na(ids)), "; sample name was set to 1 - ", sum(is.na(ids)))
      ids[is.na(ids)] <- which(is.na(ids))
    }

    sample.specs.spl <- split(sample.specs, sample.specs$ID)

    for(i in 1:length(cnv.files)){

      message("Computing SNV density for sample ", x$ID)

      if(is.na(cnv.files)[i]){
        warning("No CNV file provided for sample ", ids[i], "; sample will be excluded")
        next
      }
      if(is.na(cnv.files)[i]){
        warning("No SNV file provided for sample ", ids[i], "; sample will be excluded")
        next
      }

      x <- sample.specs.spl[[i]]
      x[,which(sapply(x, is.na)):=NULL] # remove NA entries

      if(is.null(x$ID)){
        stop()
      }
      message("Computing SNV density for sample ", x$ID)

      cnv <- readCNV(cn.info = cnv.files[i], chr.col = ifelse(is.na(cnv.chr.col[i]), NULL, cnv.chr.col[i]),
                     start.col = ifelse(is.na(cnv.start.col[i]), NULL, cnv.start.col[i]),
                     end.col = ifelse(is.na(cnv.end.col[i]), NULL, cnv.end.col[i]),
                     A.col = ifelse(is.na(cnv.A.col[i]), NULL, cnv.A.col[i]),
                     B.col = ifelse(is.na(cnv.B.col[i]), NULL, cnv.B.col[i]),
                     tcn.col = ifelse(is.na(cnv.tcn.col[i]), NULL, cnv.tcn.col[i]), tumor.id = ids[i], ...)

      snv <- readVCF(vcf = snv.files[i], vcf.source = vcf.source[i], t.sample = ids[i], ...)

      nb <- nbImport(cnv = cnv, snv = snv, purity = purity[i], ploidy = ploidy[i])
      nb.p1 <- plotNB(nb = nb, samp.name = ids[i], ...)

      if(!is.null(output.dir)){
        pdf(paste(output.dir, x$ID, "VAF_histogram.pdf", sep="/"), width = 8, height = 6)
        print(nb.p1)
        dev.off()
      }

      raw.counts <- clonalMutationCounter(nbObj = nb, ...)
      norm.counts <- normalizeCounts(countObj = raw.counts)
      mrca <- MRCA(normObj = norm.counts, ...)

      # output the result for this sample
      if(!is.null(output.dir)){
        write(attributes(mrca)[c("purity", "ploidy", "MRCA_time_mean", "MRCA_time_lower", "MRCA_time_upper", "ECA_time_mean", "ECA_time_lower", "ECA_time_upper")],
              file = paste(output.dir, x$ID, paste("MRCA_densities_", x$ID, ".txt", sep=""), sep="/"))
        write.table(mrca, file = paste(output.dir, x$ID, paste("SNV_timing_per_segment_", x$ID, ".txt", sep=""), sep="/"))
      }

      mrca.p2 <- plotMutationDensities(mrcaObj = mrca, samp.name = x$ID, ...)

      if(!is.null(output.dir)){
        pdf(paste(output.dir, x$ID, "SNV_densities.pdf", sep="/"), width = 8, height = 6)
        print(mrca.p2)
        dev.off()
      }

      this.tumor.density <- data.table::data.table(Sample_ID = ids[i],
                                                   MRCA_time_mean = attributes(mrca)$MRCA_time_mean,
                                                   MRCA_time_lower = attributes(mrca)$MRCA_time_lower,
                                                   MRCA_time_upper = attributes(mrca)$MRCA_time_upper,
                                                   ECA_time_mean = attributes(mrca)$ECA_time_mean,
                                                   ECA_time_lower = attributes(mrca)$ECA_time_lower,
                                                   ECA_time_upper = attributes(mrca)$ECA_time_upper,
                                                   Ploidy = ploidy[i],
                                                   Purity = purity[i],
                                                   Age = age[i],
                                                   OS.time = OS.time[i],
                                                   OS = OS[i],
                                                   EFS.time = EFS.time[i],
                                                   EFS = EFS[i])

      cohort.densities <- merge(cohort.densities, this.tumor.density, all=T)
    }
  }

  # plot the distribution of Mutation densities at ECA and MRCA

  return(cohort.densities)
}


#' Plot SNV densities at ECA and MRCA
#' @description
#' Visualizes results from \code{\link{LACHESIS}}. Top plot, histograms of mean mutation densities; optional, correlation with age at diagnosis; bottom plots, cumulative distribution of mean mutation densities with 95% confidence intervals.
#' @param lachesis output generated from \code{\link{LACHESIS}}
#' @param suppress.outliers shall outliers (defined as the 2.5% tumors with lowest and highest densities) be plotted? Default `TRUE`.
#' @param log.densities plot logarithmic densities. Default `FALSE`
#' @param fill.multi  optional, the color code for ECA and MRCA.
#' @param l.col, optional, the line color
#' @param binwidth optional; the bin-width in the histogram.
#' @examples
#' input.files = system.file("extdata", "Sample_template.txt")
#' lachesis <- LACHESIS(input.files = input.files)
#' plotLachesis(lachesis)
#' @export
#' @importFrom graphics abline Axis box grid hist mtext par rect text title arrows legend points polygon

plotLachesis <- function(lachesis = NULL, suppress.outliers = FALSE, log.densities = FALSE, fill.multi = NULL, l.col = NULL, binwidth = NULL){

  if(is.null(lachesis)){
    stop("Missing input. Please provide the output generated by LACHESIS()")
  }
  # graphical settings:
  if(is.null(l.col)){
    l.col <- NA
  }
  if(is.null(fill.zero)){
    fill.zero <- "#4FB12B"
  }
  if(is.null(fill.multi)){
    fill.multi <- "#176A02"
  }

  # I. plot histograms
  to.plot <- lachesis

  if(suppress.outliers){
    to.plot <- to.plot[MRCA_time_mean < quantile(MRCA_time_mean, 0.975) &
              MRCA_time_mean > quantile(MRCA_time_mean, 0.025),]
  }

  if(is.null(binwidth)){
    binwidth = (max(to.plot$MRCA_time_mean) - min(to.plot$MRCA_time_mean))/20
  }

  lo_mat <- matrix(data = c(1, 2, 3, 4, 5, 5, 6, 6), nrow = 2, ncol=4, byrow = TRUE)
  graphics::layout(mat = lo_mat, widths = c(1, 1, 1, 1, 2, 2), heights = c(1, 1, 1, 1, 1, 1))

  par(mar = c(3, 4, 3, 1))

  hist(to.plot[,MRCA_time_mean], xlim = c(0, 1.05 * max(to.plot[,MRCA_time_mean])),
       breaks = seq(0, max(to.plot[,MRCA_time_mean])*1.05, binwidth), col = fill.zero, border = l.col, main = NA,
       xlab = NA, ylab = NA)

  title(main = paste("SNV densities at MRCA"), cex.main = 1.2)
  mtext(text = "SNVs per Mb", side = 1, line = 2.5, cex = 0.8)
  mtext(text = "No. of tumors", side = 2, line = 1.8, cex = 0.8)

  par(mar = c(3, 4, 3, 1))

  binwidth = (max(to.plot$ECA_time_mean) - min(to.plot$ECA_time_mean))/20

  hist(to.plot[,ECA_time_mean], xlim = c(0, 1.05 * max(to.plot[,ECA_time_mean])),
       breaks = seq(0, max(to.plot[,ECA_time_mean])*1.05, binwidth), col = fill.multi, border = l.col, main = NA,
       xlab = NA, ylab = NA)

  title(main = paste("SNV densities at ECA"), cex.main = 1.2)
  mtext(text = "SNVs per Mb", side = 1, line = 2.5, cex = 0.8)
  mtext(text = "No. of tumors", side = 2, line = 1.8, cex = 0.8)

  # Correlation between SNV density and Age at diagnosis
  to.plot <- to.plot[!is.na(to.plot$Age),]
  if(nrow(to.plot)>0){
    par(mar = c(3, 1, 3, 1), xpd = FALSE)

    to.plot[,plot(MRCA_time_mean, Age, xlab = NA, ylab = NA, xlim = c(0, 1.05*max(MRCA_time_mean)),
                  ylim = c(0, 1.05*max(Age)))]
    title(main = "SNV densities at MRCA vs age", cex.main = 1.2)
    mtext(text = "SNVs per Mb", side = 1, line = 2.5, cex = 0.8)
    mtext(text = "Age at diagnosis", side = 2, line = 1.8, cex = 0.8)

    par(mar = c(3, 1, 3, 1), xpd = FALSE)

    to.plot[,plot(ECA_time_mean, Age, xlab = NA, ylab = NA, xlim = c(0, 1.05*max(ECA_time_mean)),
                  ylim = c(0, 1.05*max(Age)))]
    title(main = "SNV densities at ECA vs age", cex.main = 1.2)
    mtext(text = "SNVs per Mb", side = 1, line = 2.5, cex = 0.8)
    mtext(text = "Age at diagnosis", side = 2, line = 1.8, cex = 0.8)
  }else{
    # empty plots
    plot.new()
    plot.new()
  }

  # Cumulative densities at ECA/MRCA
  par(mar = c(3, 1, 3, 1), xpd = FALSE)

  to.plot <- .ecdf.stats(lachesis, mean = "MRCA_time_mean", lower = "MRCA_time_lower", upper = "MRCA_time_upper")

  x.min = 0
  x.max = max(c(lachesis$MRCA_time_upper))*1.3
  y.min = 0
  y.max = 1
  plot(NA, NA, xlim=c(x.min, x.max), ylim=c(y.min, y.max), xlab = NA, ylab = NA, main = NA, axes = FALSE, frame.plot = FALSE)
  Axis(side=1, cex = 0.7)
  Axis(side=2, cex = 0.7)
  mtext(text = "SNVs per Mb", side = 1, line = 2, cex = 0.7)
  mtext(text = "Number of tumors", side = 2, line = 2, cex = 0.7)

  polygon(c(sort(lachesis$MRCA_time_mean), sort(lachesis$MRCA_time_mean, decreasing = T)),
          c(sort(to.plot$ecdf.upper), sort(to.plot$ecdf.lower, decreasing = T)),
          col = fill.zero, border = NA)
  plot.ecdf(lachesis$MRCA_time_mean, col = "black", add=T)

  # ECA:
  par(mar = c(3, 1, 3, 1), xpd = FALSE)

  to.plot <- .ecdf.stats(lachesis, mean = "ECA_time_mean", lower = "ECA_time_lower", upper = "ECA_time_upper")

  x.min = 0
  x.max = max(c(lachesis$ECA_time_upper))*1.3
  y.min = 0
  y.max = 1
  plot(NA, NA, xlim=c(x.min, x.max), ylim=c(y.min, y.max), xlab = NA, ylab = NA, main = NA, axes = FALSE, frame.plot = FALSE)
  Axis(side=1, cex = 0.7)
  Axis(side=2, cex = 0.7)
  mtext(text = "SNVs per Mb", side = 1, line = 2, cex = 0.7)
  mtext(text = "Number of tumors", side = 2, line = 2, cex = 0.7)

  polygon(c(sort(lachesis$ECA_time_mean), sort(lachesis$ECA_time_mean, decreasing = T)),
          c(sort(to.plot$ecdf.upper), sort(to.plot$ecdf.lower, decreasing=T)), col = fill.zero, border = NA)
  plot.ecdf(lachesis$ECA_time_mean, col = "black", add = T)

  title(main = paste("SNV densities at ECA/MRCA"), cex.main = 1.2)

}


#' Compute ECDF statistics with lower and upper bounds
#' @param data data table containing mean, lower and upper bound of value of interest.
#' @param mean column name containing mean values.
#' @param lower column name containing lower bounds.
#' @param upper  column name containing upper bounds.
#' @importFrom data.table

.ecdf.stats <- function(data, mean, lower, upper){
  if(nrow(data)==0){
    return(NULL)
  }
  l <- ecdf(unlist(data[,..lower]))
  u <- ecdf(unlist(data[,..upper]))
  m <- ecdf(unlist(data[,..mean]))
  data.table::data.table(ecdf.lower = l(unlist(data[,..mean])),
              ecdf.mean = m(unlist(data[,..mean])),
              ecdf.upper = u(unlist(data[,..mean])))
}
