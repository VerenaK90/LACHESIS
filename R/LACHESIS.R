#' Run MRCA density estimation for a set of tumors
#' @description
#' Takes a set of SNV and CNV files as input and outputs per-tumor SNV densities. Input can either be a tab-delimited file containing
#' the sample specifications or vectors giving direct paths to the sample files. CNV file requires columns for the chromosome number, start and end of the segment, and either the total copy number or the number of A- and B-alleles
#' @param input.files a tab-delimited sample-specification file, it must contain the sample name, the path to the SNV file, path to CNV file, and optionally purity, ploidy, cnv.chr.col, cnv.start.col, cnv.end.col, cnv.A.col, cnv.B.col, cnv.tcn.col. A template for this spreadsheet can be downloaded from ...
#' @param ids vector of sample names, will be ignored if `input.files` is specified.
#' @param vcf.tumor.ids vector of sample names as given in the vcf file; will be ignored if `input.files` is specified.
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
#' @param min.cn minimum copy number to be included in the analysis. Default 2.
#' @param max.cn maximum copy number to be included in the analysis. Default 4.
#' @param merge.tolerance the maximum distance below which adjacent segments with equal copy number are merged. Defaults to 10^5 bp.
#' @param ignore.XY Ignore allosomes. Default TRUE
#' @param min.vaf Remove variants with vcf below threshold. Default 0.01
#' @param min.depth Minimum required depth for a variant to be considered. Default 30.
#' @param vcf.info.af The string encoding the allele frequency field in the FORMAT column of the .vcf file. Defaults to `AF`and will be ignored if `vcf.source` != `sentieon`.
#' @param vcf.info.dp The string encoding the read depth field in the FORMAT column of the .vcf file. Defaults to `DP`and will be ignored if `vcf.source` != `sentieon`.
#' @param min.seg.size the minimal segment length to be included in the quantification
#' @param fp.mean optional, the average false positive rate of clonal mutations (e.g., due to incomplete tissue sampling). Defaults to 0.
#' @param fp.sd optional, the standard deviation of the false positive rate of clonal mutations (e.g., due to incomplete tissue sampling). Defaults to 0.
#' @param excl.chr a vector of chromosomes that should be excluded from the quantification. e.g., due to reporter constructs in animal models.
#' @param ref_build Reference genome. Default `hg19`. Can be `hg18`, `hg19` or `hg38`
#' @param ... further arguments and parameters passed to `plotMutationDensities`.
#' @examples
#' #an example file with sample annotations and meta data
#' input.files = system.file("extdata", "Sample_template.txt", package = "LACHESIS")
#' input.files = data.table::fread(input.files)
#'
#' #cnv and snv files for example tumors
#' nbe11 = list.files(system.file("extdata/NBE11/", package = "LACHESIS"), full.names = TRUE)
#' nbe15 = list.files(system.file("extdata/NBE15/", package = "LACHESIS"), full.names = TRUE)
#' nbe63 = list.files(system.file("extdata/NBE63/", package = "LACHESIS"), full.names = TRUE)
#'
#' cnv.file = c(nbe11[1], nbe15[1], nbe63[1])
#' snv.file = c(nbe11[2], nbe15[2], nbe63[2])
#'
#' input.files$cnv.file = cnv.file
#' input.files$snv.file = snv.file
#'
#' # Make an axample input file with paths to cnv and snv file along with other meta data
#' lachesis_input = tempfile(pattern = "lachesis", tmpdir = tempdir(), fileext = ".tsv")
#' data.table::fwrite(x = input.files, file = lachesis_input, sep = "\t")
#'
#' #Exampele with template file with paths to multiple cnv/snv files as an input
#' lachesis <- LACHESIS(input.files = lachesis_input)
#'
#' #Example with a single sample input
#' strelka_vcf = system.file("extdata","strelka2.somatic.snvs.vcf.gz", package = "LACHESIS")
#' aceseq_cn = system.file("extdata", "ACESeq/NBE11_comb_pro_extra2.59_0.83.txt", package = "LACHESIS")
#' lachesis <- LACHESIS(ids = "NBE11", cnv.files = aceseq_cn, snv.files = strelka_vcf, vcf.source = "strelka", purity = 0.83, ploidy = 2.59)
#' @seealso \code{\link{MRCA}} \code{\link{clonalMutationCounter}} \code{\link{normalizeCounts}}
#' @import tidyr
#' @return a data.table
#' @export

LACHESIS <- function(input.files = NULL, ids = NULL, cnv.files = NULL, snv.files = NULL, vcf.source = NULL,
                     purity = NULL, ploidy = NULL,
                     cnv.chr.col = NULL, cnv.start.col = NULL, cnv.end.col = NULL, cnv.A.col = NULL,
                     cnv.B.col = NULL, cnv.tcn.col = NULL, age = NULL,
                     OS.time = NULL, OS = NULL, EFS.time = NULL, EFS = NULL, output.dir = NULL,
                     ignore.XY = TRUE, min.cn = 1, max.cn = 4, merge.tolerance = 10^5, min.vaf = 0.01, min.depth = 30,
                     vcf.info.af = "AF", vcf.info.dp = "DP", min.seg.size = 10^7, fp.mean = 0, fp.sd = 0, excl.chr = NULL,
                     ref_build = "hg19", ...){


  ID <- cnv.file <- snv.file <- write.table <- NULL

  if(is.null(input.files) & is.null(cnv.files)){
    stop("Missing input file!")
  }else if(is.null(input.files)){
    if(any(is.null(cnv.files), is.null(snv.files))){
      stop("Missing snv and cnv inputs!")
    }else if(length(cnv.files) != length(snv.files)){
      stop("Please provide snv and cnv input for every sample!")
    }
  }

  incl.chr <- setdiff(c(1:22), excl.chr)
  if(!ignore.XY){
    incl.chr <- c(incl.chr, "X", "Y")
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

    sample.specs <- data.table::fread(input.files, sep = "\t", stringsAsFactors = FALSE)

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

    if(any(lapply(sample.specs.spl, nrow) > 1)){
      stop("Duplicated IDs found!")
    }

    for(i in 1:length(sample.specs.spl)){

      x <- sample.specs.spl[[i]]
      x[,which(sapply(x, is.na)):=NULL] # remove NA entries

      if(is.null(x$ID)){
        stop("Please provide sample identifiers.")
      }
      if(is.null(x$vcf.source)){
        x$vcf.source <- x$ID
      }else if(any(is.na(x$vcf.source))){
        warning("No column identifier provided for sample ", which(is.na(x$vcf.source)), "; will be inferred.")
        x$vcf.source[is.na(x$vcf.source)] <- x$id[is.na(x$vcf.source)]
      }
      message("Computing SNV density for sample ", x$ID)

      if(!is.null(output.dir)){
        dir.create(paste(output.dir, x$ID, sep="/"), recursive = TRUE, showWarnings = FALSE) # create per-sample output directory
      }else{
        warning("No output directory specified. Per-sample output will be discarded.")
      }

      cnv <- readCNV(cn.info = x$cnv.file, chr.col = x$cnv.chr.col, start.col = x$cnv.start.col,
                     end.col = x$cnv.end.col, A.col = x$cnv.A.col, B.col = x$cnv.B.col,
                     tcn.col = x$cnv.tcn.col, tumor.id = x$ID, merge.tolerance = merge.tolerance,
                     max.cn = max.cn, ignore.XY = ignore.XY)

      snv <- readVCF(vcf = x$snv.file, vcf.source = x$vcf.source, t.sample = x$vcf.tumor.ids, min.depth = min.depth,
                     min.vaf = min.vaf, info.af = vcf.info.af, info.dp = vcf.info.dp)
      vaf.p1 <- plotVAFdistr(snv)

      nb <- nbImport(cnv = cnv, snv = snv, purity = x$purity, ploidy = x$ploidy)

      if(!is.null(output.dir)){
        plotVAFdistr(snv, output.file = paste(output.dir, x$ID, "VAF_histogram.pdf", sep="/"))
        plotNB(nb = nb, samp.name = x$ID, output.file = paste(output.dir, x$ID, "VAF_histogram_strat.pdf", sep="/"), ref_build = ref_build, min.cn = min.cn, max.cn = max.cn)
      }

      raw.counts <- clonalMutationCounter(nbObj = nb, min.cn = min.cn, max.cn = max.cn, chromosomes = incl.chr)
      norm.counts <- normalizeCounts(countObj = raw.counts)
      mrca <- MRCA(normObj = norm.counts, min.seg.size = min.seg.size, fp.mean = fp.mean, excl.chr = excl.chr)

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
        write.table(unlist(attributes(mrca)[c("purity", "ploidy", "MRCA_time_mean", "MRCA_time_lower", "MRCA_time_upper", "ECA_time_mean", "ECA_time_lower", "ECA_time_upper")]),
              file = paste(output.dir, x$ID, paste("MRCA_densities_", x$ID, ".txt", sep=""), sep="/"), quote = F, col.names = F, sep="\t")
        write.table(mrca, file = paste(output.dir, x$ID, paste("SNV_timing_per_segment_", x$ID, ".txt", sep=""), sep="/"), row.names = F, quote=F, sep="\t")
      }

      if(!is.null(output.dir)){
        plotMutationDensities(mrcaObj = mrca, samp.name = x$ID, output.file = paste(output.dir, x$ID, "SNV_densities.pdf", sep="/"), ...)
      }
    }
    rm(sample.specs.spl)
  }else{

    if(any(is.na(ids))){
      warning("No sample name provided for samples ", which(is.na(ids)), "; sample name was set to 1 - ", sum(is.na(ids)))
      ids[is.na(ids)] <- which(is.na(ids))
    }
    if(is.null(vcf.tumor.ids)){
      warning("No column identifiers provided.")
      vcf.tumor.ids <- ids
    }else if(any(is.na(vcf.tumor.ids))){
      warning("No column identifier  provided for samples ", which(is.na(vcf.tumor.ids)), "; column name will be inferred")
      vcf.tumor.ids[is.na(vcf.tumor.ids)] <- ids[is.na(vcf.tumor.ids)]
    }

    for(i in 1:length(cnv.files)){

      message("Computing SNV density for sample ", ids[i])

      if(!is.null(output.dir)){
        dir.create(paste(output.dir, ids[i], sep="/"), recursive = TRUE, showWarnings = FALSE) # create per-sample output directory
      }else{
        warning("No output directory specified. Per-sample output will be discarded.")
      }

      if(is.na(cnv.files)[i]){
        warning("No CNV file provided for sample ", ids[i], "; sample will be excluded")
        next
      }
      if(is.na(snv.files)[i]){
        warning("No SNV file provided for sample ", ids[i], "; sample will be excluded")
        next
      }


      cnv <- readCNV(cn.info = cnv.files[i], chr.col = cnv.chr.col[i], start.col = cnv.start.col[i],
                     end.col = cnv.end.col[i], A.col = cnv.A.col[i], B.col = cnv.B.col[i],
                     tcn.col = cnv.tcn.col[i], tumor.id = ids[i], merge.tolerance = merge.tolerance,
                     max.cn = max.cn, ignore.XY = ignore.XY)

      snv <- readVCF(vcf = snv.files[i], vcf.source = vcf.source[i], t.sample = vcf.tumor.ids[i], min.depth = min.depth,
                     min.vaf = min.vaf, info.af = vcf.info.af, info.dp = vcf.info.dp)
      vaf.p <- plotVAFdistr(snv)

      nb <- nbImport(cnv = cnv, snv = snv, purity = purity[i], ploidy = ploidy[i])

      if(!is.null(output.dir)){
        plotVAFdistr(snv, output.file = paste(output.dir, ids[i], "VAF_histogram.pdf", sep="/"))
        plotNB(nb = nb, samp.name = ids[i], output.file = paste(output.dir, ids[i], "VAF_histogram_strat.pdf", sep="/"), ref_build = ref_build, min.cn = min.cn, max.cn = max.cn)
      }

      raw.counts <- clonalMutationCounter(nbObj = nb, min.cn = min.cn, max.cn = max.cn, chromosomes = incl.chr)
      norm.counts <- normalizeCounts(countObj = raw.counts)
      if(nrow(norm.counts)==1){
        warning("Too few segments to estimate MRCA density for sample ", ids[i], ".")
        mrca <- ""
        attr(mrca, "MRCA_time_mean") <- NA
        attr(mrca, "MRCA_time_upper") <- NA
        attr(mrca, "MRCA_time_lower") <- NA
        attr(mrca, "ECA_time_mean") <- NA
        attr(mrca, "ECA_time_lower") <- NA
        attr(mrca, "ECA_time_upper") <- NA
      }else{
        mrca <- MRCA(normObj = norm.counts, min.seg.size = min.seg.size, fp.mean = fp.mean, excl.chr = excl.chr)

        # output the result for this sample
        if(!is.null(output.dir)){
          write(unlist(attributes(mrca)[c("purity", "ploidy", "MRCA_time_mean", "MRCA_time_lower", "MRCA_time_upper", "ECA_time_mean", "ECA_time_lower", "ECA_time_upper")]),
                file = paste(output.dir, ids[i], paste("MRCA_densities_", ids[i], ".txt", sep=""), sep="/"))
          write.table(mrca, file = paste(output.dir, ids[i], paste("SNV_timing_per_segment_", ids[i], ".txt", sep=""), sep="/"))
        }

        if(!is.null(output.dir)){
          plotMutationDensities(mrcaObj = mrca, samp.name = ids[i], output.file = paste(output.dir, ids[i], "SNV_densities.pdf", sep="/"), ...)
        }
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

      cohort.densities <- merge(cohort.densities, this.tumor.density, all=TRUE)
    }
  }

  # plot the distribution of Mutation densities at ECA and MRCA

  if(!is.null(output.dir)){
    plotLachesis(cohort.densities, output.file = paste(output.dir, "SNV_densities_cohort.pdf", sep="/"))
  }

  return(cohort.densities)
}


#' Plot SNV densities at ECA and MRCA
#' @description
#' Visualizes results from \code{\link{LACHESIS}}. Top plot, histograms of mean mutation densities; bottom plots, cumulative distribution of mean mutation densities with 95% confidence intervals.
#' @param lachesis output generated from \code{\link{LACHESIS}}
#' @param suppress.outliers shall outliers (defined as the 2.5% tumors with lowest and highest densities) be plotted? Default `TRUE`.
#' @param log.densities plot logarithmic densities. Default `FALSE`
#' @param fill.zero optional, the color code for MRCA.
#' @param fill.multi optional, the color code for ECA.
#' @param l.col, optional, the line color
#' @param binwidth optional; the bin-width in the histogram.
#' @param output.file optional; the file to which the plot will be stored.
#' @examples
#' #an example file with sample annotations and meta data
#' input.files = system.file("extdata", "Sample_template.txt", package = "LACHESIS")
#' input.files = data.table::fread(input.files)
#'
#' #cnv and snv files for example tumors
#' nbe11 = list.files(system.file("extdata/NBE11/", package = "LACHESIS"), full.names = TRUE)
#' nbe15 = list.files(system.file("extdata/NBE15/", package = "LACHESIS"), full.names = TRUE)
#' nbe63 = list.files(system.file("extdata/NBE63/", package = "LACHESIS"), full.names = TRUE)
#'
#' cnv.file = c(nbe11[1], nbe15[1], nbe63[1])
#' snv.file = c(nbe11[2], nbe15[2], nbe63[2])
#'
#' input.files$cnv.file = cnv.file
#' input.files$snv.file = snv.file
#'
#' # Make an example input file with paths to cnv and snv file along with other meta data
#' lachesis_input = tempfile(pattern = "lachesis", tmpdir = tempdir(), fileext = ".tsv")
#' data.table::fwrite(x = input.files, file = lachesis_input, sep = "\t")
#'
#' #Example with template file with paths to multiple cnv/snv files as an input
#' lachesis <- LACHESIS(input.files = lachesis_input)
#' plotLachesis(lachesis)
#' @export
#' @importFrom graphics abline Axis box grid hist mtext par rect text title arrows legend points polygon

plotLachesis <- function(lachesis = NULL, suppress.outliers = FALSE, log.densities = FALSE, fill.multi = NULL, l.col = NULL, binwidth = NULL, fill.zero = NULL, output.file = NULL){

  MRCA_time_mean <- ECA_time_mean <- NULL

  if(is.null(lachesis)){
    stop("Missing input. Please provide the output generated by LACHESIS()")
  }
  if(nrow(lachesis)==1){
    warning("Cannot produce summary statistics for a single case. Returning null.")
    return()
  }
  if(any(is.na(lachesis$MRCA_time_mean))){
    warning("Removing ", sum(is.na(lachesis$MRCA_time_mean)), " samples with missing MRCA density estimate.")
    lachesis <- lachesis[!is.na(MRCA_time_mean),]
  }
  if(nrow(lachesis)==0){
    warning("No sample with MRCA density estimate provided. Returning zero.")
    return(NULL)
  }
  if(!is.null(output.file)){
    pdf(output.file, width = 8, height = 6)
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
    if(log.densities){
      binwidth = (max(log10(to.plot$MRCA_time_mean)) - min(log10(to.plot$MRCA_time_mean)))/20
    }else{
      binwidth = (max(to.plot$MRCA_time_mean) - min(to.plot$MRCA_time_mean))/20
    }

  }

  lo_mat <- matrix(data = c(1, 2, 3, 4), nrow = 2, ncol=2, byrow = TRUE)
  graphics::layout(mat = lo_mat, widths = c(1, 2, 1, 2), heights = c(1, 1, 1, 1))

  par(mar = c(3, 4, 3, 1))

  if(log.densities){

    min.x <- floor(min(to.plot[,log10(MRCA_time_mean)]))
    max.x <- ceiling(max(to.plot[,log10(MRCA_time_mean)]))

    hist(to.plot[,log10(MRCA_time_mean)], xlim = c(min.x, max.x),
         breaks = 20,
         col = fill.zero, border = l.col, main = NA,
         xlab = NA, ylab = NA, axes = FALSE)

    Axis(side = 1, at = seq(min.x, max.x, length.out = 10),
         labels = round(10^seq(min.x, max.x, length.out = 10), digits = 2))
    Axis(side = 2)

  }else{
    hist(to.plot[,MRCA_time_mean], xlim = c(0, 1.05 * max(to.plot[,MRCA_time_mean])),
         breaks = seq(0, max(to.plot[,MRCA_time_mean])*1.05, binwidth), col = fill.zero, border = l.col, main = NA,
         xlab = NA, ylab = NA)
  }


  title(main = paste("SNV densities at MRCA"), cex.main = 1)
  mtext(text = "SNVs per Mb", side = 1, line = 2, cex = 0.7)
  mtext(text = "No. of tumors", side = 2, line = 1.8, cex = 0.7)

  # Cumulative densities at MRCA
  par(mar = c(3, 4, 3, 1), xpd = FALSE)

  x.min = 0
  x.max = max(c(lachesis$MRCA_time_upper))*1.3
  y.min = 0
  y.max = 1
  plot(NA, NA, xlim=c(x.min, x.max), ylim=c(y.min, y.max), xlab = NA, ylab = NA, main = NA, axes = FALSE, frame.plot = FALSE)
  Axis(side=1, cex = 0.7)
  Axis(side=2, cex = 0.7)
  mtext(text = "SNVs per Mb", side = 1, line = 2, cex = 0.7)
  mtext(text = "No. of tumors", side = 2, line = 2, cex = 0.7)

  to.plot <- data.frame(x.lower = rep(sort(c( lachesis$MRCA_time_mean)), each = 2)[-1],
                        x.upper = rep(sort(c( lachesis$MRCA_time_mean)), each = 2)[-2*(nrow(lachesis) )])
  to.plot$y.lower <- sapply(rep(sort(c( lachesis$MRCA_time_mean)), each = 2), function(x){sum(lachesis$MRCA_time_upper <= x)})[-2*(nrow(lachesis) )]
  to.plot$y.upper <- sapply(rep(sort(c( lachesis$MRCA_time_mean)), each = 2), function(x){sum(lachesis$MRCA_time_lower <= x)})[-1]

  polygon(c(to.plot$x.lower, rev(to.plot$x.upper)),
          c(to.plot$y.lower, rev(to.plot$y.upper))/nrow(lachesis),
          col = fill.zero, border = NA)

  plot.ecdf(lachesis$MRCA_time_mean, col = "black", add=TRUE, verticals = TRUE)

  title(main = paste("SNV densities at MRCA"), cex.main = 1)

  # Histogram of SNV density at ECA

  if(all(is.na(lachesis$ECA_time_mean)) & !is.null(output.file)){
    dev.off()
    return()
  }else if(all(is.na(lachesis$ECA_time_mean))){
    return()
  }

  to.plot <- lachesis

  if(suppress.outliers){
    to.plot <- to.plot[MRCA_time_mean < quantile(MRCA_time_mean, 0.975) &
                         MRCA_time_mean > quantile(MRCA_time_mean, 0.025),]
  }

  par(mar = c(3, 4, 3, 1))

  if(log.densities){
    min.x <- floor(min(to.plot[,log10(ECA_time_mean)], na.rm = TRUE))
    max.x <- ceiling(max(to.plot[,log10(ECA_time_mean)], na.rm = TRUE))

    hist(to.plot[,log10(ECA_time_mean)], xlim = c(min.x, max.x),
         breaks = 20,
         col = fill.zero, border = l.col, main = NA,
         xlab = NA, ylab = NA, axes = FALSE)

    Axis(side = 1, at = seq(min.x, max.x, length.out = 10),
         labels = round(10^seq(min.x, max.x, length.out = 10), digits = 2))
    Axis(side = 2)
  }else{
    binwidth = (max(to.plot$ECA_time_mean, na.rm = TRUE) - min(to.plot$ECA_time_mean, na.rm = TRUE))/20
    hist(to.plot[,ECA_time_mean], xlim = c(0, 1.05 * max(to.plot[,ECA_time_mean], na.rm = TRUE)),
         breaks = seq(0, max(to.plot[,ECA_time_mean], na.rm = TRUE)*1.05, binwidth), col = fill.multi, border = l.col, main = NA,
         xlab = NA, ylab = NA)
  }


  title(main = paste("SNV densities at ECA"), cex.main = 1)
  mtext(text = "SNVs per Mb", side = 1, line = 2, cex = 0.7)
  mtext(text = "No. of tumors", side = 2, line = 1.8, cex = 0.7)

  # Cumulative mutation densities at ECA:
  par(mar = c(3, 4, 3, 1), xpd = FALSE)

  x.min = 0
  x.max = max(lachesis$ECA_time_upper, na.rm = TRUE)*1.3
  y.min = 0
  y.max = 1
  plot(NA, NA, xlim=c(x.min, x.max), ylim=c(y.min, y.max), xlab = NA, ylab = NA, main = NA, axes = FALSE, frame.plot = FALSE)
  Axis(side=1, cex = 0.7)
  Axis(side=2, cex = 0.7)
  mtext(text = "SNVs per Mb", side = 1, line = 2, cex = 0.7)
  mtext(text = "No. of tumors", side = 2, line = 2, cex = 0.7)

  to.plot <- data.frame(x.lower = rep(sort(c( lachesis$ECA_time_mean)), each = 2)[-1],
                        x.upper = rep(sort(c( lachesis$ECA_time_mean)), each = 2)[-2*(nrow(lachesis[!is.na(ECA_time_mean),]) )])
  to.plot$y.lower <- sapply(rep(sort(c( lachesis$ECA_time_mean)), each = 2), function(x){sum(lachesis$ECA_time_upper <= x, na.rm = TRUE)})[-2*(nrow(lachesis[!is.na(ECA_time_mean),]) )]
  to.plot$y.upper <- sapply(rep(sort(c( lachesis$ECA_time_mean)), each = 2), function(x){sum(lachesis$ECA_time_lower <= x, na.rm = TRUE)})[-1]

  polygon(c(to.plot$x.lower, rev(to.plot$x.upper)),
          c(to.plot$y.lower, rev(to.plot$y.upper))/nrow(lachesis[!is.na(ECA_time_mean),]),
          col = fill.multi, border = NA)
  plot.ecdf(lachesis$ECA_time_mean, col = "black", add = TRUE, verticals = TRUE)

  title(main = paste("SNV densities at ECA"), cex.main = 1)

  if(!is.null(output.file)){
    dev.off()
  }

}


#' Correlate SNV density at ECA/MRCA with clinical parameters such as age, OS, etc.
#' @description
#' Takes SNV densities as computed by `LACHESIS` as input and correlates them with clinical data such as age at diagnosis, survival data etc.
#' @param lachesis output generated from \code{\link{LACHESIS}}
#' @param clin.par the clinical parameter used for correlation. Default `Age`.
#' @param suppress.outliers shall outliers (defined as the 2.5% tumors with lowest and highest densities) be plotted? Default `TRUE`.
#' @param log.densities plot logarithmic densities. Default `FALSE`
#' @param output.file optional; the file to which the plot will be stored.
#' @examples
#' #an example file with sample annotations and meta data
#' input.files = system.file("extdata", "Sample_template.txt", package = "LACHESIS")
#' input.files = data.table::fread(input.files)
#'
#' #cnv and snv files for example tumors
#' nbe11 = list.files(system.file("extdata/NBE11/", package = "LACHESIS"), full.names = TRUE)
#' nbe15 = list.files(system.file("extdata/NBE15/", package = "LACHESIS"), full.names = TRUE)
#' nbe63 = list.files(system.file("extdata/NBE63/", package = "LACHESIS"), full.names = TRUE)
#'
#' cnv.file = c(nbe11[1], nbe15[1], nbe63[1])
#' snv.file = c(nbe11[2], nbe15[2], nbe63[2])
#'
#' input.files$cnv.file = cnv.file
#' input.files$snv.file = snv.file
#'
#' # Make an axample input file with paths to cnv and snv file along with other meta data
#' lachesis_input = tempfile(pattern = "lachesis", tmpdir = tempdir(), fileext = ".tsv")
#' data.table::fwrite(x = input.files, file = lachesis_input, sep = "\t")
#'
#' #Exampele with template file with paths to multiple cnv/snv files as an input
#' lachesis <- LACHESIS(input.files = lachesis_input)
#' plotClinicalCorrelations(lachesis)
#' @export
#' @importFrom graphics abline Axis box grid hist mtext par rect text title arrows points
#' @importFrom stats cor

plotClinicalCorrelations <- function(lachesis = NULL, clin.par = "Age", suppress.outliers = FALSE, log.densities = FALSE,  output.file = NULL){

  ECA_time_mean <- NULL

  if(!clin.par %in% colnames(lachesis)){
    stop(clin.par, " not found!")
  }
  if(is.null(lachesis)){
    stop("Missing input. Please provide the output generated by LACHESIS()")
  }
  if(any(is.na(lachesis$MRCA_time_mean))){
    warning("Removing ", sum(is.na(lachesis$MRCA_time_mean)), " samples with missing MRCA density estimate.")
    lachesis <- lachesis[!is.na(MRCA_time_mean),]
  }
  if(nrow(lachesis)==0){
    warning("No sample with MRCA density estimate provided. Returning zero.")
    return(NULL)
  }
  if(!is.null(output.file)){
    pdf(output.file, width = 8, height = 6)
  }

  lo_mat <- matrix(data = c(1,2), nrow = 1, ncol=2, byrow = TRUE)
  graphics::layout(mat = lo_mat, widths = c(1, 1), heights = c(1, 1))

  to.plot <- lachesis

  if(suppress.outliers){
    to.plot <- to.plot[MRCA_time_mean < quantile(MRCA_time_mean, 0.975) &
                         MRCA_time_mean > quantile(MRCA_time_mean, 0.025),]
  }
  # Correlation between SNV density and the clinical parameters
  to.plot <- to.plot[!is.na(get(clin.par)),]
  if(nrow(to.plot)>0){
    par(mar = c(3, 4, 3, 1), xpd = FALSE)

    to.plot[,plot(MRCA_time_mean, get(clin.par), xlab = NA, ylab = NA, xlim = c(0, 1.05*max(MRCA_time_mean)),
                  ylim = c(0, 1.05*max(get(clin.par))), cex.axis = 0.7, log = ifelse(log.densities, "x", ""))]
    title(main = "SNV densities at MRCA vs age", cex.main = 1)
    mtext(text = "SNVs per Mb", side = 1, line = 2, cex = 0.8)
    mtext(text = clin.par, side = 2, line = 1.8, cex = 0.8)
    text(0.55*1.05*max(to.plot[,MRCA_time_mean]),
         0.75*max(to.plot[,get(clin.par)]), labels = paste("Pearson's r = ", to.plot[,round(cor(MRCA_time_mean, get(clin.par)), digits=3)]),
         cex = 0.8)

    par(mar = c(3, 4, 3, 1), xpd = FALSE)

    to.plot[,plot(ECA_time_mean, get(clin.par), xlab = NA, ylab = NA, xlim = c(0, 1.05*max(ECA_time_mean)),
                  ylim = c(0, 1.05*max(get(clin.par))), cex.axis = 0.7, log = ifelse(log.densities, "x", ""))]
    title(main = "SNV densities at ECA vs age", cex.main = 1)
    mtext(text = "SNVs per Mb", side = 1, line = 2, cex = 0.8)
    mtext(text = clin.par, side = 2, line = 1.8, cex = 0.8)
    text(0.5*1.05*max(to.plot[,ECA_time_mean]),
         0.75*max(to.plot[,get(clin.par)]), labels = paste("Pearson's r = ", to.plot[,round(cor(ECA_time_mean, get(clin.par)), digits=3)]),
         cex = 0.8)
  }else{
    # empty plots
    plot.new()
    plot.new()
  }

  if(!is.null(output.file)){
    dev.off()
  }

}

