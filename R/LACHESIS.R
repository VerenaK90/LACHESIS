#' Run MRCA density estimation for a set of tumors
#' @description
#' Takes a set of SNV and CNV files as input and outputs per-tumor SNV densities. Input can either be a tab-delimited file containing
#' the sample specifications or vectors giving direct paths to the sample files. CNV file requires columns for the chromosome number, start and end of the segment, and either the total copy number or the number of A- and B-alleles
#' @param input.files a tab-delimited sample-specification file, it must contain the sample name, the path to the SNV file, path to CNV file, and optionally purity, ploidy, cnv.chr.col, cnv.start.col, cnv.end.col, cnv.A.col, cnv.B.col, cnv.tcn.col. A template for this spreadsheet can be downloaded from ...
#' @param ids vector of sample names, will be ignored if `input.files` is specified.
#' @param vcf.tumor.ids vector of sample names as given in the vcf file; will be ignored if `input.files` is specified.
#' @param cnv.files vector of cnv files in same order as ids; should be in tab-delimited format, will be ignored if `input.files` is specified.
#' @param snv.files vector of snv files in same order as ids; should be in vcf format, will be ignored if `input.files` is specified.
#' @param vcf.source Tool used for generating VCF file. Can be `strelka` or `mutect` or `dkfz`.
#' @param purity vector tumor cell content in same order as ids; will be ignored if `input.files` is specified.
#' @param ploidy average copy number in the tumor sample in same order as ids; will be ignored if `input.files` is specified.
#' @param cnv.chr.col column index of chromosome number in cnv.files.
#' @param cnv.start.col column index of first position of the segment.
#' @param cnv.end.col column index of last position of the segment.
#' @param cnv.A.col column index of the number of A alleles. If A and B are not provided, allele configuration are assumed as 1:1 for disomic, 2:1 for trisomic and 3:1 for tetrasomic regions.
#' @param cnv.B.col column index of the number of B alleles. If A and B are not provided, allele configuration are assumed as 1:1 for disomic, 2:1 for trisomic and 3:1 for tetrasomic regions.
#' @param cnv.tcn.col column index of the total copy number. Is computed to A + B if not provided.
#' @param output.dir link to directory in which output is to be stored.
#' @param age, optional, the age at diagnosis.
#' @param OS.time, optional, overall survival time.
#' @param OS, optional, overall survival indicator variable.
#' @param EFS.time, optional, event-free survival time.
#' @param EFS, optional, event-free survival indicator variable.
#' @param min.cn minimum copy number to be included in the analysis. Default 2.
#' @param max.cn maximum copy number to be included in the analysis. Default 4.
#' @param merge.tolerance the maximum distance below which adjacent segments with equal copy number are merged. Defaults to 10^5 bp.
#' @param ignore.XY Ignore allosomes. Default TRUE.
#' @param min.vaf Remove variants with vcf below threshold. Default 0.01.
#' @param min.depth Minimum required depth for a variant to be considered. Default 30.
#' @param vcf.info.af The string encoding the allele frequency field in the FORMAT column of the .vcf file. Defaults to `AF`and will be ignored if `vcf.source` != `sentieon`.
#' @param vcf.info.dp The string encoding the read depth field in the FORMAT column of the .vcf file. Defaults to `DP`and will be ignored if `vcf.source` != `sentieon`.
#' @param min.seg.size the minimal segment length to be included in the quantification.
#' @param fp.mean optional, the average false positive rate of clonal mutations (e.g., due to incomplete tissue sampling). Defaults to 0.
#' @param fp.sd optional, the standard deviation of the false positive rate of clonal mutations (e.g., due to incomplete tissue sampling). Defaults to 0.
#' @param excl.chr a vector of chromosomes that should be excluded from the quantification. e.g., due to reporter constructs in animal models.
#' @param ref.build Reference genome. Default `hg19`. Can be `hg18`, `hg19` or `hg38`.
#' @param seed Integer. Can be user-specified or an automatically generated random seed, it will be documented in the log file.
#' @param filter.value The FILTER column value for variants that passed the filtering, defaults to PASS.
#' @param sig.assign Logical. If TRUE, each variant will be assigned to the most likely mutational signature.
#' @param assign.method Method to assign signatures: "max" to assign the signature with the highest probability, "sample" to randomly assign based on signature probabilities.
#' @param sig.file File path to the SigAssignment output file, typically named "Decomposed_MutationType_Probabilities.txt".
#' @param sig.select A character vector of specific signatures to include in the analysis (e.g., c("SBS1", "SBS5", "SBS40") to focus on clock-like mutational processes).
#' @param min.p Numeric. The minimum probability threshold from the SigAssignment output that a variant must meet to be considered as matching a specific signature.
#' @param ... further arguments and parameters passed to LACHESIS functions.
#' @examples
#' # An example file with sample annotations and meta data
#' input.files = system.file("extdata", "Sample_template.txt", package = "LACHESIS")
#' input.files = data.table::fread(input.files)
#'
#' # cnv and snv files for example tumors
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
#' # Example with template file with paths to multiple cnv/snv files as an input
#' lachesis <- LACHESIS(input.files = lachesis_input)
#'
#' # Example with a single sample input
#' strelka_vcf = system.file("extdata","strelka2.somatic.snvs.vcf.gz", package = "LACHESIS")
#' aceseq_cn = system.file("extdata", "ACESeq/NBE11_comb_pro_extra2.59_0.83.txt", package = "LACHESIS")
#' lachesis <- LACHESIS(ids = "NBE11", cnv.files = aceseq_cn, snv.files = strelka_vcf, vcf.source = "strelka", purity = 0.83, ploidy = 2.59, filter.value = "LowEVS")
#'
#' # Example with multiple sample and data frame input
#' nbe11_vcf = system.file("extdata","NBE11/snvs_NBE11_somatic_snvs_conf_8_to_10.vcf", package = "LACHESIS")
#' nbe11_cn = read.delim(system.file("extdata", "NBE11/NBE11_comb_pro_extra2.59_0.83.txt", package = "LACHESIS"), sep = "\t", header = TRUE)
#' nbe15_vcf = system.file("extdata","NBE15/snvs_NBE15_somatic_snvs_conf_8_to_10.vcf", package = "LACHESIS")
#' nbe15_cn = read.delim(system.file("extdata", "NBE15/NBE15_comb_pro_extra2.51_1.txt", package = "LACHESIS"), sep = "\t", header = TRUE)
#' lachesis <- LACHESIS(ids = c("NBE11", "NBE15"), cnv.files = list(nbe11_cn, nbe15_cn), snv.files = c(nbe11_vcf, nbe15_vcf), vcf.source = c("dkfz", "dkfz"), purity = c(0.83, 1), ploidy = c(2.59, 2.51), cnv.chr.col = c(1, 1), cnv.start.col = c(2, 2), cnv.end.col = c(3, 3), cnv.A.col = c(34, 34), cnv.B.col = c(35, 35), cnv.tcn.col = c(37, 37))
#'
#' @seealso \code{\link{MRCA}} \code{\link{clonalMutationCounter}} \code{\link{normalizeCounts}}
#' @import tidyr
#' @importFrom utils packageVersion
#' @return a data.table
#' @export

LACHESIS <- function(input.files = NULL, ids = NULL, vcf.tumor.ids = NULL, cnv.files = NULL, snv.files = NULL, vcf.source = NULL,
                     purity = NULL, ploidy = NULL,
                     cnv.chr.col = NULL, cnv.start.col = NULL, cnv.end.col = NULL, cnv.A.col = NULL,
                     cnv.B.col = NULL, cnv.tcn.col = NULL, age = NULL,
                     OS.time = NULL, OS = NULL, EFS.time = NULL, EFS = NULL, output.dir = NULL,
                     ignore.XY = TRUE, min.cn = 1, max.cn = 4, merge.tolerance = 10^5, min.vaf = 0.01, min.depth = 30,
                     vcf.info.af = "AF", vcf.info.dp = "DP", min.seg.size = 10^7, fp.mean = 0, fp.sd = 0, excl.chr = NULL,
                     ref.build = "hg19", seed = NULL, filter.value = "PASS", sig.assign = FALSE, sig.file = NULL, assign.method = "sample", sig.select = NULL, min.p = NULL, ...){


  ID <- cnv.file <- snv.file <- fwrite <- NULL

  if(is.null(input.files) & is.null(cnv.files)){
    stop("Missing input file!")
  }else if(is.null(input.files)){
    if(any(is.null(cnv.files), is.null(snv.files))){
      stop("Missing snv and cnv inputs!")
    }else if (!(is.data.frame(cnv.files) || is.data.table(cnv.files))) {
      if (length(cnv.files) != length(snv.files)) {
        stop("Please provide snv and cnv input for every sample!")
      }
    }
  }

  if (is.null(seed)) {
    seed <- sample.int(.Machine$integer.max, 1)
  }

  incl.chr <- setdiff(c(1:22), excl.chr)
  if(!ignore.XY){
    incl.chr <- c(incl.chr, "X", "Y")
  }

  # Collect ECA and MRCA densities for each tumor
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

  # Initializing datatable for logfile
  if(!is.null(input.files)){
    col.class <- fread(input.files)
    if (is.null(col.class$cnv.chr.col)) {
      col.class <- NA
    } else if (is.numeric(col.class$cnv.chr.col)) {
      col.class <- numeric()
    } else if (is.character(col.class$cnv.chr.col)) {
      col.class <- character()
    }

  log.file.data.cohort <- data.table::data.table(Sample_ID = character(),
                                                 package.version = character(),
                                                 vcf.tumor.ids = character(),
                                                 vcf.source = character(),
                                                 ploidy = numeric(),
                                                 purity = numeric(),
                                                 age = numeric(),
                                                 cnv.chr.col = col.class,
                                                 cnv.start.col = col.class,
                                                 cnv.end.col = col.class,
                                                 cnv.A.col = col.class,
                                                 cnv.B.col = col.class,
                                                 cnv.tcn.col = col.class,
                                                 OS.time = numeric(),
                                                 OS = numeric(),
                                                 EFS.time = numeric(),
                                                 EFS = numeric(),
                                                 output.dir = character(),
                                                 ignore.XY = logical(),
                                                 min.cn = numeric(),
                                                 max.cn = numeric(),
                                                 merge.tolerance = numeric(),
                                                 min.vaf = numeric(),
                                                 min.depth = numeric(),
                                                 vcf.info.af = numeric(),
                                                 vcf.info.dp = numeric(),
                                                 min.seg.size = numeric(),
                                                 fp.mean = double(),
                                                 fp.sd = double(),
                                                 excl.chr = numeric(),
                                                 ref.build = character(),
                                                 cnv.file = character(),
                                                 snv.file = character(),
                                                 seed = numeric())

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
      #x[,which(sapply(x, is.na)):=NULL] # remove NA entries

      if(is.null(x$ID)){
        stop("Please provide sample identifiers.")
      }
      if(is.null(x$vcf.source)){
        stop("Please provide vcf source.")
      }
      if(is.null(x$vcf.tumor.ids)){
        x$vcf.tumor.ids <- x$ID
      }else if(any(is.na(x$vcf.tumor.ids))){
        warning("No column identifier provided for sample ", which(is.na(x$vcf.tumor.ids)), "; will be inferred.")
        x$vcf.tumor.ids[is.na(x$vcf.tumor.ids)] <- x$id[is.na(x$vcf.tumor.ids)]
      }
      message("Computing SNV density for sample ", x$ID)

      if(!is.null(output.dir)){
        dir.create(paste(output.dir, x$ID, sep="/"), recursive = TRUE, showWarnings = FALSE) # Create per-sample output directory
      }else{
        warning("No output directory specified. LACHESIS output will be discarded.")
      }

      cnv <- readCNV(cn.info = x$cnv.file, chr.col = x$cnv.chr.col, start.col = x$cnv.start.col,
                     end.col = x$cnv.end.col, A.col = x$cnv.A.col, B.col = x$cnv.B.col,
                     tcn.col = x$cnv.tcn.col, tumor.id = x$ID, merge.tolerance = merge.tolerance,
                     max.cn = max.cn, ignore.XY = ignore.XY)

      snv <- readVCF(vcf = x$snv.file, vcf.source = x$vcf.source, t.sample = x$vcf.tumor.id, min.depth = min.depth,
                     min.vaf = min.vaf, info.af = vcf.info.af, info.dp = vcf.info.dp, filter.value = filter.value)

      nb <- nbImport(cnv = cnv, snv = snv, purity = x$purity, ploidy = x$ploidy, sig.assign = sig.assign, assign.method = assign.method, ID = x$ID, sig.file = sig.file, sig.select = sig.select, min.p = min.p, ref.build = ref.build, seed = seed)

      if(nrow(nb)==0){
        warning("Insufficient data for sample ", x$ID)
        this.tumor.density <- data.table::data.table(Sample_ID = x$ID,
                                                     MRCA_time_mean = NA,
                                                     MRCA_time_lower = NA,
                                                     MRCA_time_upper = NA,
                                                     ECA_time_mean = NA,
                                                     ECA_time_lower = NA,
                                                     ECA_time_upper = NA,
                                                     Ploidy = x$ploidy,
                                                     Purity = x$purity,
                                                     Age = x$Age,
                                                     OS.time = x$OS.time,
                                                     OS = x$OS,
                                                     EFS.time = x$EFS.time,
                                                     EFS = x$EFS)
        next
      }
      if(!is.null(output.dir)){
        plotVAFdistr(snv, output.file = paste(output.dir, x$ID, "VAF_histogram.pdf", sep="/"), ...)
        plotNB(nb = nb, samp.name = x$ID, output.file = paste(output.dir, x$ID, "VAF_histogram_strat.pdf", sep="/"), ref.build = ref.build, sig.output.file = paste(output.dir, x$ID, "VAF_histogram_strat_sig.pdf", sep="/"), ...)
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

      cohort.densities <- merge(cohort.densities, this.tumor.density, all=TRUE)

      # Output the result for this sample
      if(!is.null(output.dir)){
        mrca.densities <- transpose(data.table(unlist(attributes(mrca)[c("purity", "ploidy", "MRCA_time_mean", "MRCA_time_lower", "MRCA_time_upper", "ECA_time_mean", "ECA_time_lower", "ECA_time_upper")]))
        )
        data.table::fwrite(mrca.densities, file = file.path(output.dir, x$ID, paste0("MRCA_densities_", x$ID, ".txt")), quote = FALSE, col.names = FALSE, sep="\t")
        data.table::fwrite(mrca, file = file.path(output.dir, x$ID, paste0("SNV_timing_per_segment_", x$ID, ".txt")), row.names = FALSE, quote = FALSE, sep = "\t")
      }

      if(!is.null(output.dir)){
        plotMutationDensities(mrcaObj = mrca, samp.name = x$ID, output.file = paste(output.dir, x$ID, "SNV_densities.pdf", sep="/"), ...)
      }

      # Collecting data for log file
      package.version <- as.character(utils::packageVersion("LACHESIS"))
      log.file.data.single <- data.table::data.table(Sample_ID = x$ID,
                                                      package.version = package.version,
                                                      vcf.tumor.ids = x$vcf.tumor.ids,
                                                      vcf.source = x$vcf.source,
                                                      ploidy = x$ploidy,
                                                      purity = x$purity,
                                                      cnv.chr.col = x$cnv.chr.col,
                                                      cnv.start.col = x$cnv.start.col,
                                                      cnv.end.col = x$cnv.end.col,
                                                      cnv.A.col = x$cnv.A.col,
                                                      cnv.B.col = x$cnv.B.col,
                                                      cnv.tcn.col = x$cnv.tcn.col,
                                                      age = x$Age,
                                                      OS.time = x$OS.time,
                                                      OS = x$OS,
                                                      EFS.time = x$EFS.time,
                                                      EFS = x$EFS,
                                                      output.dir = output.dir,
                                                      ignore.XY = ignore.XY,
                                                      min.cn = min.cn,
                                                      max.cn = max.cn,
                                                      merge.tolerance = merge.tolerance,
                                                      min.vaf = min.vaf,
                                                      min.depth = min.depth,
                                                      vcf.info.af = vcf.info.af,
                                                      vcf.info.dp = vcf.info.dp,
                                                      min.seg.size = min.seg.size,
                                                      fp.mean = fp.mean,
                                                      fp.sd = fp.sd,
                                                      excl.chr = excl.chr,
                                                      ref.build = ref.build,
                                                      cnv.file = x$cnv.file,
                                                      snv.file = x$snv.file,
                                                      seed = seed)

      log.file.data.cohort <- merge(log.file.data.cohort, log.file.data.single, all=TRUE)

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
        dir.create(paste(output.dir, ids[i], sep="/"), recursive = TRUE, showWarnings = FALSE) # Create per-sample output directory
      }else{
        warning("No output directory specified. LACHESIS output will be discarded.")
      }

      if(is.na(cnv.files)[i]){
        warning("No CNV file provided for sample ", ids[i], "; sample will be excluded")
        next
      }
      if(is.na(snv.files)[i]){
        warning("No SNV file provided for sample ", ids[i], "; sample will be excluded")
        next
      }

      cnv <- readCNV(cn.info = cnv.files[[i]], chr.col = cnv.chr.col[i], start.col = cnv.start.col[i],
                     end.col = cnv.end.col[i], A.col = cnv.A.col[i], B.col = cnv.B.col[i],
                     tcn.col = cnv.tcn.col[i], tumor.id = ids[i], merge.tolerance = merge.tolerance,
                     max.cn = max.cn, ignore.XY = ignore.XY)

      snv <- readVCF(vcf = snv.files[i], vcf.source = vcf.source[i], t.sample = vcf.tumor.ids[i], min.depth = min.depth,
                     min.vaf = min.vaf, info.af = vcf.info.af, info.dp = vcf.info.dp, filter.value = filter.value)

      nb <- nbImport(cnv = cnv, snv = snv, purity = purity[i], ploidy = ploidy[i], sig.assign = sig.assign, assign.method = assign.method, ID = ids[i], sig.file = sig.file, sig.select = sig.select, min.p = min.p, ref.build = ref.build, seed = seed)

      if(nrow(nb)==0){
        warning("Insufficient data for sample ", x$ID)
        this.tumor.density <- data.table::data.table(Sample_ID = ids[i],
                                                     MRCA_time_mean = NA,
                                                     MRCA_time_lower = NA,
                                                     MRCA_time_upper = NA,
                                                     ECA_time_mean = NA,
                                                     ECA_time_lower = NA,
                                                     ECA_time_upper = NA,
                                                     Ploidy = ploidy[i],
                                                     Purity = purity[i],
                                                     Age = age[i],
                                                     OS.time = OS.time[i],
                                                     OS = OS[i],
                                                     EFS.time = EFS.time[i],
                                                     EFS = EFS[i])

        next
      }
      if(!is.null(output.dir)){
        plotVAFdistr(snv, output.file = paste(output.dir, ids[i], "VAF_histogram.pdf", sep="/"), ...)
        plotNB(nb = nb, samp.name = ids[i], output.file = paste(output.dir, ids[i], "VAF_histogram_strat.pdf", sep="/"), ref.build = ref.build, sig.output.file = paste(output.dir, x$ID, "VAF_histogram_strat_sig.pdf", sep="/"), ...)
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

        # Output the result for this sample
        if(!is.null(output.dir)){
          mrca.densities <- transpose(data.table(unlist(attributes(mrca)[c("purity", "ploidy", "MRCA_time_mean", "MRCA_time_lower", "MRCA_time_upper", "ECA_time_mean", "ECA_time_lower", "ECA_time_upper")]))
          )
          data.table::fwrite(mrca.densities, file = file.path(output.dir, ids[i], paste0("MRCA_densities_", ids[i], ".txt")), quote = F, col.names = F, sep="\t")
          data.table::fwrite(mrca, file = file.path(output.dir, ids[i], paste0("SNV_timing_per_segment_", ids[i], ".txt")), row.names = FALSE, quote = FALSE, sep = "\t")
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

  # Plot the distribution of Mutation densities at ECA and MRCA

  if(!is.null(output.dir)){
    plotLachesis(cohort.densities, output.file = paste(output.dir, "SNV_densities_cohort.pdf", sep="/"), ...)
  }

  # Save log file as tsv
  if (!is.null(output.dir) && !is.null(input.files)) {
    timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
    output.file <- paste0(output.dir, "/LACHESIS_logfile_", timestamp, ".tsv")
    fwrite(log.file.data.cohort, output.file, sep = "\t")
  }

  return(cohort.densities)
}


#' Plot SNV densities at ECA and MRCA
#' @description
#' Visualizes results from \code{\link{LACHESIS}}. Top plot, histograms of mean mutation densities; bottom plots, cumulative distribution of mean mutation densities with 95% confidence intervals.
#' @param lachesis output generated from \code{\link{LACHESIS}}
#' @param lach.suppress.outliers whether outliers (defined as the 2.5% tumors with lowest and highest densities) are to be plot. Default `TRUE`.
#' @param lach.log.densities plot logarithmic densities. Default `FALSE`.
#' @param lach.col.zero optional, bar color for single-copy SSNV densities.
#' @param lach.col.multi optional, bar color for multi-copy SSNV densities.
#' @param lach.border, optional, border color for the bars.
#' @param binwidth optional, the binwidth in the histogram.
#' @param output.file optional, the file to which the plot will be stored.
#' @param ... further arguments and parameters passed to other LACHESIS functions.
#' @examples
#' # An example file with sample annotations and meta data
#' input.files = system.file("extdata", "Sample_template.txt", package = "LACHESIS")
#' input.files = data.table::fread(input.files)
#'
#' # cnv and snv files for example tumors
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
#' # Example with template file with paths to multiple cnv/snv files as an input
#' lachesis <- LACHESIS(input.files = lachesis_input)
#' plotLachesis(lachesis)
#' @export
#' @importFrom graphics abline Axis box grid hist mtext par rect text title arrows legend points polygon

plotLachesis <- function(lachesis = NULL, lach.suppress.outliers = FALSE, lach.log.densities = FALSE, lach.col.multi = "#176A02", lach.border = NULL, binwidth = NULL, lach.col.zero = "#4FB12B", output.file = NULL, ...){

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

  # I. Plot histograms
  to.plot <- lachesis

  if(lach.suppress.outliers){
    to.plot <- to.plot[MRCA_time_mean < quantile(MRCA_time_mean, 0.975) &
              MRCA_time_mean > quantile(MRCA_time_mean, 0.025),]
  }

  if(is.null(binwidth)){
    if(lach.log.densities){
      binwidth = (max(log10(to.plot$MRCA_time_mean)) - min(log10(to.plot$MRCA_time_mean)))/20
    }else{
      binwidth = (max(to.plot$MRCA_time_mean) - min(to.plot$MRCA_time_mean))/20
    }

  }

  lo_mat <- matrix(data = c(1, 2, 3, 4), nrow = 2, ncol=2, byrow = TRUE)
  graphics::layout(mat = lo_mat, widths = c(1, 2, 1, 2), heights = c(1, 1, 1, 1))

  par(mar = c(3, 4, 3, 1))

  if(lach.log.densities){

    min.x <- floor(min(to.plot[,log10(MRCA_time_mean)]))
    max.x <- ceiling(max(to.plot[,log10(MRCA_time_mean)]))

    hist(to.plot[,log10(MRCA_time_mean)], xlim = c(min.x, max.x),
         breaks = 20,
         col = lach.col.zero, border = lach.border, main = NA,
         xlab = NA, ylab = NA, axes = FALSE)

    Axis(side = 1, at = seq(min.x, max.x, length.out = 10),
         labels = round(10^seq(min.x, max.x, length.out = 10), digits = 2))
    Axis(side = 2)

  }else{
    hist(to.plot[,MRCA_time_mean], xlim = c(0, 1.05 * max(to.plot[,MRCA_time_mean])),
         breaks = seq(0, max(to.plot[,MRCA_time_mean])*1.05, binwidth), col = lach.col.zero, border = lach.border, main = NA,
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
  mtext(text = "Fraction of tumors", side = 2, line = 2, cex = 0.7)

  to.plot <- data.frame(x.lower = rep(sort(c( lachesis$MRCA_time_mean)), each = 2)[-1],
                        x.upper = rep(sort(c( lachesis$MRCA_time_mean)), each = 2)[-2*(nrow(lachesis) )])
  to.plot$y.lower <- sapply(rep(sort(c( lachesis$MRCA_time_mean)), each = 2), function(x){sum(lachesis$MRCA_time_upper <= x)})[-2*(nrow(lachesis) )]
  to.plot$y.upper <- sapply(rep(sort(c( lachesis$MRCA_time_mean)), each = 2), function(x){sum(lachesis$MRCA_time_lower <= x)})[-1]

  polygon(c(to.plot$x.lower, rev(to.plot$x.upper)),
          c(to.plot$y.lower, rev(to.plot$y.upper))/nrow(lachesis),
          col = lach.col.zero, border = NA)

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

  if(lach.suppress.outliers){
    to.plot <- to.plot[MRCA_time_mean < quantile(MRCA_time_mean, 0.975) &
                         MRCA_time_mean > quantile(MRCA_time_mean, 0.025),]
  }

  par(mar = c(3, 4, 3, 1))

  if(lach.log.densities){
    min.x <- floor(min(to.plot[,log10(ECA_time_mean)], na.rm = TRUE))
    max.x <- ceiling(max(to.plot[,log10(ECA_time_mean)], na.rm = TRUE))

    hist(to.plot[,log10(ECA_time_mean)], xlim = c(min.x, max.x),
         breaks = 20,
         col = lach.col.zero, border = lach.border, main = NA,
         xlab = NA, ylab = NA, axes = FALSE)

    Axis(side = 1, at = seq(min.x, max.x, length.out = 10),
         labels = round(10^seq(min.x, max.x, length.out = 10), digits = 2))
    Axis(side = 2)
  }else{
    max_ECA_time_mean <- max(to.plot$ECA_time_mean, na.rm = TRUE)
    min_ECA_time_mean <- min(to.plot$ECA_time_mean, na.rm = TRUE)

    if (max_ECA_time_mean != min_ECA_time_mean) {
      binwidth <- (max_ECA_time_mean - min_ECA_time_mean) / 20
    } else {
      binwidth <- max_ECA_time_mean / 20
      if (binwidth == 0 || is.na(binwidth)) {
        binwidth <- 0.01
      }
    }

    hist(to.plot[,ECA_time_mean], xlim = c(0, 1.05 * max(to.plot[,ECA_time_mean], na.rm = TRUE)),
         breaks = seq(0, max(to.plot[,ECA_time_mean], na.rm = TRUE)*1.05, binwidth), col = lach.col.multi, border = lach.border, main = NA,
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
  mtext(text = "Fraction of tumors", side = 2, line = 2, cex = 0.7)

  to.plot <- data.frame(x.lower = rep(sort(c( lachesis$ECA_time_mean)), each = 2)[-1],
                        x.upper = rep(sort(c( lachesis$ECA_time_mean)), each = 2)[-2*(nrow(lachesis[!is.na(ECA_time_mean),]) )])
  to.plot$y.lower <- sapply(rep(sort(c( lachesis$ECA_time_mean)), each = 2), function(x){sum(lachesis$ECA_time_upper <= x, na.rm = TRUE)})[-2*(nrow(lachesis[!is.na(ECA_time_mean),]) )]
  to.plot$y.upper <- sapply(rep(sort(c( lachesis$ECA_time_mean)), each = 2), function(x){sum(lachesis$ECA_time_lower <= x, na.rm = TRUE)})[-1]

  polygon(c(to.plot$x.lower, rev(to.plot$x.upper)),
          c(to.plot$y.lower, rev(to.plot$y.upper))/nrow(lachesis[!is.na(ECA_time_mean),]),
          col = lach.col.multi, border = NA)
  plot.ecdf(lachesis$ECA_time_mean, col = "black", add = TRUE, verticals = TRUE)

  title(main = paste("SNV densities at ECA"), cex.main = 1)

  if(!is.null(output.file)){
    dev.off()
  }

}


#' Correlate SNV density at ECA/MRCA with clinical parameters such as age, OS, etc.
#' @description
#' Takes SNV densities as computed by `LACHESIS` as input and correlates them with clinical data such as age at diagnosis, survival data etc.
#' @param lachesis output generated from \code{\link{LACHESIS}}.
#' @param clin.par the clinical parameter used for correlation. Default `Age`.
#' @param clin.suppress.outliers shall outliers (defined as the 2.5% tumors with lowest and highest densities) be plotted? Default `TRUE`.
#' @param clin.log.densities plot logarithmic densities. Default `FALSE`.
#' @param output.file optional; the file to which the plot will be stored.
#' @examples
#' # An example file with sample annotations and meta data
#' input.files = system.file("extdata", "Sample_template.txt", package = "LACHESIS")
#' input.files = data.table::fread(input.files)
#'
#' # cnv and snv files for example tumors
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
#' # Example with template file with paths to multiple cnv/snv files as an input
#' lachesis <- LACHESIS(input.files = lachesis_input)
#' plotClinicalCorrelations(lachesis)
#' @export
#' @importFrom graphics abline Axis box grid hist mtext par rect text title arrows points
#' @importFrom stats cor

plotClinicalCorrelations <- function(lachesis = NULL, clin.par = "Age", clin.suppress.outliers = FALSE, clin.log.densities = FALSE,  output.file = NULL){

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

  if(clin.suppress.outliers){
    to.plot <- to.plot[MRCA_time_mean < quantile(MRCA_time_mean, 0.975) &
                         MRCA_time_mean > quantile(MRCA_time_mean, 0.025),]
  }
  # Correlation between SNV density and the clinical parameters
  to.plot <- to.plot[!is.na(get(clin.par)),]
  if(nrow(to.plot)>0){
    par(mar = c(3, 4, 3, 1), xpd = FALSE)

    to.plot[,plot(MRCA_time_mean, get(clin.par), xlab = NA, ylab = NA, xlim = c(0, 1.05*max(MRCA_time_mean)),
                  ylim = c(0, 1.05*max(get(clin.par))), cex.axis = 0.7, log = ifelse(clin.log.densities, "x", ""))]
    title(main = "SNV densities at MRCA vs age", cex.main = 1)
    mtext(text = "SNVs per Mb", side = 1, line = 2, cex = 0.8)
    mtext(text = clin.par, side = 2, line = 1.8, cex = 0.8)
    text(0.55*1.05*max(to.plot[,MRCA_time_mean]),
         0.75*max(to.plot[,get(clin.par)]), labels = paste("Pearson's r = ", to.plot[,round(cor(MRCA_time_mean, get(clin.par)), digits=3)]),
         cex = 0.8)

    par(mar = c(3, 4, 3, 1), xpd = FALSE)

    to.plot[,plot(ECA_time_mean, get(clin.par), xlab = NA, ylab = NA, xlim = c(0, 1.05*max(ECA_time_mean, na.rm=TRUE)),
                  ylim = c(0, 1.05*max(get(clin.par))), cex.axis = 0.7, log = ifelse(clin.log.densities, "x", ""))]

    title(main = "SNV densities at ECA vs age", cex.main = 1)
    mtext(text = "SNVs per Mb", side = 1, line = 2, cex = 0.8)
    mtext(text = clin.par, side = 2, line = 1.8, cex = 0.8)
    text(0.5*1.05*max(to.plot[,ECA_time_mean]),
         0.75*max(to.plot[,get(clin.par)]), labels = paste("Pearson's r = ", to.plot[,round(cor(ECA_time_mean, get(clin.par)), digits=3)]),
         cex = 0.8)
  }else{
    # Empty plots
    plot.new()
    plot.new()
  }

  if(!is.null(output.file)){
    dev.off()
  }

}

#' Correlate SNV density timing at MRCA with Survival
#' @description
#' Takes SNV density timing as computed by `LACHESIS` as input and compares survival between tumors with high and low SNV densities
#' @param lachesis output generated from \code{\link{LACHESIS}}.
#' @param mrca.cutpoint optional, MRCA density value to be used for survival stratification, will be computationally inferred to maximize survival differences if not specified by user.
#' @param output.dir the directory to which the plot will be stored.
#' @param surv.time column name containing survival time; defaults to `OS.time`.
#' @param surv.event column name containing event; defaults to `OS`.
#' @param surv.time.breaks numeric value controlling time axis breaks; defaults to `NULL`.
#' @param surv.time.scale numeric value by which survival time is to be divided (e.g., 365 for converting days into years, 30 for months), defaults to `1`.
#' @param surv.palette color palette to be used. Allowed values include "hue" for the default hue color scale; "grey" for grey color palettes; brewer palettes e.g. "RdBu", "Blues", ...; or custom color palette e.g. c("blue", "red").
#' @param surv.title main title.
#' @param surv.ylab y-axis label, defaults to `Survival`.
#' @param output.dir link to directory in which output is to be stored.
#' @examples
#' # An example file with sample annotations and meta data
#' input.files = system.file("extdata", "Sample_template.txt", package = "LACHESIS")
#' input.files = data.table::fread(input.files)
#'
#' # cnv and snv files for example tumors
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
#' # Example with template file with paths to multiple cnv/snv files as an input
#' lachesis <- LACHESIS(input.files = lachesis_input)
#' plotSurvival(lachesis, surv.time = 'EFS.time', surv.event = 'EFS')
#'
#' @export
#' @import ggplot2
#' @import survival
#' @import survminer
#' @import gridExtra
#' @importFrom stats pchisq

plotSurvival <- function(lachesis = NULL, mrca.cutpoint = NULL, output.dir = NULL, surv.time = 'OS.time', surv.event = 'OS', surv.palette = c("dodgerblue", "dodgerblue4"), surv.time.breaks = NULL, surv.time.scale = 1, surv.title = "Survival probability", surv.ylab = "Survival"){

  if (is.null(lachesis)) {
    stop("Error: 'lachesis' dataset must be provided.")
  }

  if(any(is.na(lachesis$MRCA_time_mean))){
    warning("Removing ", sum(is.na(lachesis$MRCA_time_mean)), " samples with missing MRCA density estimate.")
    lachesis <- lachesis[!is.na(MRCA_time_mean),]
  }

  if(!surv.time %in% colnames(lachesis)){
    stop("Error: please provide a valid column name for `surv.time`.")
  }

  if(!surv.event %in% colnames(lachesis)){
    stop("Error: please provide a valid column name for `surv.event`.")
  }

  if(any(is.na(lachesis[,..surv.time]))){
    warning("Removing ", sum(is.na(lachesis[,surv.time])), " samples with missing survival time.")
    lachesis <- lachesis[!is.na(get(surv.time)), .SD]
  }

  if(any(is.na(lachesis[,..surv.event]))){
    warning("Removing ", sum(is.na(lachesis[,surv.event])), " samples with missing survival event.")
    lachesis <- lachesis[!is.na(get(surv.event)), .SD]
  }

  if(nrow(lachesis)==0){
    warning("No sample with MRCA density estimate provided. Returning zero.")
    return(NULL)
  }

  if(all(lachesis[,..surv.event]==0)){
    warning("No survival events in cohort Returning zero.")
    return(NULL)
  }

  # Calculating MRCA cutpoint
  if(is.null(mrca.cutpoint)){
    mrca.cutpoint.obj <- survminer::surv_cutpoint(
      lachesis,
      time = surv.time,
      event = surv.event,
      variables = c("MRCA_time_mean")
    )

    mrca.cutpoint <- as.numeric(mrca.cutpoint.obj$cutpoint["MRCA_time_mean", "cutpoint"])
  }

  # Categorizing according to MRCA
  lachesis.categorized <- lachesis
  lachesis.categorized[[surv.time]] <- as.numeric(lachesis.categorized[[surv.time]])/surv.time.scale
  lachesis.categorized$MRCA_timing <- ifelse(lachesis.categorized$MRCA_time_mean < mrca.cutpoint, "early", "late")
  lachesis.categorized$MRCA_timing <- factor(lachesis.categorized$MRCA_timing, levels=c("early", "late"))

  # Survival analysis
  survival.fit <- survival::survfit(Surv(time = unlist(lachesis.categorized[,..surv.time]), event = unlist(lachesis[,..surv.event])) ~ MRCA_timing,
                                         data = lachesis.categorized)
  survival.diff <- survival::survdiff(Surv(time = unlist(lachesis.categorized[,..surv.time]), event = unlist(lachesis[,..surv.event])) ~ MRCA_timing,
                                           data = lachesis.categorized)

  p_value <- 1 - stats::pchisq(survival.diff$chisq, length(survival.diff$n) - 1)
  p.value.pos <- max(survival.fit$time) * (1/6)

  survival.fit.plot <- survminer::ggsurvplot_df(surv_summary(survival.fit, data = lachesis.categorized), title = surv.title, conf.int = TRUE, color = "strata", censor.shape = 124,
                                                palette = surv.palette, xlab = "Time", ylab = surv.ylab, legend.labs = c("Early MRCA", "Late MRCA"), break.time.by = surv.time.breaks) +
    annotate("text", x = p.value.pos, y = 0.2, label = paste0("p = ", ifelse(p_value > 0 & p_value < 0.0001, "< 0.0001", formatC(p_value, format = "f", digits = 4))), size = 5)

  survival.fit.risk.table <- survminer::ggrisktable(survival.fit, data = lachesis.categorized, legend.labs = c("Early MRCA", "Late MRCA"), break.time.by = surv.time.breaks)

  # Printing cutpoint txt
  if (!is.null(output.dir)) {
    mrca.cutpoint.dt <- data.table::data.table(
      cutpoint = mrca.cutpoint,
      statistic = as.numeric(mrca.cutpoint.obj$cutpoint["MRCA_time_mean", "statistic"]),
      p_value = p_value
    )

    mrca.cutpoint.rounded <- formatC(mrca.cutpoint, format = "f", digits = 2)
    data.table::fwrite(mrca.cutpoint.dt, file = file.path(output.dir, paste0("cutpoint_estimate_", mrca.cutpoint.rounded, ".txt")), sep = "\t")
  }

  # Printing pdf
  if(!is.null(output.dir)){
    pdf(paste0(output.dir, "/Stratified_", surv.event, ".pdf"), width = 9, height = 8)
  }
  gridExtra::grid.arrange(survival.fit.plot, survival.fit.risk.table,
                          ncol = 1, nrow = 2,
                          widths = c(1),
                          heights = c(3, 1))
  if(!is.null(output.dir)){
    dev.off()
  }
}


#' Classify a tumor's start of clonal outgrowth during tumorigenesis as "early" or "late" (favorable/ unfavorable prognosis) depending on the mutation density at its MRCA
#' @description
#' Takes SNV density timing as computed by `LACHESIS` as input and classifies the tumors in the cohort.
#' @param lachesis output generated from \code{\link{LACHESIS}}.
#' @param mrca.cutpoint optional; value based on SNV_densities_cohort.pdf observation, will be used as inferred from a test data set if not specified by user.
#' @param infer.cutpoint logical; should the MRCA cutpoint be inferred from the data?
#' @param entity optional; the tumor entity if classifying according to a pre-defined threshold. Currently, only "neuroblastoma" is supported.
#' @param surv.time column name containing survival time; defaults to `OS.time`.
#' @param surv.event column name containing event; defaults to `OS`.
#' @param output.dir link to directory in which output is to be stored.
#' @examples
#' # An example file with sample annotations and meta data
#' input.files = system.file("extdata", "Sample_template.txt", package = "LACHESIS")
#' input.files = data.table::fread(input.files)
#'
#' # cnv and snv files for example tumors
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
#' # Example with template file with paths to multiple cnv/snv files as an input
#' lachesis <- LACHESIS(input.files = lachesis_input)
#' classifyLACHESIS(lachesis)
#'
#' @export
#' @import survminer

classifyLACHESIS <- function(lachesis, mrca.cutpoint = NULL, output.dir = NULL, infer.cutpoint = FALSE, entity = "neuroblastoma", surv.time = 'OS.time', surv.event = 'OS'){

  if (is.null(lachesis)) {
    stop("Error: 'lachesis' dataset must be provided.")
  }
  entities <- c("neuroblastoma")
  entity <- match.arg(arg = entity, choices = entities, several.ok = FALSE)

  if(infer.cutpoint == TRUE & sum(!(is.na(lachesis[,..surv.time]))) < 2){
    stop("Please provide survival time if inferring cutpoint de novo.")
  }

  if(infer.cutpoint == TRUE & (sum(!(is.na(lachesis[,..surv.event]))) < 2 | sum(lachesis[,..surv.event]!=0, na.rm = T) < 2)){
    stop("Please provide survival information if inferring cutpoint de novo.")
  }
  message("Classifying ", entity, " samples.")

  if( infer.cutpoint==TRUE){
    message("MRCA cutpoint will be newly inferred.")
  }else if(is.null(mrca.cutpoint)){
    message("Samples will be classified according to established MRCA cutpoint for ", entity, ".")
  }else if(infer.cutpoint==FALSE & is.null(mrca.cutpoint)){
    message("Please provide cutpoint or set `infer.cutpoint`=`TRUE`")
  }else if(infer.cutpoint==FALSE){
    message("MRCA cutpoint taken as ", mrca.cutpoint, ".")
  }

  # Calculating MRCA cutpoint
  if(infer.cutpoint){
    mrca.cutpoint <- survminer::surv_cutpoint(
      lachesis[!is.na(get(surv.time)), .SD],
      time = surv.time,
      event = surv.event,
      variables = c("MRCA_time_mean")
    )
    mrca.cutpoint <- as.numeric(mrca.cutpoint$cutpoint["MRCA_time_mean", "cutpoint"])
  }else if(is.null(mrca.cutpoint)){
    mrca.cutpoint <- .getCutpoint(entity)
  }

  # Categorizing according to MRCA
  lachesis.categorized <- lachesis
  lachesis.categorized$MRCA_timing <- ifelse(lachesis.categorized$MRCA_time_mean < mrca.cutpoint, "early", "late")
  lachesis.categorized$MRCA_timing <- factor(lachesis.categorized$MRCA_timing, levels=c("early", "late"))

  attr(lachesis.categorized, "MRCA Cutpoint") <- mrca.cutpoint
  attr(lachesis.categorized, "Entity") <- entity

  if (!is.null(output.dir)) {
    data.table::fwrite(lachesis.categorized, file = file.path(output.dir, "Lachesis_classifier.txt"), sep = "\t")
  }

  return(lachesis.categorized)
}

# Cutpoint for neuroblastoma
.getCutpoint <- function(entity = "neuroblastoma"){

  if(entity == 'neuroblastoma'){
    cut.point = 0.05
  }else{
    stop('Available entities: neuroblastoma')
  }

  cut.point

}
