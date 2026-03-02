#' Run MRCA density estimation for a set of tumors
#' @description
#' Takes a set of SNV and CNV files as input and outputs per-tumor SNV
#' densities. Input can either be a tab-delimited file containing
#' the sample specifications or vectors giving direct paths to the sample files.
#'  CNV file requires columns for the chromosome number, start and end of the
#'  segment, and either the total copy number or the number of A- and B-alleles
#' @param input.files a tab-delimited sample-specification file, it must contain
#'  the sample name, the path to the SNV file, path to CNV file, and optionally
#'  purity, ploidy, cnv.chr.col, cnv.start.col, cnv.end.col, cnv.A.col,
#'  cnv.B.col, cnv.tcn.col. A template for this spreadsheet can be downloaded
#'  from ...
#' @param ids vector of sample names, will be ignored if `input.files` is
#' specified.
#' @param vcf.tumor.ids vector of sample names as given in the vcf file; will be
#'  ignored if `input.files` is specified.
#' @param cnv.files vector of cnv files in same order as ids; should be in
#' tab-delimited format, will be ignored if `input.files` is specified.
#' @param snv.files vector of snv files in same order as ids; should be in vcf
#' format, will be ignored if `input.files` is specified.
#' @param vcf.source Tool used for generating VCF file. Can be `strelka` or
#' `mutect` or `dkfz` or `sentieon`.
#' @param purity vector tumor cell content in same order as ids; will be ignored
#'  if `input.files` is specified.
#' @param ploidy average copy number in the tumor sample in same order as ids;
#' will be ignored if `input.files` is specified.
#' @param cnv.chr.col column index of chromosome number in cnv.files.
#' @param cnv.start.col column index of first position of the segment.
#' @param cnv.end.col column index of last position of the segment.
#' @param cnv.A.col column index of the number of A alleles. If A and B are not
#' provided, allele configuration are assumed as 1:1 for disomic, 2:1 for
#' trisomic and 3:1 for tetrasomic regions.
#' @param cnv.B.col column index of the number of B alleles. If A and B are not
#' provided, allele configuration are assumed as 1:1 for disomic, 2:1 for
#' trisomic and 3:1 for tetrasomic regions.
#' @param cnv.tcn.col column index of the total copy number. Is computed to A +
#' B if not provided.
#' @param output.dir link to directory in which output is to be stored.
#' @param age, optional, the age at diagnosis.
#' @param OS.time, optional, overall survival time.
#' @param OS, optional, overall survival indicator variable.
#' @param EFS.time, optional, event-free survival time.
#' @param EFS, optional, event-free survival indicator variable.
#' @param min.cn minimum copy number to be included in the analysis. Default 2.
#' @param max.cn maximum copy number to be included in the analysis. Default 4.
#' @param merge.tolerance the maximum distance below which adjacent segments
#' with equal copy number are merged. Defaults to 10^5 bp.
#' @param ignore.XY Ignore allosomes. Default TRUE.
#' @param min.vaf Remove variants with vcf below threshold. Default 0.01.
#' @param min.depth Minimum required depth for a variant to be considered.
#' Default 30.
#' @param vcf.info.af The string encoding the allele frequency field in the
#' FORMAT column of the .vcf file. Defaults to `AF`and will be ignored if
#' `vcf.source` != `sentieon`.
#' @param vcf.info.dp The string encoding the read depth field in the FORMAT
#' column of the .vcf file. Defaults to `DP`and will be ignored if `vcf.source`
#'  != `sentieon`.
#' @param min.seg.size the minimal segment length to be included in the
#' quantification.
#' @param fp.mean optional, the average false positive rate of clonal mutations
#' (e.g., due to incomplete tissue sampling). Defaults to 0.
#' @param fp.sd optional, the standard deviation of the false positive rate of
#' clonal mutations (e.g., due to incomplete tissue sampling). Defaults to 0.
#' @param excl.chr a vector of chromosomes that should be excluded from the
#' quantification. e.g., due to reporter constructs in animal models.
#' @param ref.build Reference genome. Default `hg19`. Can be `hg18`, `hg19` or
#' `hg38`.
#' @param filter.value The FILTER column value for variants that passed the
#' filtering, defaults to PASS.
#' @param sig.assign Logical. If TRUE, each variant will be assigned to the most
#'  likely mutational signature.
#' @param assign.method Method to assign signatures: "max" to assign the
#' signature with the highest probability, "sample" to randomly assign based on
#' signature probabilities.
#' @param sig.file The path to the output file from `SigProfilerAssignment`,
#' typically named "Decomposed_MutationType_Probabilities.txt". If `NULL` and
#' `sig.assign = TRUE`, signatures will be assigned using functions from
#' `MutationalPatterns`.
#' @param sig.select A character vector of specific signatures to include in the
#'  analysis (e.g., c("SBS1", "SBS5", "SBS40") to focus on clock-like mutational
#'   processes).
#' @param min.p Numeric. The minimum probability threshold from the
#' SigAssignment output that a variant must meet to be considered as matching a
#' specific signature.
#' @param driver.file optional, path to file with "chrom", "snv_start", "ref",
#' "alt", "gene" column containing known driver SNVs.
#' @param ... further arguments and parameters passed to LACHESIS functions.
#' @examples
#' # An example file with sample annotations and meta data
#' input.files <- system.file("extdata", "Sample_template.txt",
#'     package = "LACHESIS"
#' )
#' input.files <- data.table::fread(input.files)
#'
#' # cnv and snv files for example tumors
#' nbe11 <- list.files(system.file("extdata/NBE11/", package = "LACHESIS"),
#'     full.names = TRUE
#' )
#' nbe15 <- list.files(system.file("extdata/NBE15/", package = "LACHESIS"),
#'     full.names = TRUE
#' )
#' nbe26 <- list.files(system.file("extdata/NBE26/", package = "LACHESIS"),
#'     full.names = TRUE
#' )
#'
#' cnv.file <- c(nbe11[1], nbe15[1], nbe26[1])
#' snv.file <- c(nbe11[2], nbe15[2], nbe26[2])
#'
#' input.files$cnv.file <- cnv.file
#' input.files$snv.file <- snv.file
#'
#' # Make an example input file with paths to cnv and snv file along with other
#' # meta data
#' lachesis_input <- tempfile(
#'     pattern = "lachesis", tmpdir = tempdir(),
#'     fileext = ".tsv"
#' )
#' data.table::fwrite(x = input.files, file = lachesis_input, sep = "\t")
#'
#' # Example with template file with paths to multiple cnv/snv files as an input
#' lachesis <- LACHESIS(input.files = lachesis_input)
#'
#' # Example with a single sample input
#' strelka_vcf <- system.file("extdata", "strelka2.somatic.snvs.vcf.gz",
#'     package = "LACHESIS"
#' )
#' aceseq_cn <- system.file("extdata",
#'     "ACESeq/NBE11_comb_pro_extra2.59_0.83.txt",
#'     package = "LACHESIS"
#' )
#' lachesis <- LACHESIS(
#'     ids = "NBE11", cnv.files = aceseq_cn,
#'     snv.files = strelka_vcf, vcf.source = "strelka", purity = 0.83,
#'     ploidy = 2.59
#' )
#'
#' # Example with multiple sample and data frame input
#' nbe11_vcf <- system.file("extdata",
#'     "NBE11/snvs_NBE11_somatic_snvs_conf_8_to_10.vcf",
#'     package = "LACHESIS"
#' )
#' nbe11_cn <- read.delim(
#'     system.file("extdata",
#'         "NBE11/NBE11_comb_pro_extra2.59_0.83.txt",
#'         package = "LACHESIS"
#'     ),
#'     sep = "\t",
#'     header = TRUE
#' )
#' nbe15_vcf <- system.file("extdata",
#'     "NBE15/snvs_NBE15_somatic_snvs_conf_8_to_10.vcf",
#'     package = "LACHESIS"
#' )
#' nbe15_cn <- read.delim(
#'     system.file("extdata",
#'         "NBE15/NBE15_comb_pro_extra2.51_1.txt",
#'         package = "LACHESIS"
#'     ),
#'     sep = "\t",
#'     header = TRUE
#' )
#' lachesis <- LACHESIS(
#'     ids = c("NBE11", "NBE15"), cnv.files =
#'         list(nbe11_cn, nbe15_cn), snv.files = c(nbe11_vcf, nbe15_vcf),
#'     vcf.source = c("dkfz", "dkfz"), purity = c(0.83, 1), ploidy = c(2.59, 2.51),
#'     cnv.chr.col = c(1, 1), cnv.start.col = c(2, 2), cnv.end.col = c(3, 3),
#'     cnv.A.col = c(34, 34), cnv.B.col = c(35, 35), cnv.tcn.col = c(37, 37)
#' )
#'
#' @seealso \code{\link{MRCA}} \code{\link{clonalMutationCounter}}
#' \code{\link{normalizeCounts}}
#' @import tidyr
#' @import ggplot2
#' @importFrom utils packageVersion
#' @importFrom stats setNames
#' @return a data.table
#' @export

LACHESIS <- function(input.files = NULL, ids = NULL, vcf.tumor.ids = NULL,
                     cnv.files = NULL, snv.files = NULL, vcf.source = NULL,
                     purity = NULL, ploidy = NULL, cnv.chr.col = NULL,
                     cnv.start.col = NULL, cnv.end.col = NULL, cnv.A.col = NULL,
                     cnv.B.col = NULL, cnv.tcn.col = NULL, age = NULL,
                     OS.time = NULL, OS = NULL, EFS.time = NULL, EFS = NULL,
                     output.dir = NULL, ignore.XY = TRUE, min.cn = 1,
                     max.cn = 4, merge.tolerance = 10^5, min.vaf = 0.01,
                     min.depth = 30, vcf.info.af = "AF", vcf.info.dp = "DP",
                     min.seg.size = 10^7, fp.mean = 0, fp.sd = 0,
                     excl.chr = NULL, ref.build = "hg19",
                     filter.value = "PASS", sig.assign = FALSE, sig.file = NULL,
                     assign.method = "sample", sig.select = NULL, min.p = NULL,
                     driver.file = NULL, ...) {
    ID <- cnv.file <- snv.file <- fwrite <- known_driver_gene <- Sample <-
        Clonality <- NULL

    if (is.null(input.files) & is.null(cnv.files)) {
        stop("Missing input file!")
    } else if (is.null(input.files)) {
        if (any(is.null(cnv.files), is.null(snv.files))) {
            stop("Missing snv and cnv inputs!")
        } else if (!(is.data.frame(cnv.files) || is.data.table(cnv.files))) {
            if (length(cnv.files) != length(snv.files)) {
                stop("Please provide snv and cnv input for every sample!")
            }
        }
    }

    ref.build <- match.arg(
        arg = ref.build, choices = c("hg19", "hg18", "hg38"),
        several.ok = FALSE
    )

    incl.chr <- setdiff(seq_len(22), excl.chr)
    if (!ignore.XY) {
        incl.chr <- c(incl.chr, "X", "Y")
    }

    clonality_list <- list()

    # Collect ECA and MRCA densities for each tumor
    cohort.densities <- data.table::data.table(
        Sample_ID = character(),
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
        EFS = numeric()
    )

    # Initializing datatable for logfile
    if (!is.null(input.files)) {
        col.class <- fread(input.files)
        if (is.null(col.class$cnv.chr.col)) {
            col.class <- NA
        } else if (is.numeric(col.class$cnv.chr.col)) {
            col.class <- numeric()
        } else if (is.character(col.class$cnv.chr.col)) {
            col.class <- character()
        }

        log.file.data.cohort <- data.table::data.table(
            Sample_ID = character(),
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
            snv.file = character()
        )

        sample.specs <- data.table::fread(input.files,
            sep = "\t",
            stringsAsFactors = FALSE
        )

        if (any(is.na(sample.specs[, ID]))) {
            tmp1 <- sample.specs[, which(is.na(ID))]
            tmp2 <- sample.specs[, sum(is.na(ID))]
            warning(
                "No sample name provided for samples ", tmp1, ";
                sample name was set to 1 - ", tmp2
            )
            rm(tmp1, tmp2)
            sample.specs[, ID := as.character(ID)][
                is.na(ID),
                ID := which(is.na(ID))
            ]
        }

        if (any(is.na(sample.specs[, cnv.file]))) {
            tmp1 <- toString(sample.specs[, ID[which(is.na(cnv.file))]])
            warning(
                "No CNV file provided for sample(s) ", tmp1, ";
                sample(s) will be excluded"
            )
            rm(tmp1)
            sample.specs[!is.na(cnv.file), ]
            if (nrow(sample.specs) == 0) {
                stop("No files retained! Stopping analysis.")
            }
        }

        if (any(is.na(sample.specs[, snv.file]))) {
            tmp1 <- toString(sample.specs[, ID[which(is.na(snv.file))]])
            warning(sprintf(
                "No SNV file provided for sample(s) %s; will be excluded",
                tmp1
            ))
            rm(tmp1)
            sample.specs[!is.na(snv.file), ]
            if (nrow(sample.specs) == 0) {
                stop("No files retained! Stopping analysis.")
            }
        }

        sample.specs.spl <- split(sample.specs, sample.specs$ID)

        if (any(lapply(sample.specs.spl, nrow) > 1)) {
            stop("Duplicated IDs found!")
        }

        for (i in seq_len(length(sample.specs.spl))) {
            x <- sample.specs.spl[[i]]

            if (is.null(x$ID)) {
                stop("Please provide sample identifiers.")
            }
            if (is.null(x$vcf.source)) {
                stop("Please provide vcf source.")
            }
            vcf.source <- match.arg(
                arg = x$vcf.source, choices = c(
                    "strelka", "mutect",
                    "sentieon", "dkfz"
                ),
                several.ok = FALSE
            )
            if (is.null(x$vcf.tumor.ids)) {
                x$vcf.tumor.ids <- x$ID
            } else if (any(is.na(x$vcf.tumor.ids))) {
                tmp1 <- which(is.na(x$vcf.tumor.ids))
                warning(sprintf(
                    "No column ID provided for sample %s; will be inferred.",
                    tmp1
                ))
                x$vcf.tumor.ids[is.na(x$vcf.tumor.ids)] <- x$id[
                    is.na(x$vcf.tumor.ids)
                ]
            }
            message("Computing SNV density for sample ", x$ID)

            if (!is.null(output.dir)) {
                dir.create(paste(output.dir, x$ID, sep = "/"),
                    recursive = TRUE,
                    showWarnings = FALSE
                ) # Create per-sample output directory
            } else {
                warning("No output directory specified.
                LACHESIS output will be discarded.")
            }

            cnv <- readCNV(
                cn.info = x$cnv.file, chr.col = x$cnv.chr.col,
                start.col = x$cnv.start.col, end.col = x$cnv.end.col,
                A.col = x$cnv.A.col, B.col = x$cnv.B.col,
                tcn.col = x$cnv.tcn.col, tumor.id = x$ID,
                merge.tolerance = merge.tolerance,
                max.cn = max.cn, ignore.XY = ignore.XY
            )

            snv <- readVCF(
                vcf = x$snv.file, vcf.source = x$vcf.source,
                t.sample = x$vcf.tumor.id, min.depth = min.depth,
                min.vaf = min.vaf, info.af = vcf.info.af,
                ignore.XY = ignore.XY, info.dp = vcf.info.dp,
                filter.value = filter.value, ...
            )

            nb <- nbImport(
                cnv = cnv, snv = snv, purity = x$purity, ploidy = x$ploidy,
                sig.assign = sig.assign, assign.method = assign.method,
                ID = x$ID, sig.file = sig.file, sig.select = sig.select,
                min.p = min.p, ref.build = ref.build, ...
            )

            if (nrow(nb) == 0) {
                warning("Insufficient data for sample ", x$ID)
                this.tumor.density <- data.table::data.table(
                    Sample_ID = x$ID,
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
                    EFS = x$EFS
                )
                next
            }

            if (!is.null(output.dir)) {
                plotVAFdistr(snv, output.file = paste(output.dir, x$ID,
                    "01_VAF_histogram.pdf",
                    sep = "/"
                ), ...)
            }

            raw.counts <- clonalMutationCounter(
                nbObj = nb, min.cn = min.cn,
                max.cn = max.cn, chromosomes = incl.chr
            )
            norm.counts <- normalizeCounts(countObj = raw.counts)
            mrca <- MRCA(
                normObj = norm.counts, min.seg.size = min.seg.size,
                fp.mean = fp.mean, excl.chr = excl.chr
            )

            this.tumor.density <- data.table::data.table(
                Sample_ID = x$ID,
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
                EFS = x$EFS
            )

            cohort.densities <- merge(cohort.densities, this.tumor.density,
                all = TRUE
            )

            # Output the result for this sample
            if (!is.null(output.dir)) {
                mrca_colnames <- c(
                    "purity", "ploidy", "MRCA_time_mean", "MRCA_time_lower",
                    "MRCA_time_upper", "ECA_time_mean", "ECA_time_lower",
                    "ECA_time_upper"
                )
                mrca.densities <- transpose(data.table(unlist(attributes(mrca)[
                    mrca_colnames
                ])))
                setnames(mrca.densities, mrca_colnames)
                data.table::fwrite(mrca.densities,
                    file = file.path(
                        output.dir, x$ID,
                        paste0(
                            "03_MRCA_densities_",
                            x$ID, ".txt"
                        )
                    ),
                    quote = FALSE, col.names = TRUE, sep = "\t"
                )
                data.table::fwrite(mrca,
                    file = file.path(
                        output.dir, x$ID,
                        paste0(
                            "04_SNV_timing_per_segment_",
                            x$ID, ".txt"
                        )
                    ),
                    row.names = FALSE, quote = FALSE, sep = "\t"
                )
            }

            if (!is.null(output.dir)) {
                plotMutationDensities(
                    mrcaObj = mrca, samp.name = x$ID,
                    output.file = paste(output.dir, x$ID,
                        "05_Evolutionary_timeline.pdf",
                        sep = "/"
                    ), ...
                )
            }

            snvClonality <- estimateClonality(
                nbObj = nb, mrcaObj = mrca, ID = x$ID,
                purity = x$purity,
                driver.file = driver.file,
                ref.build = ref.build
            )
            clonality_list[[i]] <- snvClonality
            if (!is.null(output.dir)) {
                data.table::fwrite(snvClonality,
                    file = file.path(
                        output.dir, x$ID,
                        paste0(
                            "06_SNV_timing_per_SNV_",
                            x$ID, ".txt"
                        )
                    ),
                    quote = FALSE, col.names = TRUE, sep = "\t"
                )
                plotNB(
                    nb = nb, snvClonality = snvClonality, samp.name = x$ID,
                    output.file = paste(output.dir, x$ID,
                        "02_VAF_histogram_strat.pdf",
                        sep = "/"
                    ),
                    ref.build = ref.build, max.cn = max.cn, min.cn = min.cn, ...
                )
                plotClonality(
                    snvClonality = snvClonality, nbObj = nb,
                    sig.assign = sig.assign,
                    output.file = paste(output.dir, x$ID,
                        "07_SNV_timing_per_SNV.pdf",
                        sep = "/"
                    ),
                    ...
                )
            }

            # Collecting data for log file
            package.version <- as.character(utils::packageVersion("LACHESIS"))
            log.file.data.single <- data.table::data.table(
                Sample_ID = x$ID,
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
                snv.file = x$snv.file
            )

            log.file.data.cohort <- merge(log.file.data.cohort,
                log.file.data.single,
                all = TRUE
            )
        }
        rm(sample.specs.spl)
    } else {
        if (any(is.na(ids))) {
            tmp1 <- which(is.na(ids))
            tmp2 <- sum(is.na(ids))
            warning(sprintf(
                "No sample name provided for samples %s;
                sample name was set to 1 - %s", tmp1, tmp2
            ))
            ids[is.na(ids)] <- tmp1
            rm(tmp1, tmp2)
        }
        if (is.null(vcf.tumor.ids)) {
            warning("No column identifiers provided.")
            vcf.tumor.ids <- ids
        } else if (any(is.na(vcf.tumor.ids))) {
            tmp1 <- which(is.na(vcf.tumor.ids))
            warning(
                "No column ID provided for samples %s;
                column name will be inferred"
            )
            rm(tmp1)
            vcf.tumor.ids[is.na(vcf.tumor.ids)] <- ids[is.na(vcf.tumor.ids)]
        }

        for (i in seq_len(length(cnv.files))) {
            message("Computing SNV density for sample ", ids[i])

            if (!is.null(output.dir)) {
                dir.create(paste(output.dir, ids[i], sep = "/"),
                    recursive = TRUE,
                    showWarnings = FALSE
                ) # Create per-sample output directory
            } else {
                warning("No output directory specified.
                LACHESIS output will be discarded.")
            }

            if (is.na(cnv.files)[i]) {
                tmp1 <- ids[1]
                warning(sprintf(
                    "No CNV file provided for sample %s; sample will be excluded",
                    tmp1
                ))
                rm(tmp1)
                next
            }
            if (is.na(snv.files)[i]) {
                tmp1 <- ids[1]
                warning(sprintf(
                    "No SNV file provided for sample %s; will be excluded",
                    tmp1
                ))
                rm(tmp1)
                next
            }
            vcf.source[i] <- match.arg(
                arg = vcf.source[i], choices = c(
                    "strelka", "mutect",
                    "sentieon", "dkfz"
                ),
                several.ok = FALSE
            )

            cnv <- readCNV(
                cn.info = cnv.files[[i]], chr.col = cnv.chr.col[i],
                start.col = cnv.start.col[i], end.col = cnv.end.col[i],
                A.col = cnv.A.col[i], B.col = cnv.B.col[i],
                tcn.col = cnv.tcn.col[i], tumor.id = ids[i],
                merge.tolerance = merge.tolerance,
                max.cn = max.cn, ignore.XY = ignore.XY
            )

            snv <- readVCF(
                vcf = snv.files[i], vcf.source = vcf.source[i],
                t.sample = vcf.tumor.ids[i], min.depth = min.depth,
                min.vaf = min.vaf, info.af = vcf.info.af, ignore.XY = ignore.XY,
                info.dp = vcf.info.dp, filter.value = filter.value, ...
            )

            nb <- nbImport(
                cnv = cnv, snv = snv, purity = purity[i], ploidy = ploidy[i],
                sig.assign = sig.assign, assign.method = assign.method,
                ID = ids[i], sig.file = sig.file, sig.select = sig.select,
                min.p = min.p, ref.build = ref.build, ...
            )

            if (nrow(nb) == 0) {
                warning("Insufficient data for sample ", ids[i])
                this.tumor.density <- data.table::data.table(
                    Sample_ID = ids[i],
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
                    EFS = EFS[i]
                )

                next
            }

            if (!is.null(output.dir)) {
                plotVAFdistr(snv,
                    output.file = paste(output.dir, ids[i],
                        "01_VAF_histogram.pdf",
                        sep = "/"
                    ),
                    ...
                )
            }

            raw.counts <- clonalMutationCounter(
                nbObj = nb, min.cn = min.cn,
                max.cn = max.cn, chromosomes = incl.chr
            )
            norm.counts <- normalizeCounts(countObj = raw.counts)
            if (nrow(norm.counts) == 1) {
                tmp1 <- ids[i]
                warning(sprintf(
                    "Too few segments to estimate MRCA density for sample %s.",
                    tmp1
                ))
                rm(tmp1)
                mrca <- ""
                attr(mrca, "MRCA_time_mean") <- NA
                attr(mrca, "MRCA_time_upper") <- NA
                attr(mrca, "MRCA_time_lower") <- NA
                attr(mrca, "ECA_time_mean") <- NA
                attr(mrca, "ECA_time_lower") <- NA
                attr(mrca, "ECA_time_upper") <- NA
            } else {
                mrca <- MRCA(
                    normObj = norm.counts, min.seg.size = min.seg.size,
                    fp.mean = fp.mean, excl.chr = excl.chr
                )

                # Output the result for this sample
                if (!is.null(output.dir)) {
                    mrca_colnames <- c(
                        "purity", "ploidy", "MRCA_time_mean",
                        "MRCA_time_lower", "MRCA_time_upper", "ECA_time_mean",
                        "ECA_time_lower", "ECA_time_upper"
                    )
                    mrca.densities <- transpose(data.table(unlist(
                        attributes(mrca)[mrca_colnames]
                    )))
                    setnames(mrca.densities, mrca_colnames)
                    data.table::fwrite(mrca.densities,
                        file = file.path(
                            output.dir, ids[i],
                            paste0(
                                "03_MRCA_densities_",
                                ids[i], ".txt"
                            )
                        ),
                        quote = FALSE, col.names = TRUE, sep = "\t"
                    )
                    data.table::fwrite(mrca,
                        file = file.path(
                            output.dir, ids[i],
                            paste0(
                                "04_SNV_timing_per_segment_",
                                ids[i], ".txt"
                            )
                        ),
                        row.names = FALSE, quote = FALSE, sep = "\t"
                    )
                }

                if (!is.null(output.dir)) {
                    plotMutationDensities(
                        mrcaObj = mrca, samp.name = ids[i],
                        output.file = paste(output.dir, ids[i],
                            "05_Evolutionary_timeline.pdf",
                            sep = "/"
                        ), ...
                    )
                }
            }

            snvClonality <- estimateClonality(
                nbObj = nb, mrcaObj = mrca, ID = ids[i],
                purity = purity[i],
                driver.file = driver.file,
                ref.build = ref.build
            )
            clonality_list[[i]] <- snvClonality
            if (!is.null(output.dir)) {
                data.table::fwrite(snvClonality,
                    file = file.path(
                        output.dir, ids[i],
                        paste0(
                            "06_SNV_timing_per_SNV_",
                            ids[i], ".txt"
                        )
                    ),
                    quote = FALSE, col.names = TRUE, sep = "\t"
                )
                plotNB(
                    nb = nb, snvClonality = snvClonality, samp.name = ids[i],
                    output.file = paste(output.dir, ids[i],
                        "02_VAF_histogram_strat.pdf",
                        sep = "/"
                    ),
                    ref.build = ref.build, max.cn = max.cn, min.cn = min.cn, ...
                )
                plotClonality(
                    snvClonality = snvClonality, nbObj = nb,
                    sig.assign = sig.assign,
                    output.file = paste(output.dir, ids[i],
                        "07_SNV_timing_per_SNV.pdf",
                        sep = "/"
                    ),
                    ...
                )
            }

            this.tumor.density <- data.table::data.table(
                Sample_ID = ids[i],
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
                EFS = EFS[i]
            )

            cohort.densities <- merge(cohort.densities, this.tumor.density,
                all = TRUE
            )
        }
    }

    # Plot clonality distribution of SNVs
    clonality_cohort <- rbindlist(clonality_list, use.names = TRUE, fill = TRUE)
    if (!is.null(output.dir)) {
        output.file <- paste(output.dir, "SNV_timing_per_SNV_cohort.txt",
            sep = "/"
        )
        fwrite(clonality_cohort, output.file, sep = "\t")

        clonality_colors <- c(
            "Precnv" = "#66c2a5", "Postcnv" = "#fc8d62",
            "C" = "#8da0cb", "SC" = "#e78ac3"
        )

        driver_dt <- clonality_cohort[!is.na(known_driver_gene) &
            trimws(known_driver_gene) != ""]
        driver_dt[, Sample := factor(Sample)]

        p1 <- ggplot(driver_dt, aes(
            x = Sample, y = known_driver_gene,
            fill = Clonality
        )) +
            geom_tile(color = "white") +
            scale_fill_manual(
                values = clonality_colors,
                labels = c(
                    "Precnv" = "Clonal\n- Pre-CNV",
                    "Postcnv" = "Clonal\n- Post-CNV",
                    "C" = "Clonal\n-NOS",
                    "SC" = "Subclonal"
                )
            ) +
            labs(
                title = "Clonality of Driver Mutations",
                x = "Patient",
                y = "Gene"
            ) +
            theme_classic() +
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1),
                axis.text.y = element_text(size = 8, face = "italic")
            )

        pdf(paste(output.dir, "Driver_mutations_cohort.pdf", sep = "/"))
        print(p1)
        dev.off()
    }

    # Plot the distribution of Mutation densities at ECA and MRCA

    if (!is.null(output.dir)) {
        plotLachesis(cohort.densities,
            output.file = paste(output.dir, "SNV_densities_cohort.pdf",
                sep = "/"
            ), ...
        )
    }

    # Save log file as tsv
    if (!is.null(output.dir) && !is.null(input.files)) {
        timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
        output.file <- paste0(
            output.dir, "/LACHESIS_logfile_", timestamp,
            ".tsv"
        )
        fwrite(log.file.data.cohort, output.file, sep = "\t")
    }

    return(cohort.densities)
}


#' Plot SNV densities at ECA and MRCA
#' @description
#' Visualizes results from \code{\link{LACHESIS}}. Top plot, histograms of mean
#' mutation densities; bottom plots, cumulative distribution of mean mutation
#' densities with 95% confidence intervals.
#' @param lachesis output generated from \code{\link{LACHESIS}}
#' @param lach.suppress.outliers whether outliers (defined as the 2.5% tumors
#' with lowest and highest densities) are to be plot. Default `TRUE`.
#' @param lach.log.densities plot logarithmic densities. Default `FALSE`.
#' @param lach.col.zero optional, bar color for single-copy SSNV densities.
#' @param lach.col.multi optional, bar color for multi-copy SSNV densities.
#' @param lach.border, optional, border color for the bars.
#' @param binwidth optional, the binwidth in the histogram.
#' @param output.file optional, the file to which the plot will be stored.
#' @param ... further arguments and parameters passed to other LACHESIS
#' functions.
#' @return graph with cohort overview of SNV densities at ECA/ MRCA
#' @examples
#' # An example file with sample annotations and meta data
#' input.files <- system.file("extdata", "Sample_template.txt",
#'     package = "LACHESIS"
#' )
#' input.files <- data.table::fread(input.files)
#'
#' # cnv and snv files for example tumors
#' nbe11 <- list.files(system.file("extdata/NBE11/", package = "LACHESIS"),
#'     full.names = TRUE
#' )
#' nbe15 <- list.files(system.file("extdata/NBE15/", package = "LACHESIS"),
#'     full.names = TRUE
#' )
#' nbe26 <- list.files(system.file("extdata/NBE26/", package = "LACHESIS"),
#'     full.names = TRUE
#' )
#'
#' cnv.file <- c(nbe11[1], nbe15[1], nbe26[1])
#' snv.file <- c(nbe11[2], nbe15[2], nbe26[2])
#'
#' input.files$cnv.file <- cnv.file
#' input.files$snv.file <- snv.file
#'
#' # Make an example input file with paths to cnv and snv file along with other
#' # meta data
#' lachesis_input <- tempfile(
#'     pattern = "lachesis", tmpdir = tempdir(),
#'     fileext = ".tsv"
#' )
#' data.table::fwrite(x = input.files, file = lachesis_input, sep = "\t")
#'
#' # Example with template file with paths to multiple cnv/snv files as an input
#' lachesis <- LACHESIS(input.files = lachesis_input)
#' plotLachesis(lachesis)
#' @export
#' @importFrom graphics abline Axis box grid hist mtext par rect text title
#' @importFrom grDevices adjustcolor
#' arrows legend points polygon

plotLachesis <- function(lachesis = NULL, lach.suppress.outliers = FALSE,
                         lach.log.densities = FALSE, lach.col.multi = "#176A02",
                         lach.border = NULL, binwidth = NULL,
                         lach.col.zero = "#4FB12B", output.file = NULL, ...) {
    MRCA_time_mean <- ECA_time_mean <- NULL

    if (is.null(lachesis)) {
        stop("Missing input. Please provide the output generated by LACHESIS()")
    }
    if (nrow(lachesis) == 1) {
        warning(
            "Cannot produce summary statistics for a single case.
            Returning null."
        )
        return()
    }
    if (any(is.na(lachesis$MRCA_time_mean))) {
        tmp1 <- sum(is.na(lachesis$MRCA_time_mean))
        warning(sprintf(
            "Removing %s samples with missing MRCA density estimate.", tmp1
        ))
        lachesis <- lachesis[!is.na(MRCA_time_mean), ]
        rm(tmp1)
    }
    if (nrow(lachesis) == 0) {
        warning(
            "No sample with MRCA density estimate provided. Returning zero."
        )
        return(NULL)
    }
    if (!is.null(output.file)) {
        pdf(output.file, width = 8, height = 6)
    }

    # I. Plot histograms
    to.plot <- lachesis

    if (lach.suppress.outliers) {
        to.plot <- to.plot[MRCA_time_mean < quantile(MRCA_time_mean, 0.975) &
            MRCA_time_mean > quantile(MRCA_time_mean, 0.025), ]
    }

    if (is.null(binwidth)) {
        if (lach.log.densities) {
            binwidth <- (max(log10(to.plot$MRCA_time_mean)) -
                min(log10(to.plot$MRCA_time_mean))) / 20
        } else {
            binwidth <- (max(to.plot$MRCA_time_mean) -
                min(to.plot$MRCA_time_mean)) / 20
        }
    }

    lo_mat <- matrix(data = c(1, 3, 2, 3), nrow = 2, ncol = 2, byrow = TRUE)
    graphics::layout(
        mat = lo_mat, widths = c(1, 2),
        heights = c(1, 1)
    )

    par(mar = c(3, 4, 3, 1))

    if (lach.log.densities) {
        min.x <- floor(min(to.plot[, log10(MRCA_time_mean)]))
        max.x <- ceiling(max(to.plot[, log10(MRCA_time_mean)]))

        hist(to.plot[, log10(MRCA_time_mean)],
            xlim = c(min.x, max.x),
            breaks = 20,
            col = lach.col.zero, border = lach.border, main = NA,
            xlab = NA, ylab = NA, axes = FALSE
        )

        Axis(
            side = 1, at = seq(min.x, max.x, length.out = 10),
            labels = round(10^seq(min.x, max.x, length.out = 10), digits = 2)
        )
        Axis(side = 2)
    } else {
        hist(to.plot[, MRCA_time_mean],
            xlim = c(0, 1.05 * max(to.plot[, MRCA_time_mean])),
            breaks = seq(0, max(to.plot[, MRCA_time_mean]) * 1.05, binwidth),
            col = lach.col.zero, border = lach.border, main = NA,
            xlab = NA, ylab = NA
        )
    }


    title(main = paste("SNV densities at MRCA"), cex.main = 1)
    mtext(text = "SNVs per Mb", side = 1, line = 2, cex = 0.7)
    mtext(text = "No. of tumors", side = 2, line = 1.8, cex = 0.7)

    # Histogram of SNV density at ECA

    if (all(is.na(lachesis$ECA_time_mean)) & !is.null(output.file)) {
        dev.off()
        return()
    } else if (all(is.na(lachesis$ECA_time_mean))) {
        return()
    }

    to.plot <- lachesis

    if (lach.suppress.outliers) {
        to.plot <- to.plot[MRCA_time_mean < quantile(MRCA_time_mean, 0.975) &
            MRCA_time_mean > quantile(MRCA_time_mean, 0.025), ]
    }

    par(mar = c(3, 4, 3, 1))

    if (lach.log.densities) {
        min.x <- floor(min(to.plot[, log10(ECA_time_mean)], na.rm = TRUE))
        max.x <- ceiling(max(to.plot[, log10(ECA_time_mean)], na.rm = TRUE))

        hist(to.plot[, log10(ECA_time_mean)],
            xlim = c(min.x, max.x),
            breaks = 20,
            col = lach.col.zero, border = lach.border, main = NA,
            xlab = NA, ylab = NA, axes = FALSE
        )

        Axis(
            side = 1, at = seq(min.x, max.x, length.out = 10),
            labels = round(10^seq(min.x, max.x, length.out = 10), digits = 2)
        )
        Axis(side = 2)
    } else {
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

        hist(to.plot[, ECA_time_mean],
            xlim = c(0, max(0.01, 1.05 * max(to.plot[, ECA_time_mean],
                na.rm = TRUE
            ))),
            breaks = seq(
                0, max(0.01, max(to.plot[, ECA_time_mean], na.rm = TRUE)) *
                  1.05,
                binwidth
            ),
            col = lach.col.multi, border = lach.border, main = NA,
            xlab = NA, ylab = NA
        )
    }


    title(main = paste("SNV densities at ECA"), cex.main = 1)
    mtext(text = "SNVs per Mb", side = 1, line = 2, cex = 0.7)
    mtext(text = "No. of tumors", side = 2, line = 1.8, cex = 0.7)

    # Cumulative mutation densities at ECA and MRCA
    par(mar = c(3, 4, 3, 1), xpd = FALSE)

    x.min <- 0
    x.max <- max(lachesis$MRCA_time_upper,
        lachesis$ECA_time_upper,
        na.rm = TRUE
    ) * 1.3
    y.min <- 0
    y.max <- 1
    plot(NA, NA,
        xlim = c(x.min, x.max), ylim = c(y.min, y.max), xlab = NA, ylab = NA,
        main = NA, axes = FALSE, frame.plot = FALSE
    )
    Axis(side = 1, cex = 0.7)
    Axis(side = 2, cex = 0.7)
    mtext(text = "SNVs per Mb", side = 1, line = 2, cex = 0.7)
    mtext(text = "Fraction of tumors", side = 2, line = 2, cex = 0.7)

    to.plot.MRCA <- data.frame(
        x.lower = rep(sort(c(lachesis$MRCA_time_mean)),
            each = 2
        )[-1],
        x.upper = rep(sort(c(lachesis$MRCA_time_mean)),
            each = 2
        )[-2 * (nrow(lachesis))]
    )
    to.plot.MRCA$y.lower <- vapply(
        rep(sort(c(lachesis$MRCA_time_mean)), each = 2),
        function(x) {
            sum(lachesis$MRCA_time_upper <= x)
        },
        numeric(1)
    )[-2 * (nrow(lachesis))]
    to.plot.MRCA$y.upper <- vapply(
        rep(sort(c(lachesis$MRCA_time_mean)), each = 2),
        function(x) {
            sum(lachesis$MRCA_time_lower <= x)
        },
        numeric(1)
    )[-1]

    polygon(c(to.plot.MRCA$x.lower, rev(to.plot.MRCA$x.upper)),
        c(to.plot.MRCA$y.lower, rev(to.plot.MRCA$y.upper)) / nrow(lachesis),
        col = adjustcolor(lach.col.zero, alpha.f = 0.3),
        border = NA
    )

    plot.ecdf(lachesis$MRCA_time_mean,
        col = lach.col.zero, add = TRUE,
        verticals = TRUE
    )

    to.plot.ECA <- data.frame(
        x.lower = rep(sort(c(lachesis$ECA_time_mean)),
            each = 2
        )[-1],
        x.upper = rep(sort(c(lachesis$ECA_time_mean)),
            each = 2
        )[-2 * (nrow(lachesis[!is.na(ECA_time_mean), ]))]
    )
    to.plot.ECA$y.lower <- vapply(
        rep(sort(c(lachesis$ECA_time_mean)), each = 2),
        function(x) {
            sum(lachesis$ECA_time_upper <= x,
                na.rm = TRUE
            )
        }, numeric(1)
    )[
        -2 * (nrow(lachesis[
            !is.na(ECA_time_mean),
        ]))
    ]
    to.plot.ECA$y.upper <- vapply(
        rep(sort(c(lachesis$ECA_time_mean)), each = 2),
        function(x) {
            sum(lachesis$ECA_time_lower <= x,
                na.rm = TRUE
            )
        }, numeric(1)
    )[-1]

    polygon(c(to.plot.ECA$x.lower, rev(to.plot.ECA$x.upper)),
        c(to.plot.ECA$y.lower, rev(to.plot.ECA$y.upper)) /
            nrow(lachesis[!is.na(ECA_time_mean), ]),
        col = adjustcolor(lach.col.multi, alpha.f = 0.3),
        border = NA
    )

    legend("topright",
        legend = c("ECA", "MRCA"),
        fill = c(lach.col.multi, lach.col.zero),
        border = NA,
        bty = "o",
        cex = 0.7,
        inset = c(0.05, 0.1)
    )

    plot.ecdf(lachesis$ECA_time_mean,
        col = lach.col.multi, add = TRUE,
        verticals = TRUE
    )

    title(main = paste("Cumulative SNV densities at ECA and MRCA"),
          cex.main = 1)

    if (!is.null(output.file)) {
        dev.off()
    }
}



#' Classify a tumor's start of clonal outgrowth during tumorigenesis as "early"
#'  or "late" (favorable/ unfavorable prognosis) depending on the mutation
#'  density at its MRCA
#' @description
#' Takes SNV density timing as computed by `LACHESIS` as input and classifies
#' the tumors in the cohort.
#' @param lachesis output generated from \code{\link{LACHESIS}}.
#' @param mrca.cutpoint optional; value based on SNV_densities_cohort.pdf
#' observation, will be used as inferred from a test data set if not specified
#' by user.
#' @param infer.cutpoint logical; should the MRCA cutpoint be inferred from the
#' data?
#' @param entity optional; the tumor entity if classifying according to a
#' pre-defined threshold. Currently, only "neuroblastoma" is supported.
#' @param surv.time column name containing survival time; defaults to `OS.time`.
#' @param surv.event column name containing event; defaults to `OS`.
#' @param lach.col.zero optional, bar color for single-copy SSNV densities.
#' @param lach.col.multi optional, bar color for multi-copy SSNV densities.
#' @param surv.time.scale numeric value by which survival time is to be divided
#' (e.g., 365 for converting days into years, 30 for months), defaults to `1`.
#' @param output.dir link to directory in which output is to be stored.
#' @return data.table with binary assignment early/ late
#' @examples
#' # An example file with sample annotations and meta data
#' input.files <- system.file("extdata", "Sample_template.txt",
#'     package = "LACHESIS"
#' )
#' input.files <- data.table::fread(input.files)
#'
#' # cnv and snv files for example tumors
#' nbe11 <- list.files(system.file("extdata/NBE11/", package = "LACHESIS"),
#'     full.names = TRUE
#' )
#' nbe15 <- list.files(system.file("extdata/NBE15/", package = "LACHESIS"),
#'     full.names = TRUE
#' )
#' nbe26 <- list.files(system.file("extdata/NBE26/", package = "LACHESIS"),
#'     full.names = TRUE
#' )
#'
#' cnv.file <- c(nbe11[1], nbe15[1], nbe26[1])
#' snv.file <- c(nbe11[2], nbe15[2], nbe26[2])
#'
#' input.files$cnv.file <- cnv.file
#' input.files$snv.file <- snv.file
#'
#' # Make an example input file with paths to cnv and snv file along with other
#' # meta data
#' lachesis_input <- tempfile(
#'     pattern = "lachesis", tmpdir = tempdir(),
#'     fileext = ".tsv"
#' )
#' data.table::fwrite(x = input.files, file = lachesis_input, sep = "\t")
#'
#' # Example with template file with paths to multiple cnv/snv files as an input
#' lachesis <- LACHESIS(input.files = lachesis_input)
#' classifyLACHESIS(lachesis)
#'
#' @export
#' @import survminer

classifyLACHESIS <- function(lachesis, mrca.cutpoint = NULL,
                             infer.cutpoint = FALSE, entity = "neuroblastoma",
                             lach.col.multi = "#176A02", lach.col.zero = "#4FB12B",
                             surv.time = "OS.time", surv.event = "OS", surv.time.scale = 1,
                             output.dir = NULL) {
    MRCA_time_mean <- NULL

    if (is.null(lachesis)) {
        stop("'lachesis' dataset must be provided.")
    }

    entities <- c("neuroblastoma")
    entity <- match.arg(arg = entity, choices = entities, several.ok = FALSE)

    if (infer.cutpoint == TRUE & sum(!(is.na(lachesis[[surv.time]]))) < 2) {
        stop("Please provide survival time if inferring cutpoint de novo.")
    }

    if (infer.cutpoint == TRUE & (sum(!(is.na(lachesis[[surv.event]]))) < 2 |
        sum(lachesis[[surv.event]] != 0, na.rm = TRUE) < 2)) {
        stop(
            "Please provide survival information if inferring cutpoint de novo."
        )
    }
    message("Classifying ", entity, " samples.")

    if (infer.cutpoint == TRUE) {
        message("MRCA cutpoint will be newly inferred.")
    } else if (is.null(mrca.cutpoint)) {
        message("Samples will be classified according to established
        MRCA cutpoint for ", entity, ".")
    } else if (infer.cutpoint == FALSE & is.null(mrca.cutpoint)) {
        message("Please provide cutpoint or set `infer.cutpoint`=`TRUE`")
    } else if (infer.cutpoint == FALSE) {
        message("MRCA cutpoint taken as ", mrca.cutpoint, ".")
    }

    # Calculating MRCA cutpoint
    if (infer.cutpoint) {
        mrca.cutpoint <- survminer::surv_cutpoint(
            lachesis[!is.na(get(surv.time)), .SD],
            time = surv.time,
            event = surv.event,
            variables = c("MRCA_time_mean")
        )
        mrca.cutpoint <- as.numeric(mrca.cutpoint$cutpoint[
            "MRCA_time_mean",
            "cutpoint"
        ])
    } else if (is.null(mrca.cutpoint)) {
        mrca.cutpoint <- .getCutpoint(entity)
    }

    # Categorizing according to MRCA
    lachesis.categorized <- lachesis
    lachesis.categorized$MRCA_timing <-
        ifelse(lachesis.categorized$MRCA_time_mean <
            mrca.cutpoint, "Early MRCA", "Late MRCA")
    lachesis.categorized$MRCA_timing <- factor(lachesis.categorized$MRCA_timing,
        levels = c("Early MRCA", "Late MRCA")
    )

    attr(lachesis.categorized, "MRCA Cutpoint") <- mrca.cutpoint
    attr(lachesis.categorized, "Entity") <- entity

    if (!is.null(output.dir)) {
        pdf(file = file.path(output.dir, "MRCA_Classification.pdf"), width = 10, height = 5)
    }

    lachesis.categorized[, MRCA_timing := fifelse(
        MRCA_time_mean < mrca.cutpoint,
        "Early MRCA",
        "Late MRCA"
    )]

    par(mfrow = c(1, 2), mar = c(3, 4, 4, 1), xpd = FALSE)

    x.min <- 0
    x.max <- max(
        lachesis.categorized$MRCA_time_mean,
        lachesis.categorized$ECA_time_mean,
        na.rm = TRUE
    ) * 1.3

    plot_early_late_ecdf <- function(dt, title) {
        plot(
            NA, NA,
            xlim = c(x.min, x.max),
            ylim = c(0, 1),
            xlab = NA,
            ylab = NA,
            axes = FALSE,
            frame.plot = FALSE
        )

        Axis(side = 1, cex = 0.7)
        Axis(side = 2, cex = 0.7)
        mtext("SNVs per Mb", side = 1, line = 2, cex = 0.7)
        mtext("Fraction of tumors", side = 2, line = 2, cex = 0.7)
        mtext(title, side = 3, line = 2, cex = 0.9, font = 2)
        mtext(paste("cutpoint =", mrca.cutpoint, "SNVs per Mb"), side = 3, line = 1, cex = 0.7, font = 1)

        ## MRCA
        to.plot.MRCA <- data.frame(
            x.lower = rep(sort(c(dt$MRCA_time_mean)),
                each = 2
            )[-1],
            x.upper = rep(sort(c(dt$MRCA_time_mean)),
                each = 2
            )[-2 * (nrow(dt))]
        )
        to.plot.MRCA$y.lower <- vapply(
            rep(sort(c(dt$MRCA_time_mean)), each = 2),
            function(x) {
                sum(dt$MRCA_time_upper <= x)
            },
            numeric(1)
        )[-2 * (nrow(dt))]
        to.plot.MRCA$y.upper <- vapply(
            rep(sort(c(dt$MRCA_time_mean)), each = 2),
            function(x) {
                sum(dt$MRCA_time_lower <= x)
            },
            numeric(1)
        )[-1]

        polygon(c(to.plot.MRCA$x.lower, rev(to.plot.MRCA$x.upper)),
            c(to.plot.MRCA$y.lower, rev(to.plot.MRCA$y.upper)) / nrow(dt),
            col = adjustcolor(lach.col.zero, alpha.f = 0.3), border = NA
        )
        plot.ecdf(
            dt$MRCA_time_mean,
            add = TRUE,
            verticals = TRUE,
            col = "#4FB12B",
            lwd = 2,
            pch = 19,
            cex = 0.8
        )

        ## ECA
        if ("ECA_time_mean" %in% colnames(dt)) {
            eca_dt <- dt[!is.na(ECA_time_mean)]
            if (nrow(eca_dt) > 1) {
                to.plot.ECA <- data.frame(
                    x.lower = rep(sort(c(eca_dt$ECA_time_mean)),
                        each = 2
                    )[-1],
                    x.upper = rep(sort(c(eca_dt$ECA_time_mean)),
                        each = 2
                    )[-2 * (nrow(eca_dt[!is.na(ECA_time_mean), ]))]
                )
                to.plot.ECA$y.lower <- vapply(
                    rep(sort(c(eca_dt$ECA_time_mean)), each = 2),
                    function(x) {
                        sum(eca_dt$ECA_time_upper <= x,
                            na.rm = TRUE
                        )
                    }, numeric(1)
                )[
                    -2 * (nrow(eca_dt[
                        !is.na(ECA_time_mean),
                    ]))
                ]
                to.plot.ECA$y.upper <- vapply(
                    rep(sort(c(eca_dt$ECA_time_mean)), each = 2),
                    function(x) {
                        sum(eca_dt$ECA_time_lower <= x,
                            na.rm = TRUE
                        )
                    }, numeric(1)
                )[-1]

                polygon(c(to.plot.ECA$x.lower, rev(to.plot.ECA$x.upper)),
                    c(to.plot.ECA$y.lower, rev(to.plot.ECA$y.upper)) /
                        nrow(eca_dt[!is.na(ECA_time_mean), ]),
                    col = adjustcolor(lach.col.multi, alpha.f = 0.3), border = NA
                )
                plot.ecdf(
                    eca_dt$ECA_time_mean,
                    add = TRUE,
                    verticals = TRUE,
                    col = "#176A02",
                    lwd = 2,
                    pch = 19,
                    cex = 0.8
                )
            }
        }
    }
    plot_early_late_ecdf(
        lachesis.categorized[MRCA_timing == "Early MRCA"],
        "Early MRCA"
    )

    plot_early_late_ecdf(
        lachesis.categorized[MRCA_timing == "Late MRCA"],
        "Late MRCA"
    )

    if (!is.null(output.dir)) {
        dev.off()
    }

    if (!is.null(output.dir)) {
        data.table::fwrite(lachesis.categorized,
            file = file.path(output.dir, "Lachesis_classifier.txt"),
            sep = "\t"
        )
    }

    return(lachesis.categorized)
}

# Cutpoint for neuroblastoma
.getCutpoint <- function(entity = "neuroblastoma") {
    if (entity == "neuroblastoma") {
        cut.point <- 0.05
    } else {
        stop("Available entities: neuroblastoma")
    }

    cut.point
}
