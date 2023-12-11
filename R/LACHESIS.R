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
#' snvs <- system.file("extdata", "snvs_NBE15_somatic_snvs_conf_8_to_10.vcf", package = "LACHESIS")
#' s_data <- readVCF(vcf = snvs, vcf.source = "dkfz")
#' aceseq_cn <- system.file("extdata", "NBE15_comb_pro_extra2.51_1.txt", package = "LACHESIS")
#' c_data <- readCNV(aceseq_cn)
#' nb <- nbImport(cnv = c_data, snv = s_data, purity = 1, ploidy = 2.51)
#' @seealso \code{\link{plotNB}}
#' @return a data.table
#' @export

LACHESIS <- function(input.files = NULL, ids = NULL, cnv.files = NULL, snv.files = NULL, vcf.source = "strelka",
                     purity = NULL, ploidy = NULL,
                     cnv.chr.col = NULL, cnv.start.col = NULL, cnv.end.col = NULL, cnv.A.col = NULL,
                     cnv.B.col = NULL, cnv.tcn.col = NULL, ...){

  if(is.null(input.files) & is.null(cnf.files)){
    stop("Missing input file!")
  }else if(is.null(input.files)){
    if(any(is.null(cnv.files), is.null(snv.files))){
      stop("Missing snv and cnv inputs!")
    }else if(length(cnv.files) != length(snv.files)){
      stop("Please provide snv and cnv input for every sample!")
    }
  }

  if(!is.null(input.files)){

    sample.specs <- data.table::fread(input.files, sep = "\t", header = T, stringsAsFactors = F)

    sample.specs.spl <- split(sample.specs, sample.specs$ID)

    for(i in 1:length(sample.specs.spl)){

      x <- sample.specs.spl[[i]]
      x[,which(sapply(x, is.na)):=NULL] # remove NA entries

      message("Computing SNV density for sample ", x$ID)

      cnv <- readCNV(cn.info = x$cnv.file, chr.col = x$cnv.chr.col, start.col = x$cnv.start.col,
                     end.col = x$cnv.end.col, A.col = x$cnv.A.col, B.col = x$cnv.B.col,
                     tcn.col = x$cnv.tcn.col, tumor.id = x$ID, ...)

      snv <- readVCF(vcf = x$snv.file, ignore.XY = ignore.XY, vcf.source = x$vcf.source, t.sample = x$ID, ...)

      nb <- nbImport(cnv = cnv, snv = snv, purity = x$purity, ploidy = x$ploidy)
      nb.p1 <- plotNB(nb = nb, samp.name = x$ID, ...)

      raw.counts <- clonalMutationCounter(nbObj = nb, ...)
      norm.counts <- normalizeCounts(countObj = raw.counts)
      mrca <- MRCA(normObj = norm.counts, ...)

      mrca.p2 <- plotMutationDensities(mrcaObj = mrca, samp.name = x$ID, ...)

    }
    rm(sample.specs.spl)
  }

  # which information to return? generate per-sample output file?
}
