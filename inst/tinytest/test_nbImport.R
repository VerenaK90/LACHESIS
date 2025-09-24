# Test nbImport function

# Test parameter validation
expect_error(nbImport(), "Missing snv and cnv inputs!")
expect_error(nbImport(cnv = data.frame(), snv = NULL), "Missing snv and cnv inputs!")
expect_error(nbImport(cnv = data.frame(), snv = data.frame()), "Missing purity and ploidy inputs!")

# Test BSgenome availability checks for hg18 and hg38
cnv_mock <- data.table::data.table(chr = 1, start = 1000, end = 2000, A = 1, B = 1, TCN = 2)
snv_mock <- data.table::data.table(chr = 1, pos = 1500, ref = "A", alt = "T",
                                  t_ref_count = 10, t_alt_count = 5, t_depth = 15, t_vaf = 0.33)

# Test hg19 (should work - it's in Imports)
if (requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)) {
  result_hg19 <- nbImport(cnv = cnv_mock, snv = snv_mock,
                         purity = 1, ploidy = 2, ref.build = "hg19")
  expect_true(data.table::is.data.table(result_hg19))
}

# Test hg18 error message when package not available
if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg18", quietly = TRUE)) {
  expect_error(
    nbImport(cnv = cnv_mock, snv = snv_mock, purity = 1, ploidy = 2,
            sig.assign = TRUE, ref.build = "hg18", ID = "test",
            sig.file = tempfile()),
    "BSgenome.Hsapiens.UCSC.hg18 is required for hg18"
  )
}

# Test hg38 error message when package not available
if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
  expect_error(
    nbImport(cnv = cnv_mock, snv = snv_mock, purity = 1, ploidy = 2,
            sig.assign = TRUE, ref.build = "hg38", ID = "test",
            sig.file = tempfile()),
    "BSgenome.Hsapiens.UCSC.hg38 is required for hg38"
  )
}

# Test with example data if available
snvs_file <- system.file("extdata", "NBE15", "snvs_NBE15_somatic_snvs_conf_8_to_10.vcf", package = "LACHESIS")
cnv_file <- system.file("extdata", "NBE15", "NBE15_comb_pro_extra2.51_1.txt", package = "LACHESIS")

if (file.exists(snvs_file) && file.exists(cnv_file)) {
  snv_data <- readVCF(vcf = snvs_file, vcf.source = "dkfz")
  cnv_data <- readCNV(cnv_file)

  nb_result <- nbImport(cnv = cnv_data, snv = snv_data, purity = 1, ploidy = 2.51)

  expect_true(is.data.table(nb_result))
  expect_true(all(c("snv_start", "snv_end", "cn_start", "cn_end") %in% colnames(nb_result)))
  expect_true(attr(nb_result, "purity") == 1)
  expect_true(attr(nb_result, "ploidy") == 2.51)
  expect_true(nrow(nb_result) > 0)
}

# Test utility functions
expect_true(is.numeric(.getContigLens("hg19")))
expect_true(is.numeric(.getContigLens("hg18")))
expect_true(is.numeric(.getContigLens("hg38")))
expect_error(.getContigLens("invalid"), "Available reference builds")

# Test expected clonal VAF calculation
vafs <- .expectedClVAF(CN = 2, purity = 1)
expect_true(is.numeric(vafs))
expect_true(all(vafs > 0 & vafs <= 1))
expect_equal(length(vafs), 2)
