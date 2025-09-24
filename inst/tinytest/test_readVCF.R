# Test readVCF function

# Test parameter validation
expect_error(readVCF(), "Missing input VCF file!")
expect_error(readVCF(vcf = "nonexistent.vcf"), pattern = "does not exist")
expect_error(readVCF(vcf = tempfile(), vcf.source = "invalid"), pattern = "should be one of")

# Test with example data files
mutect_vcf <- system.file("extdata", "mutect.somatic.vcf.gz", package = "LACHESIS")
strelka_vcf <- system.file("extdata", "strelka2.somatic.snvs.vcf.gz", package = "LACHESIS")
dkfz_vcf <- system.file("extdata", "NBE15", "snvs_NBE15_somatic_snvs_conf_8_to_10.vcf", package = "LACHESIS")

# Test mutect VCF if file exists
if (file.exists(mutect_vcf)) {
  mutect_data <- readVCF(vcf = mutect_vcf, vcf.source = "mutect", filter.value = ".")
  expect_true(is.data.table(mutect_data))
  expect_true(all(c("chrom", "pos", "ref", "alt", "t_vaf", "t_depth") %in% colnames(mutect_data)))
  expect_true(all(mutect_data$t_vaf >= 0 & mutect_data$t_vaf <= 1))
  expect_true(all(mutect_data$t_depth >= 0))
}

# Test strelka VCF if file exists
if (file.exists(strelka_vcf)) {
  strelka_data <- readVCF(vcf = strelka_vcf, vcf.source = "strelka")
  expect_true(is.data.table(strelka_data))
  expect_true(all(c("chrom", "pos", "ref", "alt", "t_vaf", "t_depth") %in% colnames(strelka_data)))
  expect_true(all(strelka_data$t_vaf >= 0 & strelka_data$t_vaf <= 1))
}

# Test DKFZ VCF if file exists
if (file.exists(dkfz_vcf)) {
  dkfz_data <- readVCF(vcf = dkfz_vcf, vcf.source = "dkfz")
  expect_true(is.data.table(dkfz_data))
  expect_true(all(c("chrom", "pos", "ref", "alt", "t_vaf", "t_depth") %in% colnames(dkfz_data)))
}

# Test filtering parameters
if (file.exists(strelka_vcf)) {
  # Test VAF filtering
  data_low_vaf <- readVCF(vcf = strelka_vcf, vcf.source = "strelka", min.vaf = 0.3)
  expect_true(all(data_low_vaf$t_vaf >= 0.3))

  # Test depth filtering
  data_high_depth <- readVCF(vcf = strelka_vcf, vcf.source = "strelka", min.depth = 100)
  expect_true(all(data_high_depth$t_depth >= 100))
}