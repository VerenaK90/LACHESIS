# Load required packages
library(LACHESIS)
library(data.table)

# Test readVCF function

# Test 1: Read Strelka VCF file
test_readVCF_strelka <- function() {
    strelka_vcf <- system.file("extdata", "strelka2.somatic.snvs.vcf.gz",
                               package = "LACHESIS")
    result <- readVCF(vcf = strelka_vcf, vcf.source = "strelka", filter.value = "PASS")

    expect_true(is.data.table(result),
                info = "readVCF should return a data.table")
    expect_true(nrow(result) > 0,
                info = "Strelka VCF should contain variants")
    expect_true(all(c("chrom", "pos", "ref", "alt", "t_vaf", "t_depth") %in% names(result)),
                info = "Output should contain required columns")
}

# Test 2: Read Mutect VCF file
test_readVCF_mutect <- function() {
    mutect_vcf <- system.file("extdata", "mutect.somatic.vcf.gz",
                              package = "LACHESIS")
    result <- readVCF(vcf = mutect_vcf, vcf.source = "mutect", filter.value = ".")

    expect_true(is.data.table(result),
                info = "readVCF should return a data.table")
    expect_true(nrow(result) > 0,
                info = "Mutect VCF should contain variants")
}

# Test 3: Read DKFZ VCF file
test_readVCF_dkfz <- function() {
    dkfz_vcf <- system.file("extdata", "NBE15",
                            "snvs_NBE15_somatic_snvs_conf_8_to_10.vcf",
                            package = "LACHESIS")
    result <- readVCF(vcf = dkfz_vcf, vcf.source = "dkfz")

    expect_true(is.data.table(result),
                info = "readVCF should return a data.table")
    expect_true(nrow(result) > 0,
                info = "DKFZ VCF should contain variants")
    expect_true(all(c("chrom", "pos", "ref", "alt", "t_vaf") %in% names(result)),
                info = "Output should contain variant information")
}

# Test 4: Invalid VCF source
test_readVCF_invalid_source <- function() {
    dkfz_vcf <- system.file("extdata", "NBE15",
                            "snvs_NBE15_somatic_snvs_conf_8_to_10.vcf",
                            package = "LACHESIS")

    expect_error(readVCF(vcf = dkfz_vcf, vcf.source = "invalid_caller"),
                 info = "Invalid VCF source should raise error")
}

# Test 5: Output column types
test_readVCF_column_types <- function() {
    strelka_vcf <- system.file("extdata", "strelka2.somatic.snvs.vcf.gz",
                               package = "LACHESIS")
    result <- readVCF(vcf = strelka_vcf, vcf.source = "strelka")

    expect_true(is.character(result$chrom) || is.numeric(result$chrom),
                info = "chrom should be character or numeric")
    expect_true(is.integer(result$pos) || is.numeric(result$pos),
                info = "pos should be numeric")
    expect_true(is.character(result$ref),
                info = "ref should be character")
    expect_true(is.character(result$alt),
                info = "alt should be character")
    expect_true(is.numeric(result$t_vaf),
                info = "t_vaf should be numeric")
}

# Test 6: VAF and depth filtering
test_readVCF_vaf_filtering <- function() {
    strelka_vcf <- system.file("extdata", "strelka2.somatic.snvs.vcf.gz",
                               package = "LACHESIS")

    # Strict filtering
    result_strict <- readVCF(vcf = strelka_vcf, vcf.source = "strelka",
                             min.vaf = 0.2, min.depth = 50)

    # Lenient filtering
    result_lenient <- readVCF(vcf = strelka_vcf, vcf.source = "strelka",
                              min.vaf = 0.01, min.depth = 10)

    # Lenient should have more or equal variants
    expect_true(nrow(result_lenient) >= nrow(result_strict),
                info = "Lenient filtering should have more/equal variants than strict")

    # Check that all VAFs meet threshold
    expect_true(all(result_strict$t_vaf >= 0.2),
                info = "All variants should meet min.vaf threshold")
}

# Test 7: Ignore XY chromosomes
test_readVCF_ignore_XY <- function() {
    strelka_vcf <- system.file("extdata", "strelka2.somatic.snvs.vcf.gz",
                               package = "LACHESIS")

    result_with_XY <- readVCF(vcf = strelka_vcf, vcf.source = "strelka",
                              ignore.XY = FALSE)
    result_no_XY <- readVCF(vcf = strelka_vcf, vcf.source = "strelka",
                            ignore.XY = TRUE)

    # Should have fewer or equal variants when ignoring XY
    expect_true(nrow(result_no_XY) <= nrow(result_with_XY),
                info = "Ignoring XY should not increase variant count")
}

# Test 8: Non-existent file
test_readVCF_nonexistent_file <- function() {
    expect_error(readVCF(vcf = "/nonexistent/path/file.vcf", vcf.source = "strelka"),
                 info = "Non-existent file should raise error")
}

# Test 9: Multiple samples in VCF
test_readVCF_multiple_samples <- function() {
    # Strelka hardcodes TUMOR sample
    strelka_vcf <- system.file("extdata", "strelka2.somatic.snvs.vcf.gz",
                               package = "LACHESIS")
    result <- readVCF(vcf = strelka_vcf, vcf.source = "strelka")

    expect_true(is.data.table(result),
                info = "Should handle VCF with multiple samples")
}

# Test 10: Empty result handling
test_readVCF_no_passing_variants <- function() {
    # Use very stringent filtering to get few/no variants
    strelka_vcf <- system.file("extdata", "strelka2.somatic.snvs.vcf.gz",
                               package = "LACHESIS")

    # This should either return empty DT or handle gracefully
    result <- readVCF(vcf = strelka_vcf, vcf.source = "strelka",
                      min.vaf = 0.9)  # Very high VAF threshold

    # Should return data.table even if empty
    expect_true(is.data.table(result),
                info = "Should return data.table even with no passing variants")
}

# Run all tests in this file
for (test_func in ls(pattern = "^test_")) {
    get(test_func)()
}
