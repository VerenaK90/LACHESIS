# Load required packages
library(LACHESIS)
library(data.table)

# Test nbImport function - Core integration and signature assignment

# Helper to load standard test data
load_test_data <- function(sample = "NBE15") {
    if (sample == "NBE15") {
        cnv_file <- system.file("extdata", "NBE15",
                                "NBE15_comb_pro_extra2.51_1.txt",
                                package = "LACHESIS")
        vcf_file <- system.file("extdata", "NBE15",
                                "snvs_NBE15_somatic_snvs_conf_8_to_10.vcf",
                                package = "LACHESIS")
    } else if (sample == "NBE11") {
        cnv_file <- system.file("extdata", "NBE11",
                                "NBE11_comb_pro_extra2.59_0.83.txt",
                                package = "LACHESIS")
        vcf_file <- system.file("extdata", "NBE11",
                                "snvs_NBE11_somatic_snvs_conf_8_to_10.vcf",
                                package = "LACHESIS")
    } else if (sample == "NBE26") {
        cnv_file <- system.file("extdata", "NBE26",
                                "NBE26_comb_pro_extra3.25_0.88.txt",
                                package = "LACHESIS")
        vcf_file <- system.file("extdata", "NBE26",
                                "snvs_NBE26_somatic_snvs_conf_8_to_10.vcf",
                                package = "LACHESIS")
    }

    cnv <- readCNV(cn.info = cnv_file)
    snv <- readVCF(vcf = vcf_file, vcf.source = "dkfz")

    list(cnv = cnv, snv = snv)
}

# Test 1: Basic nbImport functionality
test_nbImport_basic <- function() {
    data <- load_test_data("NBE15")

    result <- nbImport(cnv = data$cnv, snv = data$snv,
                       purity = 1.0, ploidy = 2.51)

    expect_true(is.data.table(result),
                info = "nbImport should return a data.table")
    expect_true(nrow(result) > 0,
                info = "Result should contain variants")
    expect_true(all(c("snv_start", "cn_start", "cn_end", "TCN", "t_vaf") %in% names(result)),
                info = "Output should contain required columns")
}

# Test 2: Output attributes
test_nbImport_attributes <- function() {
    data <- load_test_data("NBE15")

    result <- nbImport(cnv = data$cnv, snv = data$snv,
                       purity = 0.8, ploidy = 2.51)

    expect_true(!is.null(attr(result, "purity")),
                info = "Should have purity attribute")
    expect_true(!is.null(attr(result, "ploidy")),
                info = "Should have ploidy attribute")
    expect_true(!is.null(attr(result, "cnv")),
                info = "Should have cnv attribute")
    expect_equal(attr(result, "purity"), 0.8,
                 info = "Purity attribute should match input")
}

# Test 3: Purity and ploidy coercion to numeric
test_nbImport_numeric_coercion <- function() {
    data <- load_test_data("NBE15")

    # Pass as strings
    result <- nbImport(cnv = data$cnv, snv = data$snv,
                       purity = "0.8", ploidy = "2.51")

    expect_true(is.numeric(attr(result, "purity")),
                info = "Purity should be coerced to numeric")
    expect_true(is.numeric(attr(result, "ploidy")),
                info = "Ploidy should be coerced to numeric")
}

# Test 4: Missing CNV input
test_nbImport_missing_cnv <- function() {
    data <- load_test_data("NBE15")

    expect_error(nbImport(cnv = NULL, snv = data$snv,
                          purity = 1.0, ploidy = 2.51),
                 info = "Missing CNV should raise error")
}

# Test 5: Missing SNV input
test_nbImport_missing_snv <- function() {
    data <- load_test_data("NBE15")

    expect_error(nbImport(cnv = data$cnv, snv = NULL,
                          purity = 1.0, ploidy = 2.51),
                 info = "Missing SNV should raise error")
}

# Test 6: Missing purity/ploidy
test_nbImport_missing_purity_ploidy <- function() {
    data <- load_test_data("NBE15")

    expect_error(nbImport(cnv = data$cnv, snv = data$snv,
                          purity = NULL, ploidy = 2.51),
                 info = "Missing purity should raise error")

    expect_error(nbImport(cnv = data$cnv, snv = data$snv,
                          purity = 1.0, ploidy = NULL),
                 info = "Missing ploidy should raise error")
}

# Test 7: No overlapping SNVs (should return empty or minimal result)
test_nbImport_no_overlap <- function() {
    # Create CNV with chromosome that doesn't match SNVs
    cnv <- data.table::data.table(
        chrom = 99, start = 1000, end = 2000,
        A = 1, B = 1, TCN = 2
    )
    snv <- data.table::data.table(
        chrom = 1, pos = 1500, ref = "A", alt = "T",
        t_ref_count = 50, t_alt_count = 50, t_depth = 100, t_vaf = 0.5
    )

    # nbImport handles non-overlapping data gracefully (returns result with no rows)
    result <- nbImport(cnv = cnv, snv = snv,
                       purity = 1.0, ploidy = 2.0)

    # Result should be a data.table (even if empty)
    expect_true(is.data.table(result),
                info = "nbImport should return data.table even with no overlap")
}

# Test 8: Invalid ref.build parameter
test_nbImport_invalid_ref_build <- function() {
    data <- load_test_data("NBE15")

    expect_error(nbImport(cnv = data$cnv, snv = data$snv,
                          purity = 1.0, ploidy = 2.51,
                          ref.build = "hg99"),
                 info = "Invalid ref.build should raise error")
}

# Test 9: Invalid assign.method parameter
test_nbImport_invalid_assign_method <- function() {
    data <- load_test_data("NBE15")

    expect_error(nbImport(cnv = data$cnv, snv = data$snv,
                          purity = 1.0, ploidy = 2.51,
                          sig.assign = TRUE,
                          assign.method = "invalid"),
                 info = "Invalid assign.method should raise error")
}

# Test 10: Signature assignment with file - "max" method
test_nbImport_sig_assign_max <- function() {
    data <- load_test_data("NBE15")
    sig_file <- system.file("extdata",
                            "NBE15_Decomposed_MutationType_Probabilities.txt",
                            package = "LACHESIS")

    result <- nbImport(cnv = data$cnv, snv = data$snv,
                       purity = 1.0, ploidy = 2.51,
                       sig.assign = TRUE, ID = "NBE15",
                       sig.file = sig_file, assign.method = "max")

    expect_true("Signature" %in% names(result),
                info = "Should have Signature column")
    expect_true("Probability" %in% names(result),
                info = "Should have Probability column")
    expect_true(nrow(result) > 0,
                info = "Should return variants")
}

# Test 11: Signature assignment with file - "sample" method
test_nbImport_sig_assign_sample <- function() {
    data <- load_test_data("NBE15")
    sig_file <- system.file("extdata",
                            "NBE15_Decomposed_MutationType_Probabilities.txt",
                            package = "LACHESIS")

    set.seed(123)
    result <- nbImport(cnv = data$cnv, snv = data$snv,
                       purity = 1.0, ploidy = 2.51,
                       sig.assign = TRUE, ID = "NBE15",
                       sig.file = sig_file, assign.method = "sample")

    expect_true("Signature" %in% names(result),
                info = "Should have Signature column")
    expect_true(nrow(result) > 0,
                info = "Should return variants")
}

# Test 12: Signature filtering with sig.select
test_nbImport_sig_select <- function() {
    data <- load_test_data("NBE15")
    sig_file <- system.file("extdata",
                            "NBE15_Decomposed_MutationType_Probabilities.txt",
                            package = "LACHESIS")

    result <- nbImport(cnv = data$cnv, snv = data$snv,
                       purity = 1.0, ploidy = 2.51,
                       sig.assign = TRUE, ID = "NBE15",
                       sig.file = sig_file,
                       sig.select = c("SBS1", "SBS5"))

    # Should only have selected signatures (plus NA for unmatched)
    expected_sigs <- c("SBS1", "SBS5", NA)
    unique_sigs <- unique(result$Signature)
    expect_true(all(unique_sigs %in% expected_sigs) || all(is.na(unique_sigs)),
                info = "Should only contain selected signatures")
}

# Test 13: Signature probability filtering with min.p
test_nbImport_sig_min_prob <- function() {
    data <- load_test_data("NBE15")
    sig_file <- system.file("extdata",
                            "NBE15_Decomposed_MutationType_Probabilities.txt",
                            package = "LACHESIS")

    result <- nbImport(cnv = data$cnv, snv = data$snv,
                       purity = 1.0, ploidy = 2.51,
                       sig.assign = TRUE, ID = "NBE15",
                       sig.file = sig_file, min.p = 0.5)

    # Check that probabilities meet threshold (excluding NA)
    valid_probs <- result$Probability[!is.na(result$Probability)]
    expect_true(all(valid_probs >= 0.5),
                info = "All probabilities should meet min.p threshold")
}

# Test 14: Signature colors attribute
test_nbImport_sig_colors_attribute <- function() {
    data <- load_test_data("NBE15")
    sig_file <- system.file("extdata",
                            "NBE15_Decomposed_MutationType_Probabilities.txt",
                            package = "LACHESIS")

    result <- nbImport(cnv = data$cnv, snv = data$snv,
                       purity = 1.0, ploidy = 2.51,
                       sig.assign = TRUE, ID = "NBE15",
                       sig.file = sig_file)

    expect_true(!is.null(attr(result, "sig.colors")),
                info = "Should have sig.colors attribute")
    expect_true(!is.null(names(attr(result, "sig.colors"))),
                info = "Signature colors should be named")
}

# Test 15: Multiple samples without signature assignment
test_nbImport_multiple_samples <- function() {
    data_nbe11 <- load_test_data("NBE11")
    data_nbe15 <- load_test_data("NBE15")

    result_11 <- nbImport(cnv = data_nbe11$cnv, snv = data_nbe11$snv,
                          purity = 0.83, ploidy = 2.59)
    result_15 <- nbImport(cnv = data_nbe15$cnv, snv = data_nbe15$snv,
                          purity = 1.0, ploidy = 2.51)

    expect_true(nrow(result_11) > 0 && nrow(result_15) > 0,
                info = "Should work with multiple samples")
}

# Test 16: Different reference genomes
test_nbImport_different_ref_builds <- function() {
    data <- load_test_data("NBE15")

    for (ref in c("hg19", "hg18", "hg38")) {
        result <- nbImport(cnv = data$cnv, snv = data$snv,
                           purity = 1.0, ploidy = 2.51,
                           ref.build = ref)
        expect_true(is.data.table(result),
                    info = paste("Should work with ref.build =", ref))
    }
}

# Test 17: Variant filtering via SNV input
test_nbImport_variant_filtering <- function() {
    data <- load_test_data("NBE15")

    # Use all variants
    result_all <- nbImport(cnv = data$cnv, snv = data$snv,
                           purity = 1.0, ploidy = 2.51)

    # Subset SNVs to high VAF only
    data$snv_filtered <- data$snv[t_vaf > 0.3]

    result_filtered <- nbImport(cnv = data$cnv, snv = data$snv_filtered,
                                purity = 1.0, ploidy = 2.51)

    # Filtered should have fewer or equal variants
    expect_true(nrow(result_filtered) <= nrow(result_all),
                info = "Filtered SNVs should result in fewer variants")
}

# Test 18: Copy number info in output
test_nbImport_cn_info_preserved <- function() {
    data <- load_test_data("NBE15")

    result <- nbImport(cnv = data$cnv, snv = data$snv,
                       purity = 1.0, ploidy = 2.51)

    expect_true("A" %in% names(result),
                info = "Should preserve A allele count")
    expect_true("B" %in% names(result),
                info = "Should preserve B allele count")
    expect_true("TCN" %in% names(result),
                info = "Should preserve total copy number")
}

# Run all tests in this file
for (test_func in ls(pattern = "^test_")) {
    get(test_func)()
}
