# Load required packages
library(LACHESIS)
library(data.table)

# Test readCNV function

# Test 1: Read ACESeq CNV file
test_readCNV_aceseq <- function() {
    aceseq_cn <- system.file("extdata", "ACESeq/NBE11_comb_pro_extra2.59_0.83.txt",
        package = "LACHESIS"
    )
    result <- readCNV(cn.info = aceseq_cn)

    expect_true(is.data.table(result),
        info = "readCNV should return a data.table"
    )
    expect_true(nrow(result) > 0,
        info = "ACESeq CNV should contain segments"
    )
    expect_true(all(c("Chr", "Start", "End", "TCN", "A", "B") %in% names(result)),
        info = "Output should contain required columns"
    )
}

# Test 2: Read ASCAT CNV file
test_readCNV_ascat <- function() {
    ascat_cn <- system.file("extdata", "ASCAT/S98.segments.txt",
        package = "LACHESIS"
    )
    result <- readCNV(cn.info = ascat_cn)

    expect_true(is.data.table(result),
        info = "readCNV should return a data.table"
    )
    expect_true(nrow(result) > 0,
        info = "ASCAT CNV should contain segments"
    )
}

# Test 3: Read PURPLE CNV file (NEW format)
test_readCNV_purple <- function() {
    purple_cn <- system.file("extdata", "PURPLE/NB-S-599-T.purple.cnv.somatic.tsv",
        package = "LACHESIS"
    )
    result <- readCNV(cn.info = purple_cn)

    expect_true(is.data.table(result),
        info = "readCNV should return a data.table"
    )
    expect_true(nrow(result) > 0,
        info = "PURPLE CNV should contain segments"
    )
}

# Test 4: Output column types
test_readCNV_column_types <- function() {
    aceseq_cn <- system.file("extdata", "ACESeq/NBE11_comb_pro_extra2.59_0.83.txt",
        package = "LACHESIS"
    )
    result <- readCNV(cn.info = aceseq_cn)

    expect_true(is.character(result$Chr) || is.numeric(result$Chr),
        info = "Chr should be character or numeric"
    )
    expect_true(is.integer(result$Start) || is.numeric(result$Start),
        info = "Start should be numeric"
    )
    expect_true(is.integer(result$End) || is.numeric(result$End),
        info = "End should be numeric"
    )
    expect_true(is.integer(result$TCN) || is.numeric(result$TCN),
        info = "TCN should be numeric"
    )
    expect_true(is.integer(result$A) || is.numeric(result$A),
        info = "A should be numeric"
    )
    expect_true(is.integer(result$B) || is.numeric(result$B),
        info = "B should be numeric"
    )
}

# Test 5: Segment length validation
test_readCNV_segment_length <- function() {
    aceseq_cn <- system.file("extdata", "ACESeq/NBE11_comb_pro_extra2.59_0.83.txt",
        package = "LACHESIS"
    )
    result <- readCNV(cn.info = aceseq_cn)

    # End should be >= start
    expect_true(all(result$End >= result$Start),
        info = "Segment End should be >= Start"
    )

    # Segment length should be positive
    expect_true(all((result$End - result$Start) > 0),
        info = "Segment length should be positive"
    )
}

# Test 6: Copy number range filtering
test_readCNV_cn_filtering <- function() {
    aceseq_cn <- system.file("extdata", "ACESeq/NBE11_comb_pro_extra2.59_0.83.txt",
        package = "LACHESIS"
    )

    # Strict copy number range
    result_strict <- readCNV(cn.info = aceseq_cn, max.cn = 3)

    # Lenient copy number range
    result_lenient <- readCNV(cn.info = aceseq_cn, max.cn = 5)

    # Lenient should have more or equal segments
    expect_true(nrow(result_lenient) >= nrow(result_strict),
        info = "Lenient CN filtering should have more/equal segments"
    )

    # Check that all TCN values in strict are <= 3
    expect_true(all(result_strict$TCN <= 3),
        info = "All segments should have TCN <= max.cn threshold"
    )
}

# Test 7: Ignore XY chromosomes
test_readCNV_ignore_XY <- function() {
    aceseq_cn <- system.file("extdata", "ACESeq/NBE11_comb_pro_extra2.59_0.83.txt",
        package = "LACHESIS"
    )

    result_with_XY <- readCNV(cn.info = aceseq_cn, ignore.XY = FALSE)
    result_no_XY <- readCNV(cn.info = aceseq_cn, ignore.XY = TRUE)

    # Should have fewer or equal segments when ignoring XY
    expect_true(nrow(result_no_XY) <= nrow(result_with_XY),
        info = "Ignoring XY should not increase segment count"
    )

    # Should not have X or Y chromosomes
    expect_true(!any(result_no_XY$Chr %in% c("X", "Y", 23, 24)),
        info = "Should not contain X or Y chromosomes"
    )
}

# Test 8: Allele ratio validation
test_readCNV_allele_ratios <- function() {
    aceseq_cn <- system.file("extdata", "ACESeq/NBE11_comb_pro_extra2.59_0.83.txt",
        package = "LACHESIS"
    )
    result <- readCNV(cn.info = aceseq_cn)

    # A + B should equal TCN
    expect_true(all(result$A + result$B == result$TCN),
        info = "A + B should equal TCN"
    )
}

# Test 9: Merge tolerance (adjacent segments)
test_readCNV_merge_tolerance <- function() {
    aceseq_cn <- system.file("extdata", "ACESeq/NBE11_comb_pro_extra2.59_0.83.txt",
        package = "LACHESIS"
    )

    # No merging
    result_no_merge <- readCNV(cn.info = aceseq_cn, merge.tolerance = 0)

    # With merging
    result_merge <- readCNV(cn.info = aceseq_cn, merge.tolerance = 1e6)

    # Merging should reduce or maintain segment count
    expect_true(nrow(result_merge) <= nrow(result_no_merge),
        info = "Merging should reduce or maintain segment count"
    )
}

# Test 10: Non-existent file
test_readCNV_nonexistent_file <- function() {
    expect_error(readCNV(cn.info = "/nonexistent/path/file.txt"),
        info = "Non-existent file should raise error"
    )
}

# Test 11: Multiple samples
test_readCNV_multiple_formats <- function() {
    # Test that different CNV formats work
    aceseq_cn <- system.file("extdata", "ACESeq/NBE11_comb_pro_extra2.59_0.83.txt",
        package = "LACHESIS"
    )
    ascat_cn <- system.file("extdata", "ASCAT/S98.segments.txt",
        package = "LACHESIS"
    )

    result_aceseq <- readCNV(cn.info = aceseq_cn)
    result_ascat <- readCNV(cn.info = ascat_cn)

    expect_true(nrow(result_aceseq) > 0 && nrow(result_ascat) > 0,
        info = "Should read both ACESeq and ASCAT formats"
    )
}

# Test 12: Tumor ID parameter
test_readCNV_tumor_id <- function() {
    aceseq_cn <- system.file("extdata", "ACESeq/NBE11_comb_pro_extra2.59_0.83.txt",
        package = "LACHESIS"
    )
    result <- readCNV(cn.info = aceseq_cn, tumor.id = "NBE11")

    expect_true(is.data.table(result),
        info = "Should work with tumor.id parameter"
    )
}

# Run all tests in this file
for (test_func in ls(pattern = "^test_")) {
    get(test_func)()
}
