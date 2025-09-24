# Test readCNV function

# Test parameter validation
expect_error(readCNV(), "Missing cn.info!")

# Test with example data files
aceseq_cn <- system.file("extdata", "ACESeq/NBE11_comb_pro_extra2.59_0.83.txt", package = "LACHESIS")
ascat_cn <- system.file("extdata", "ASCAT/S98.segments.txt", package = "LACHESIS")
nb15_cn <- system.file("extdata", "NBE15", "NBE15_comb_pro_extra2.51_1.txt", package = "LACHESIS")

# Test ACESeq CNV if file exists
if (file.exists(aceseq_cn)) {
  aceseq_data <- readCNV(aceseq_cn)
  expect_true(is.data.frame(aceseq_data))
  expect_true(all(c("Chr", "Start", "End", "A", "B", "TCN") %in% colnames(aceseq_data)))
  expect_true(all(aceseq_data$TCN >= 0))
  expect_true(all(aceseq_data$A >= 0))
  expect_true(all(aceseq_data$B >= 0))
  expect_true(all(aceseq_data$TCN == aceseq_data$A + aceseq_data$B))
}

# Test ASCAT CNV if file exists
if (file.exists(ascat_cn)) {
  ascat_data <- readCNV(ascat_cn)
  expect_true(is.data.frame(ascat_data))
  expect_true(all(c("Chr", "Start", "End", "A", "B", "TCN") %in% colnames(ascat_data)))
}

# Test NBE15 CNV if file exists
if (file.exists(nb15_cn)) {
  nb15_data <- readCNV(nb15_cn)
  expect_true(is.data.frame(nb15_data))
  expect_true(all(c("Chr", "Start", "End", "A", "B", "TCN") %in% colnames(nb15_data)))
}

# Test max copy number filtering
if (file.exists(aceseq_cn)) {
  filtered_data <- readCNV(aceseq_cn, max.cn = 2)
  expect_true(all(filtered_data$TCN <= 2))
}

# Test ignore XY parameter
if (file.exists(aceseq_cn)) {
  data_with_xy <- readCNV(aceseq_cn, ignore.XY = FALSE)
  data_without_xy <- readCNV(aceseq_cn, ignore.XY = TRUE)

  # Should have fewer or equal rows when ignoring XY
  expect_true(nrow(data_without_xy) <= nrow(data_with_xy))
}