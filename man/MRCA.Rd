% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MRCA.R
\name{MRCA}
\alias{MRCA}
\title{Compute mutation densities at ECA and MRCA}
\usage{
MRCA(
  normObj = NULL,
  min.seg.size = 10^7,
  fp.mean = 0,
  fp.sd = 0,
  excl.chr = NULL
)
}
\arguments{
\item{normObj}{normalized clonal SNV counts stratified by copy number as generated by \code{\link{normalizeCounts}}.}

\item{min.seg.size}{the minimal segment length to be included in the quantification}

\item{fp.mean}{optional, the average false positive rate of clonal mutations (e.g., due to incomplete tissue sampling). Defaults to 0.}

\item{fp.sd}{optional, the standard deviation of the false positive rate of clonal mutations (e.g., due to incomplete tissue sampling). Defaults to 0.}

\item{excl.chr}{a vector of chromosomes that should be excluded from the quantification. e.g., due to reporter constructs in animal models.}
}
\value{
a data table reporting the assignment of individual segments to ECA or MRCA. Mutation densities at ECA and MRCA, and the bootstrapped 95\% CIs are stored as attributes. The columns in the data table report the following information:
\itemize{
\item{\code{chrom}}{Chromsoome}
\item{\code{TCN}}{Total copy number}
\item{\code{A}}{Number of A alleles}
\item{\code{B}}{Number of B alleles}
\item{\code{Seglength}}{Number of bps with the given copy number configuration on this chromosome}
\item{\code{n_mut_A}}{Normalized number of mutations present on all A alleles}
\item{\code{n_mut_B}}{Normalized number of mutations present on all B alleles}
\item{\code{n_mut_total}}{Normalized number of mutations per single copy of the segment}
\item{\code{density_A_mean}}{Normalized mean density of mutations present on all A alleles (1/Mb)}
\item{\code{density_B_mean}}{Normalized mean density of mutations present on all B alleles (1/Mb)}
\item{\code{density_total_mean}}{Normalized mean density of mutations per single copy of the segment (1/Mb)}
\item{\code{density_total_lower}}{Lower bound (95\% CI) of normalized density of mutations per single copy of the segment (1/Mb)}
\item{\code{density_total_upper}}{Upper bound (95\% CI) of normalized density of mutations per single copy of the segment (1/Mb)}
\item{\code{density_A_lower}}{Lower bound (95\% CI) of normalized density of mutations on all A alleles of the segment (1/Mb)}
\item{\code{density_A_upper}}{Upper bound (95\% CI) of normalized density of mutations on all A alleles of the segment (1/Mb)}
\item{\code{density_B_lower}}{Lower bound (95\% CI) of normalized density of mutations on all B alleles of the segment (1/Mb)}
\item{\code{density_B_upper}}{Upper bound (95\% CI) of normalized density of mutations on all B alleles of the segment (1/Mb)}
\item{\code{p_total_to_mrca}}{Probability that the density of mutations per single copy of the segment agrees with the mutation density at MRCA.}
\item{\code{p_A_to_to_mrca}}{Probability that the density of mutations on all A alleles of the segment agrees with the mutation density at MRCA.}
\item{\code{p_B_to_to_mrca}}{Probability that the density of mutations on all B alleles of the segment agrees with the mutation density at MRCA.}
\item{\code{p_adj_total_to_mrca}}{Probability that the density of mutations on all B alleles of the segment agrees with the mutation density at MRCA, adjusted for multiple sampling (Holm correction).}
\item{\code{p_adj_A_to_mrca}}{Probability that the density of mutations on all A alleles of the segment agrees with the mutation density at MRCA, adjusted for multiple sampling (Holm correction).}
\item{\code{p_adj_B_to_mrca}}{Probability that the density of mutations on all B alleles of the segment agrees with the mutation density at MRCA, adjusted for multiple sampling (Holm correction).}
\item{\code{MRCA_qual}}{Quality control. \code{PASS}, if the density of mutations on single copies agrees with the density at the MRCA.}
\item{\code{p_total_to_eca}}{Probability that the density of mutations per single copy of the segment agrees with the mutation density at ECA.}
\item{\code{p_A_to_to_eca}}{Probability that the density of mutations on all A alleles of the segment agrees with the mutation density at ECA.}
\item{\code{p_B_to_to_eca}}{Probability that the density of mutations on all B alleles of the segment agrees with the mutation density at ECA.}
\item{\code{p_adj_total_to_eca}}{Probability that the density of mutations on all B alleles of the segment agrees with the mutation density at ECA, adjusted for multiple sampling (Holm correction).}
\item{\code{p_adj_A_to_eca}}{Probability that the density of mutations on all A alleles of the segment agrees with the mutation density at ECA, adjusted for multiple sampling (Holm correction).}
\item{\code{p_adj_B_to_eca}}{Probability that the density of mutations on all B alleles of the segment agrees with the mutation density at ECA, adjusted for multiple sampling (Holm correction).}
\item{\code{A_time}}{Time of A allele gain (can be "ECA", "MRCA", "ECA/MRCA" if assignment is unclear, or "not mapped to ECA or MRCA" if density does not agree with either ECA or MRCA).}
\item{\code{B_time}}{Time of B allele gain (can be "ECA", "MRCA", "ECA/MRCA" if assignment is unclear, or "not mapped to ECA or MRCA" if density does not agree with either ECA or MRCA).}
}
snvs <- system.file("extdata", "NBE15", "snvs_NBE15_somatic_snvs_conf_8_to_10.vcf", package = "LACHESIS")
s_data <- readVCF(vcf = snvs, vcf.source = "dkfz")
aceseq_cn <- system.file("extdata", "NBE15", "NBE15_comb_pro_extra2.51_1.txt", package = "LACHESIS")
c_data <- readCNV(aceseq_cn)
nb <- nbImport(cnv = c_data, snv = s_data, purity = 1, ploidy = 2.51)
cl_muts <- clonalMutationCounter(nb)
norm_muts <- normalizeCounts(cl_muts)
mrca <- MRCA(norm_muts)
}
\description{
This function takes the normalized clonal mutation counts obtained with \code{\link{normalizeCounts}}to estimate mutation densities at MRCA and an earlier common ancestor, ECA.
}
