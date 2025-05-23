% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nbImport.R
\name{nbImport}
\alias{nbImport}
\title{Combine CNVs and SNVs}
\usage{
nbImport(cnv = NULL, snv = NULL, purity = NULL, ploidy = NULL)
}
\arguments{
\item{cnv}{CNV data from \code{\link{readCNV}}}

\item{snv}{SNV data from \code{\link{readVCF}}}

\item{purity}{tumor cell content}

\item{ploidy}{average copy number in the tumor sample}
}
\value{
a data.table
}
\description{
Merges CNVs and SNVs into a single data.table. Each variant is assigned to its corresponding copy number segment and status.
}
\examples{
snvs <- system.file("extdata", "NBE15", "snvs_NBE15_somatic_snvs_conf_8_to_10.vcf", package = "LACHESIS")
s_data <- readVCF(vcf = snvs, vcf.source = "dkfz")
aceseq_cn <- system.file("extdata", "NBE15", "NBE15_comb_pro_extra2.51_1.txt", package = "LACHESIS")
c_data <- readCNV(aceseq_cn)
nb <- nbImport(cnv = c_data, snv = s_data, purity = 1, ploidy = 2.51)
}
\seealso{
\code{\link{plotNB}}
}
