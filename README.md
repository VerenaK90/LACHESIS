# &Lambda;&Alpha;&Xi;&Epsilon;&Sigma;&Iota;&Sigma; (LACHESIS)

<img align="left" src="man/figures/LACHESIS_logo.svg" alt="LACHESIS logo" width="180"></img>

[![R-CMD-check](https://github.com/VerenaK90/LACHESIS/workflows/R-CMD-check/badge.svg)](https://github.com/VerenaK90/LACHESIS/actions)
[![Bioc-check](https://github.com/VerenaK90/LACHESIS/workflows/Bioconductor-check/badge.svg)](https://github.com/VerenaK90/LACHESIS/actions)

&Lambda;&Alpha;&Xi;&Epsilon;&Sigma;&Iota;&Sigma; (LACHESIS) is an R package to infer tumor evolution during malignant tranformation from bulk whole genome sequencing data. It uses single nucleotide variants and (ideally allele-specific) copy number information as input to compute mutation densities at clonal chromosomal gains and at the tumor's most recent common ancestor. Additionally, &Lambda;&Alpha;&Xi;&Epsilon;&Sigma;&Iota;&Sigma; tests whether there is evidence for an early common ancestor, oncogenic events acquired prior to malignant transformation. &Lambda;&Alpha;&Xi;&Epsilon;&Sigma;&Iota;&Sigma; also provides modalities to compare tumor evolution across patient collectives, and to correlate evolutionary timings with outcome. Please refer to our [vignette](/vignettes/vignette_LACHESIS.Rmd) for further information. 

<br clear="all"/>

## Documentation

A [vignette](vignettes/vignette_LACHESIS.Rmd) explaining how to use &Lambda;&Alpha;&Xi;&Epsilon;&Sigma;&Iota;&Sigma;, input and output formats alongside with example data can be installed by setting `build_vignettes = TRUE` and `dependencies = TRUE`.

## Demo 

Please refer to our [vignette](vignettes/vignette_LACHESIS.Rmd) for an example analysis with &Lambda;&Alpha;&Xi;&Epsilon;&Sigma;&Iota;&Sigma;.

## System requirements

### Hardware requirements

&Lambda;&Alpha;&Xi;&Epsilon;&Sigma;&Iota;&Sigma; runs on a standard computer.

### Software requirements

#### OS Requirements

This package is supported for macOS, Windows and Linux. The package has been tested on the following systems:

MacBook Pro: Sonoma (14.6.1)
Windows 10, Version 22H2

#### R dependencies

&Lambda;&Alpha;&Xi;&Epsilon;&Sigma;&Iota;&Sigma; was tested on R v4.5.1 and requires installation of the packages data.table (1.15.4), ggplot (≥3.5.1), tidyr (1.3.1), vcfR (1.15.0).

## Version

0.99.0

## Citation

Körber et al., Neuroblastoma arises in early fetal development and its evolutionary duration predicts outcome, *Nature Genetics* **55**:619-630 (2023).

## Installation

```r
devtools::install_github("VerenaK90/LACHESIS")
````

To install the vignette, run

```r
devtools::install_github("VerenaK90/LACHESIS", build_vignettes = TRUE, dependencies = TRUE)
````

The vignette can then be viewed by typing `vignette("vignette_LACHESIS", "LACHESIS")`. The package including its vignette can be installed within a few minutes.

## License

&Lambda;&Alpha;&Xi;&Epsilon;&Sigma;&Iota;&Sigma; is run under an **GPL (>= 3)** License.

## Contact

Verena Körber (verena.korber@ndcls.ox.ac.uk), Anand Mayakonda (a.mayakonda@kitz-heidelberg.de), Maximilia Eggle (maximilia.eggle@kitz-heidelberg.de).
