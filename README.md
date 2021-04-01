# hypogen <img src="man/figures/logo.png" align="right" alt="" width="120" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/k-hench/hypogen/workflows/R-CMD-check/badge.svg)](https://github.com/k-hench/hypogen/actions)
<!-- badges: end -->

The **R** package **hypogen** provides resources for population genetic analysis based on the hamlet reference genome.

The genome was published in our **Inter-chromosomal coupling between vision and pigmentation genes during genomic divergence** in [*Nature Ecology & Evolotion*](https://www.nature.com/natecolevol/) (DOI: 10.1038/s41559-019-0814-5).

The original genome can be downloaded from the [European Nucleotide Archive](https://www.ebi.ac.uk/ena/data/view/GCA_900610375.1) (project accession number PRJEB27858).

The annotation is also deposited at [dryad](https://doi.org/10.5061/dryad.pg8q56g).

Please also look at the [documentation](https://k-hench.github.io/hypogen/index.html).

## Dependencies 

**hypogen** depends on a bioconductor (non-CRAN) R-package.
To be able to install the package successfully, the following package will also need to be installed:

```r
install.packages("remotes")
remotes::install_bioc("rtracklayer")
```

## Install

To install **hypogen** please run:
```
remotes::install_github("k-hench/hypogen")
```
