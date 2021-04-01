#' Hamlet Genome Linkage Group Start Positions.
#'
#' A dataset containing the genomic start positions of the linkage groups of the Hamlet reference genome.
#'
#' @format A tibble with 24 rows and 2 variables:
#' \describe{
#'   \item{CHROM}{string, linkage group ID, eg ("LG01")}
#'   \item{GSTART}{integer, genomic start positions of the linkage groups (bp)}
#' }
#' @source \url{https://doi.org/10.1038/s41559-019-0814-5}
"hypo_chrom_start"

#' Hamlet Genome Color Scheme.
#'
#' A color scheme for the hamlet genome consisting of 24 colors for the 24 linkage groups.
"hypo_clr_LGs"

#' Hamlet Genome LG Background Color.
#'
#' A default background color for "even" linkage groups of the hamlet genome.
"hypo_clr_lg"

#' Hamlet Genome Karyotype.
#'
#' A dataset containing the hamlet karyotype according to the first reference genome.
#'
#' @format A tibble with 24 rows and 5 variables:
#' \describe{
#'   \item{CHROM}{string, linkage group ID, eg ("LG01")}
#'   \item{LENGTH}{integer, length of the linkage groups (bp)}
#'   \item{GSTART}{integer, genomic start positions of the linkage groups (bp)}
#'   \item{GEND}{integer, genomic end positions of the linkage groups (bp)}
#'   \item{GROUP}{string ("odd"/"even"), alternating across linagke groups}
#' }
#' @source \url{https://doi.org/10.1038/s41559-019-0814-5}
"hypo_karyotype"
