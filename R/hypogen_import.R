#' Import SNP based statistics
#'
#' \code{hypo_import_snps} imports SNP based statistics (such as VCFtools output).
#'
#' This function imports SNP based statistics. The input is
#' expected to be tab separated, to contain a CHROM and POS column
#' and can optionally be gz compressed.
#' After import, the CHROM based position is transposed to a genome
#' wide position for continuous visualization.
#'
#' @param file string skalar (mandatory), the input file
#' @param gz logical skalar (optional), is the input file gz compressed?
#'
#' @seealso \code{\link{hypo_import_windows}}
#'
#' @examples
#'
#' file_snps <- system.file("extdata", "example.weir.fst.gz", package = "hypogen")
#'
#' hypo_import_snps(file = file_snps, gz = TRUE)
#'
#' @export
hypo_import_snps <- function(file, gz=FALSE,...){
  if(gz){
    import <- function(...,file){
      readr::read_delim(file = gzfile(file), delim='\t', ...)
      }
    } else {
      import <- function(...,file){
        readr::read_delim(file = file, delim='\t', ...)
      }
    }

  data <- import(file=file,...) %>%
    dplyr::left_join(hypo_chrom_start,by="CHROM") %>%
    dplyr::mutate(GPOS = GSTART + POS)
  data

}

#' Import window based statistics
#'
#' \code{hypo_import_windows} imports window based statistics (such as VCFtools output).
#'
#' This function imports window based statistics. The input is
#' expected to be tab separated, to contain a CHROM, BIN_START
#' and BIN END column and can optionally be gz compressed.
#' After import, the CHROM based position (center of the window)
#' is transposed to a genome wide position for continuous
#' visualization.
#'
#' @param file string skalar (mandatory), the input file
#' @param gz logical skalar (optional), is the input file gz compressed?
#'
#' @seealso \code{\link{hypo_import_snps}}
#'
#' @examples
#'
#' file_windowed <- system.file("extdata", "example.windowed.weir.fst.gz", package = "hypogen")
#'
#' hypo_import_windows(file = file_windowed, gz = TRUE)
#'
#' @export
hypo_import_windows <- function(file, gz=FALSE,...){
  if(gz){
    import <- function(...,file){
      readr::read_delim(file = gzfile(file), delim='\t', ...)
    }
  } else {
    import <- function(...,file){
      readr::read_delim(file = file, delim='\t', ...)
    }
  }

  data <- import(file=file,...) %>%
    dplyr::left_join(hypo_chrom_start,by="CHROM") %>%
    dplyr::mutate(POS = (BIN_START + BIN_END)/2,
                  GPOS = GSTART + POS)
  data

}
