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
#' @param file string scalar (mandatory), the input file
#' @param gz   logical scalar (optional), is the input file gz compressed?
#' @param run  string scalar (optional), appends a column RUN as ID if several files should be merged later
#' @param ...  parameters passed to vroom::vroom()
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
hypo_import_snps <- function(file, gz=FALSE, run, ...){
  if(gz){
    import <- function(...,file){
      vroom::vroom(file = gzfile(file), delim='\t', ...)
      }
    } else {
      import <- function(...,file){
        vroom::vroom(file = file, delim='\t', ...)
      }
    }

  data <- import(file = file,...) %>%
    dplyr::left_join(hypo_chrom_start,by="CHROM") %>%
    dplyr::mutate(GPOS = GSTART + POS)

  if(!missing("run")){
    stopifnot(length(run) == 1)
    stopifnot(is.character(run))

    data <- data %>% mutate(RUN = run)
  }
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
#' @param file string scalar (mandatory), the input file
#' @param gz   logical scalar (optional), is the input file gz compressed?
#' @param run  string scalar (optional), appends a column RUN as ID if several files should be merged later
#' @param ...  parameters passed to vroom::vroom()
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
hypo_import_windows <- function(file, gz=FALSE, run,...){
  if(gz){
    import <- function(...,file){
      vroom::vroom(file = gzfile(file), delim='\t', ...)
    }
  } else {
    import <- function(...,file){
      vroom::vroom(file = file, delim='\t', ...)
    }
  }

  data <- import(file=file,...) %>%
    dplyr::left_join(hypo_chrom_start,by="CHROM") %>%
    dplyr::mutate(POS = (BIN_START + BIN_END)/2,
                  GPOS = GSTART + POS)

  if(!missing("run")){
    stopifnot(length(run) == 1)
    stopifnot(is.character(run))

    data <- data %>% mutate(RUN = run)
  }
  data
}

#' Import genotype frequencies
#'
#' \code{hypo_import_genotype_freq} takes a genotype table (such as VCFtools --012 output) and computes genotype frequencies.
#'
#' This function takes a genotype table (a transposed VCFtools --012
#' output with samples as columns and SNPs as rows) and computes
#' genotype frequencies from it. Genotypes are assumed to be encoded as:
#' \itemize{
#'   \item{0: Homozygous (reference allele)}
#'   \item{1: Heterozygous}
#'   \item{2: Homozygous (alternative allele)}
#' }
#'
#' The genotype frequencies can later be visualized
#' using the {ggtern} package.
#'
#' This approach assumes biallelic SNPs!
#'
#' @param AA        string scalar ('ref' or 'major', optional), should the reference or the major allele be encoded as A?
#' @param delim     string scalar (optional), delimiter of the input file
#' @param file_path string scalar (mandatory), the input file
#'
#' @export
hypo_import_genotype_freq <- function(file_path, AA = 'ref', delim = ' '){
  stopifnot(AA %in% c('ref','major'))

  df <- vroom::vroom(file_path, delim = delim)

  if(AA == 'ref')  {
    df <- df %>% dplyr::mutate(AAc = rowSums(.==0,na.rm = TRUE),
                Aac = rowSums(.==1,na.rm = TRUE),
                aac = rowSums(.==2,na.rm = TRUE),
                n = rowSums(!is.na(.)),
                AA = AAc/n, Aa = Aac/n, aa = aac/n) %>%
    dplyr::select(AA,aa,Aa,n)
    return(df)
  } else if(AA == 'major') {
    df <- df %>% dplyr::mutate(AAc = rowSums(.==0,na.rm = TRUE),
                  Aac = rowSums(.==1,na.rm = TRUE),
                  aac = rowSums(.==2,na.rm = TRUE),
                  n = rowSums(!is.na(.)),
                  AAp = AAc/n, Aa = Aac/n, aap = aac/n) %>%
      dplyr::select(AAp,aap,Aa,n) %>%
      dplyr::mutate(AA = ifelse(AAp>aap,AAp,aap),
             aa = ifelse(AAp<aap,AAp,aap))%>%
      dplyr::select(AA,aa,Aa,n)
    return(df)
  }
}


#' Compute Hardy-Weinberg genopye frequencies
#'
#' \code{hypo_hwe} computes Hardy-Weinberg genopye frequencies for p ranging from 0 to 1.
#'
#' This function creates a range of allele frequencies ranging from p = 0 to p = 1
#' and computes the respective genotype frequencies according to Hardy-Weinberg:
#'
#' \itemize{
#'   \item{\eqn{q = 1 - p}}
#'   \item{\eqn{AA = p^2}}
#'   \item{\eqn{Aa = 2pq}}
#'   \item{\eqn{aa = q^2}}
#' }
#'
#' The function returns the results as table.
#'
#' @param n integer scalar (>= 2, mandatory), length of the
#'    allele frequency range
#'
#' @examples
#'
#' genotype_hwe <- hypo_hwe(100)
#'
#' @export
hypo_hwe <- function(n = 100L){
  stopifnot(is.numeric(n) & n %%1 == 0)
  stopifnot(n > 2)

  tibble::tibble(p=seq(0,1,length.out = n),
                 q=1-p,
                 AA = p^2,
                 Aa = 2*p*q,
                 aa = q^2)
  }
