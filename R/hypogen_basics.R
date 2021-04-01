#' hypogen: Provides Hypoplectrus PopGen Data and Functions
#'
#' Basic reference data for population genomic studies in Caribbean
#' hamlets (Hypoplectrus sp.). It includes summary data of the Hamlet reference
#' genome (accession number: PRJEB27858, Hench et al. 2018) as well its annotation.
#' Furthermore it includes some functions to use the reference data  to put new data
#' into context.
#'
#' @docType package
#' @name hypogen
#'
#' @import dplyr
#' @import ggplot2
#' @import purrr
#' @import stringr
#' @import tibble
#' @importFrom rlang is_missing
NULL

.onAttach <- function(libname, pkgname) {
  # cat("--- Welcome to", crayon::blue("hypogen"), "---\n")
}
#' Shaded LG backgrounds
#'
#' \code{geom_hypo_LG} adds the linkage group (LG) boundaries to background of a ggplot.
#'
#' To add context to genome wide plots of a statistic (e.g. Fst) this function
#' adds a alternating shaded background to the plot.
#'
#' @param ... catch all parameter (legacy)
#'
#' @seealso \code{\link{scale_fill_hypo_LG_bg}}
#'
#' @examples
#' ggplot() +
#'   geom_hypo_LG() +
#'   scale_fill_hypo_LG_bg() +
#'   scale_x_hypo_LG() +
#'   coord_cartesian(ylim = c(0, 1)) +
#'   theme_hypo()
#'
#' @export
geom_hypo_LG <- function(...){
  geom_hypobg(inherit.aes = FALSE,
              data = hypogen::hypo_karyotype,
            aes(xmin = GSTART,
                xmax = GEND,
                ymin = -Inf,
                ymax = Inf,
                hypobg = GROUP))
}

#' Format genome wide x axis
#'
#' \code{scale_x_hypo_LG} adjusts the x axis to match the hamet genome.
#'
#' Each of the LGs is labelled with on tick which is placed in the
#' center of the LG.
#' The should be combined with \code{hypogen::geom_hypo_LG}
#'
#' @param ...      parameters passed to ggplot2::scale_x_continuous()
#' @param name     string/ expression, x axis title
#' @param expand   numeric vector of length 2
#' @param breaks   numeric vector, x axis breaks
#' @param labels   string vector, x axis labels
#' @param position string ("top"/ "bottom"), axis placement
#'
#' @seealso \code{\link{geom_hypo_LG}}, \code{\link{scale_color_hypo_LG}}
#'
#' @examples
#' ggplot() +
#'   geom_hypo_LG() +
#'   scale_fill_hypo_LG_bg() +
#'   scale_x_hypo_LG() +
#'   coord_cartesian(ylim = c(0, 1)) +
#'   theme_hypo()
#'
#' @export
scale_x_hypo_LG <- function(..., name = '', expand = c(0, 0),
                            breaks = (hypo_karyotype$GSTART + hypo_karyotype$GEND) / 2,
                            labels = 1:24, position = "top"){
  ggplot2::scale_x_continuous(name = name, expand = expand,breaks = breaks,
                     labels = labels, position = position, ...)
}

#' Set LG color palette
#'
#' \code{scale_color_hypo_LG} adjusts the color scale with default color for each LG.
#'
#' Each of the LGs is assigned to one default color.
#'
#' \strong{! CAREFUL:} should not be used if additional color levels exiest.
#'
#' @param ...  parameters passed to ggplot2::scale_color_manual()
#' @param name string/ expression, color guide title
#'
#' @seealso \code{\link{geom_hypo_LG}}, \code{\link{scale_x_hypo_LG}}
#'
#' @export
scale_color_hypo_LG <- function(..., name = 'LG'){
  ggplot2::scale_color_manual(..., name = name, values = hypo_clr_LGs)
}

#' Set LG fill
#'
#' \code{scale_fill_hypo_LG_bg} sets the color scheme for the LG background shading.
#'
#' Fill values are set for tnhe two LG classes (even/odd LG nr).
#'
#' \strong{! CAREFUL:} LG fill values are set to "zz1" $ zz2.
#' If additional fill levels exist, the fill values of the LGs should be appended
#' to the desired fill scale.
#'
#' @param ...    parameters passed to scale_hypobg_manual()
#' @param values named string vector `c(odd = ..., even = ...)`, color scheme for alternating linkage group background
#'
#' @seealso \code{\link{geom_hypo_LG}}, \code{\link{scale_x_hypo_LG}}, \code{\link{scale_color_hypo_LG}}

#' @examples
#' ggplot() +
#'   geom_hypo_LG() +
#'   scale_fill_hypo_LG_bg(values = c(even = "#274263", odd = "#C09C60")) +
#'   scale_x_hypo_LG() +
#'   coord_cartesian(ylim = c(0, 1)) +
#'   theme_hypo()
#'
#' @export
scale_fill_hypo_LG_bg <- function(..., values = c(odd = NA, even = hypo_clr_lg)){
  scale_hypobg_manual(..., values = values, guide = F)
}

#' Set hypogen plot theme
#'
#' \code{theme_hypo} sets the default hypogen plot theme.
#'
#' This theme is optimized for genome wide plots.
#'
#' @param ... parameters passed to ggplot2::theme()
#'
#' @seealso \code{\link{theme_hypo_anno_extra}}
#'
#' @export
theme_hypo <- function (...) {
  ggplot2::theme_bw(base_size = 10, base_family = "Helvetica") %+replace%
    ggplot2::theme(...,
          plot.background = element_blank(),
          panel.background  = element_blank(),
          panel.grid=element_blank(),
          panel.border = element_blank(),
          axis.title.x = element_blank(),
          axis.line = element_line(),
          strip.background = element_rect(fill = NA, color = hypo_clr_lg),
          legend.background = element_rect(fill = "transparent", color = NA),
          legend.key = element_rect(fill = "transparent", color = NA)
    )
}

#' Set hypogen plot theme
#'
#' \code{theme_hypo_anno_extra} adjusts the theme for annotation plots.
#'
#' @param ... catch all parameter (legacy)
#'
#' @seealso \code{\link{theme_hypo}}
#' @export
theme_hypo_anno_extra <- function(...){
  ggplot2::theme(axis.title.x = element_text(),
        panel.grid.major.x = element_line(color = hypo_clr_lg),
        panel.grid.minor.x = element_line(color = hypo_clr_lg),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.placement = 'outside')
}


#' Determining the LG from the genomic position
#'
#' \code{hypo_which_CHROM_s} determines the LG from the genomic position.
#'
#' The function finds the respective LG on which a specific genomic position (1 - 559,649,677)
#' resides.
#'
#' @param x integer scalar (mandatory), a bp position on the hamlet genome.
#'
#' @seealso \code{\link{hypo_which_CHROM}}
#'
#' @examples
#' hypo_which_CHROM_s(123581321)
#'
#' @export
hypo_which_CHROM_s <- function(x){
  CHROM <- hypo_karyotype$CHROM[x >= hypo_karyotype$GSTART & x <= hypo_karyotype$GEND]
  CHROM
}

#' Determining the LG from the genomic positions
#'
#' \code{hypo_which_CHROM} determines the LG from a vector of genomic position.
#'
#' The function vectorizes the \code{hypo_which_CHROM_s} function.
#'
#' @param x integer vector (mandatory), a vector of bp positions on the hamlet genome.
#'
#' @seealso \code{\link{hypo_which_CHROM_s}}
#'
#' @examples
#'
#' hypo_which_CHROM(c(11235813, 123581321))
#'
#' @export
hypo_which_CHROM <- function(x){
  purrr::map_chr(x,hypo_which_CHROM_s)
}
