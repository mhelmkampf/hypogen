#' Loading the annotation dataframe
#'
#' \code{hypo_annotation_get} loads the genome annotation.
#'
#' Hypogen comes with the annotation data of the halmet genome.
#' To load the data into R \code{hypo_annotation_get} makes use of
#' the \strong{rtracklayer} package.
#' The function returns a list containg two data frames: the first one
#' holds the general mRNA extends, the second ond the detailed Exon boundaries.
#'
#' Only the data of a specified range on a specific LG is loaded.
#' Genes of interes as well as genes of minor interest can be specified and
#' the number of lines that the mRNAs can be specified to avoid overlapping
#' when plotting the annotations.
#'
#' @param searchLG string skalar (mandatory), should be on of "LG01" - "LG24"
#' @param xrange integer vector (mandatory), data range to be loaded(start bp - end bp).
#'   Positions are defined with respect to the LG, NOT the overall genomic position.
#' @param genes_of_interest string vector (optional), tags specific gene names for highlighting when plotting.
#'   (needs to exactly match the gene name of the original gff file)
#' @param genes_of_sec_interest string vector (optional), tags specific gene names for secondary highlighting
#'  when plotting. (needs to exactly match the gene name of the original gff file)
#' @param anno_rown integer skalar (optional), defines the number of gene rows to avoid overlapping.
#'
#' @seealso \code{\link{hypo_annotation_baseplot}}
#'
#' @export
hypo_annotation_get <- function(searchLG, xrange, genes_of_interest=c(),
                                    genes_of_sec_interest=c(), anno_rown = 3){
  gfffile <- system.file("extdata", "HP.annotation.named.gff.gz", package="hypogen")
  gff_filter <- list(seqid = searchLG)
  data <- as.data.frame(rtracklayer::readGFF(gzfile(gfffile), filter = gff_filter)) %>%
    dplyr::mutate(Parent = as.character(Parent))

  mRNAs <- data %>%
    dplyr::filter(type == 'mRNA', end > xrange[1], start < xrange[2]) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(yl = row_number() %% anno_rown+2) %>%
    dplyr::rowwise()%>%
    dplyr::mutate(checkStart =ifelse(start < xrange[1], -Inf, start),
           checkEnd = ifelse(end > xrange[2], Inf, end),
           ps = ifelse(strand == '-', checkEnd, checkStart),
           pe = ifelse(strand == '-', checkStart, checkEnd),
           labelx = mean(c(sort(c(xrange[1],ps))[2],
                         sort(c(xrange[2],pe))[1])),
           window = 'bold(Gene)',
           clr = ifelse(Parentgenename %in% genes_of_interest, "zzz3",
                      ifelse(Parentgenename %in% genes_of_sec_interest, "zzz2", "zzz1")),
           label = unlist(strsplit(tolower(Parentgenename), "_"))[1]) %>%
    dplyr::select(seqid, type, start, end, strand, Name, ID,
           oId, Parentgene, Parentgenename,  yl, checkStart, checkEnd,
          ps, pe, labelx, window, clr,label);

  names(mRNAs)[names(mRNAs)=='ID'] <- 'Parent'

  exons <- data %>%
    dplyr::filter(type=='exon',end>xrange[1],start<xrange[2]) %>%
    merge(.,mRNAs %>% dplyr::select(Parent,yl,clr),by='Parent',all.x=T) %>%
    dplyr::mutate(ps=ifelse(strand=='-',end,start),
           pe=ifelse(strand=='-',start,end),
           window='bold(Gene)') %>%
    dplyr::select(Parent, seqid, type, start, end,  strand,
           Name, ID, oId, Parentgene, Parentgenename, yl,
           clr, ps, pe, window)


  return(list(mRNAs,exons))}

#' Plotting the annotation dataframe
#'
#' \code{hypo_annotation_baseplot} initalizes the annotation plot.
#'
#' This is wrapper that includes loding the annotations, and initalizing the ggplot.
#' Individual data tracks should be realized using facetting over the "window" column.
#'
#' @param searchLG string skalar (mandatory), should be on of "LG01" - "LG24"
#' @param xrange integer vector (mandatory), data range to be loaded(start bp - end bp).
#'   Positions are defined with respect to the LG, NOT the overall genomic position.
#' @param genes_of_interest string vector (optional), tags specific gene names for highlighting when plotting.
#'   (needs to exactly match the gene name of the original gff file)
#' @param genes_of_sec_interest string vector (optional), tags specific gene names for secondary highlighting
#'  when plotting. (needs to exactly match the gene name of the original gff file)
#' @param anno_rown integer skalar (optional), defines the number of gene rows to avoid overlapping.
#' @param width float skalar (0-1, optional), defines the height of the exon boxes.
#'
#' @seealso \code{\link{hypo_annotation_get}}, \code{\link{theme_hypo_anno_extra}}
#'
#' @examples
#'
#' hypo_annotation_baseplot(searchLG = "LG03", xrange = c(0,100000),
#'   genes_of_interest = 'SUN1_4') +
#'   geom_line(data = dataframe,
#'             aes(x = x, y = y, color = window)) +
#'   facet_grid(window~., scales = 'free_y',
#'              switch = 'y',labeller = label_parsed,
#'              as.table = TRUE)
#'
#' @export
hypo_annotation_baseplot <- function(..., searchLG, xrange, genes_of_interest = c(),
                                 genes_of_sec_interest = c(), anno_rown = 3,
                                 width = .1){

  df_list <- hypo_annotation_get(searchLG = searchLG,xrange = xrange,
                                 genes_of_interest = genes_of_interest,
                                 genes_of_sec_interest = genes_of_sec_interest,
                                 anno_rown = anno_rown)

  ggplot2::ggplot()+
    ggplot2::geom_rect(data = df_list[[2]],
            aes(xmin = ps, xmax = pe, ymin = yl-(width/2), ymax = yl+(width/2),
                fill = clr, col = clr, group = Parent),
            lwd = .1)+
    ggplot2::geom_segment(data = (df_list[[1]] %>% filter(strand %in% c('+', '-'))),
               aes(x = ps, xend = pe, y = yl, yend = yl, group = Parent),
               lwd=.2,arrow = arrow(length = unit(2,"pt"),type = 'closed'),
               color = 'black')+
    ggplot2::geom_segment(data = (df_list[[1]] %>% filter(!strand %in% c('+', '-'))),
               aes(x = ps, xend = pe, y = yl, yend = yl, group = Parent),
               lwd = .2, color = 'black')+
    ggplot2::geom_text(data = df_list[[1]],
            aes(x = labelx, label = gsub('000000', '...', label),y = yl-.5))
  }

