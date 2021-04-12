## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  dev = "ragg_png"
)

library("ragg")
library("hypogen")
library("svglite")

## ----load_libs, eval = TRUE, message = FALSE, warning = FALSE-----------------
library(tidyverse)
library(ggtern)

## ----data_import, message = FALSE, warning = FALSE----------------------------
file_genotypes <- system.file("extdata", "genotype_table.012.trans.txt.gz", package = "hypogen")

genotype_freqs <- hypo_import_genotype_freq(file_genotypes)

## ----plot, eval = TRUE, fig.width=10,fig.height=6, out.width="100%", message = FALSE, warning = FALSE----
ggplot() +
  coord_tern()+
  geom_point(data = genotype_freqs,
             aes(x = AA, y = Aa, z = aa ),
             col = hypo_clr_LGs[18]) +
  geom_path(data = hypo_hwe(n = 100),
            aes(x = AA, y = Aa, z = aa ),
            color ='black',alpha = .6)

## ---- echo = FALSE, fig.asp = 1, out.height = "150pt", out.width = "150pt", fig.align = "center"----
grImport2::grid.picture(grImport2::readPicture("logo.c.svg"))

