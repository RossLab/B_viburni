rm(list=ls())
ls()
library(plyr)
library(tidyverse)
library(gridExtra)
library(patchwork)
library(lattice)
library(grid)
library(gridExtra)
library(reshape2)

##### setwd and import data

setwd("E:/agdel/Documents/projects_rosslab/B_viburni")

# ------------------------------------------------------
# Make a compatible GO annotation file
Pcitri_genes_with_GO <- read_table2("output/pviburni.gene.GO", 
                                    col_names = FALSE)
colnames(Pcitri_genes_with_GO) <- c("gene","go")
new_annotations <- separate_rows(Pcitri_genes_with_GO, go, sep =';')
#write.table(new_annotations, file="P_citri_GO_terms.txt", sep="\t", quote = F,
#            col.names = T, row.names = F)

# ------------------------------------------------------
# Make the expression gene lists
logFC_DEgenes <- read_delim("output/sex_diff_expr/fdr_log/FvsM_results_anno.csv", 
                            ",", escape_double = FALSE, col_names = T, 
                            trim_ws = TRUE)
