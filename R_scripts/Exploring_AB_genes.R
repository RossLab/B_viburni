m(list=ls())
ls()
library(plyr)
library(tidyverse)
library(gridExtra)
library(patchwork)
library(lattice)
library(grid)
library(gridExtra)
library(reshape2)

# setwd and import data

setwd("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/data")

# master anno
freeze.v0.genes.anno <- read_delim("output/freeze.v0.genes.anno.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE)

# B chromosome assignment
scaffolds.final.assignment <- read_delim("output/scaffolds.final.assignment.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE)

# differentially expressed genes -- from Isabelle (received 03.11.20)
fb.vs.fnob <- read_delim("output/rsem_gene_femaleB.femalenoB.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE)
mb.vs.fb <- read_delim("output/rsem_gene_maleB.femaleB.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE)
mb.vs.fnob <- read_delim("output/rsem_gene_maleB.femalenoB.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE)
mb.vs.mnob <- read_delim("output/rsem_gene_maleB.noB_annot.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE)
mnob.vs.fnb <- read_delim("output/rsem_gene_malenoB.femalenoB.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE)