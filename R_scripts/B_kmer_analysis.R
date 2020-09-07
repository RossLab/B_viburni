rm(list=ls())
ls()
library(tidyverse)
library(gridExtra)
#devtools::install_github("thomasp85/patchwork")
library(patchwork)
library(lattice)
library(grid)
library(gridExtra)
library(reshape2)
library(plyr)

setwd("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/data/coverage_analysis")

PV04_kmer_counts <- read_delim("PV_18-04_kmer_counts.2.hist","\t", col_names = FALSE, escape_double = FALSE, trim_ws = TRUE)
PV13_kmer_counts <- read_delim("PV_18-13_kmer_counts.2.hist","\t", col_names = FALSE, escape_double = FALSE, trim_ws = TRUE)
PV21_kmer_counts <- read_delim("PV_18-21_kmer_counts.2.hist","\t", col_names = FALSE, escape_double = FALSE, trim_ws = TRUE)
PV23_kmer_counts <- read_delim("PV_18-23_kmer_counts.2.hist","\t", col_names = FALSE, escape_double = FALSE, trim_ws = TRUE)

PV04_kmer_counts_filt <- PV04_kmer_counts[5:200,]
PV13_kmer_counts_filt <- PV13_kmer_counts[5:200,]
PV21_kmer_counts_filt <- PV21_kmer_counts[5:200,]
PV23_kmer_counts_filt <- PV23_kmer_counts[5:200,]

p1 <- qplot(X1,X2,data=PV04_kmer_counts_filt, geom="line")
p2 <- qplot(X1,X2,data=PV13_kmer_counts_filt, geom="line")
p3 <- qplot(X1,X2,data=PV21_kmer_counts_filt, geom="line")
p4 <- qplot(X1,X2,data=PV23_kmer_counts_filt, geom="line")

kmer_hist <- p1 + p2 + p3 + p4

