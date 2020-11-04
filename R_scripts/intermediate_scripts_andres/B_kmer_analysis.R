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

setwd("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/data/coverage_analysis/kmer")

PV04_kmer_counts <- read_delim("PV_18-04_kmer_counts.2.hist","\t", col_names = FALSE, escape_double = FALSE, trim_ws = TRUE)
PV13_kmer_counts <- read_delim("PV_18-13_kmer_counts.2.hist","\t", col_names = FALSE, escape_double = FALSE, trim_ws = TRUE)
PV21_kmer_counts <- read_delim("PV_18-21_kmer_counts.2.hist","\t", col_names = FALSE, escape_double = FALSE, trim_ws = TRUE)
PV23_kmer_counts <- read_delim("PV_18-23_kmer_counts.2.hist","\t", col_names = FALSE, escape_double = FALSE, trim_ws = TRUE)

PV04_kmer_counts_filt <- PV04_kmer_counts[5:200,]
PV13_kmer_counts_filt <- PV13_kmer_counts[5:200,]
PV21_kmer_counts_filt <- PV21_kmer_counts[5:200,]
PV23_kmer_counts_filt <- PV23_kmer_counts[5:200,]

p1 <- qplot(X1,X2,data=PV04_kmer_counts_filt, geom="line") + scale_y_continuous(limits = c(0,1.2e+07), labels = scales::scientific) + labs(title="PV04",x="log10(kmer cov)", y="Count") + theme_bw()
p2 <- qplot(X1,X2,data=PV13_kmer_counts_filt, geom="line") + scale_y_continuous(limits = c(0,1.2e+07), labels = scales::scientific) + labs(title="PV13",x="log10(kmer cov)", y="Count") + theme_bw()
p3 <- qplot(X1,X2,data=PV21_kmer_counts_filt, geom="line") + scale_y_continuous(limits = c(0,1.2e+07), labels = scales::scientific) + labs(title="PV21",x="log10(kmer cov)", y="Count") + theme_bw()
p4 <- qplot(X1,X2,data=PV23_kmer_counts_filt, geom="line") + scale_y_continuous(limits = c(0,1.2e+07), labels = scales::scientific) + labs(title="PV23",x="log10(kmer cov)", y="Count") + theme_bw()
kmer_hist <- p1 + p2 + p3 + p4

# import histograms with higher upper limit

PV04.kmer.counts.max <- read_delim("PV_18-04_kmer_counts_round2.max.hist","\t", col_names = FALSE, escape_double = FALSE, trim_ws = TRUE)
PV13.kmer.counts.max <- read_delim("PV_18-13_kmer_counts_round2.max.hist","\t", col_names = FALSE, escape_double = FALSE, trim_ws = TRUE)
PV21.kmer.counts.max <- read_delim("PV_18-21_kmer_counts_round2.max.hist","\t", col_names = FALSE, escape_double = FALSE, trim_ws = TRUE)
PV23.kmer.counts.max <- read_delim("PV_18-23_kmer_counts_round2.max.hist","\t", col_names = FALSE, escape_double = FALSE, trim_ws = TRUE)

PV04.kmer.counts.max.no.0 <- PV04.kmer.counts.max[PV04.kmer.counts.max$X2 > 0,]
PV13.kmer.counts.max.no.0 <- PV13.kmer.counts.max[PV13.kmer.counts.max$X2 > 0,]
PV21.kmer.counts.max.no.0 <- PV21.kmer.counts.max[PV21.kmer.counts.max$X2 > 0,]
PV23.kmer.counts.max.no.0 <- PV23.kmer.counts.max[PV23.kmer.counts.max$X2 > 0,]
tail(PV04.kmer.counts.max.no.0)
tail(PV13.kmer.counts.max.no.0)
tail(PV21.kmer.counts.max.no.0)
tail(PV23.kmer.counts.max.no.0)

p1 <- qplot(log10(X1),X2,data=PV04.kmer.counts.max.no.0, geom="line") + scale_y_continuous(limits = c(0,1.2e+07), labels = scales::scientific) + labs(title="PV04",x="log10(kmer cov)", y="count") + theme_bw()
p2 <- qplot(log10(X1),X2,data=PV13.kmer.counts.max.no.0, geom="line") + scale_y_continuous(limits = c(0,1.2e+07), labels = scales::scientific) + labs(title="PV13",x="log10(kmer cov)", y="count") + theme_bw()
p3 <- qplot(log10(X1),X2,data=PV21.kmer.counts.max.no.0, geom="line") + scale_y_continuous(limits = c(0,1.2e+07), labels = scales::scientific) + labs(title="PV21",x="log10(kmer cov)", y="count") + theme_bw()
p4 <- qplot(log10(X1),X2,data=PV23.kmer.counts.max.no.0, geom="line") + scale_y_continuous(limits = c(0,1.2e+07), labels = scales::scientific) + labs(title="PV23",x="log10(kmer cov)", y="count") + theme_bw()
q1 <- qplot(log10(X1),log10(X2),data=PV04.kmer.counts.max.no.0, geom="line") + labs(title="PV04",x="log10(kmer cov)", y="log10(count)") + theme_bw()
q2 <- qplot(log10(X1),log10(X2),data=PV13.kmer.counts.max.no.0, geom="line") + labs(title="PV13",x="log10(kmer cov)", y="log10(count)") + theme_bw()
q3 <- qplot(log10(X1),log10(X2),data=PV21.kmer.counts.max.no.0, geom="line") + labs(title="PV21",x="log10(kmer cov)", y="log10(count)") + theme_bw()
q4 <- qplot(log10(X1),log10(X2),data=PV23.kmer.counts.max.no.0, geom="line") + labs(title="PV23",x="log10(kmer cov)", y="log10(count)") + theme_bw()
 
kmer_hist_100000 <- p1 + p2 + p3 + p4 + q1 + q2 + q3 + q4 + plot_layout(ncol = 4)

# edit for genomescope

PV04.kmer.counts.1e6 <- PV04.kmer.counts.max
PV13.kmer.counts.1e6 <- PV13.kmer.counts.max
PV21.kmer.counts.1e6 <- PV21.kmer.counts.max
PV23.kmer.counts.1e6 <- PV23.kmer.counts.max







