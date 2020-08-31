rm(list=ls())
ls()
library(tidyverse)
library(gridExtra)
#devtools::install_github("thomasp85/patchwork")
library(patchwork)
library(lattice)
library(grid)
library(gridExtra)

setwd("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/data/coverage_analysis")

# import number of reads mapped to scaffold

PV04.reads.mapped <- read_delim("PV_18-04.reads.mapped.count", 
                                    "\t", escape_double = FALSE, col_names = FALSE, 
                                    trim_ws = TRUE)
PV13.reads.mapped <- read_delim("PV_18-13.reads.mapped.count", 
                                "\t", escape_double = FALSE, col_names = FALSE, 
                                trim_ws = TRUE)
PV21.reads.mapped <- read_delim("PV_18-21.reads.mapped.count", 
                                "\t", escape_double = FALSE, col_names = FALSE, 
                                trim_ws = TRUE)
PV23.reads.mapped <- read_delim("PV_18-23.reads.mapped.count", 
                                "\t", escape_double = FALSE, col_names = FALSE, 
                                trim_ws = TRUE)

colnames(PV04.reads.mapped) <- c("seq","length","mapped","unmapped")
colnames(PV13.reads.mapped) <- c("seq","length","mapped","unmapped")
colnames(PV21.reads.mapped) <- c("seq","length","mapped","unmapped")
colnames(PV23.reads.mapped) <- c("seq","length","mapped","unmapped")

# explore differences in coverage. Let's always do +1 to avoid 0s. We can plot this in all pairs of samples. 

cov.13v21 <- log2((PV13.reads.mapped$mapped+1)/(PV21.reads.mapped$mapped+1))
cov.13v23 <- log2((PV13.reads.mapped$mapped+1)/(PV23.reads.mapped$mapped+1))
cov.04v21 <- log2((PV04.reads.mapped$mapped+1)/(PV21.reads.mapped$mapped+1))
cov.04v23 <- log2((PV04.reads.mapped$mapped+1)/(PV23.reads.mapped$mapped+1))
cov.04v13 <- log2((PV04.reads.mapped$mapped+1)/(PV13.reads.mapped$mapped+1))
cov.21v23 <- log2((PV21.reads.mapped$mapped+1)/(PV23.reads.mapped$mapped+1))

cov.diff <- data.frame(cov.13v21,cov.13v23,cov.04v21,cov.04v23,cov.04v13,cov.21v23)

p1 <- ggplot(cov.diff, aes(cov.13v21)) + geom_histogram(bins=150)
p2 <- ggplot(cov.diff, aes(cov.13v23)) + geom_histogram(bins=150)
p3 <- ggplot(cov.diff, aes(cov.04v21)) + geom_histogram(bins=150)
p4 <- ggplot(cov.diff, aes(cov.04v23)) + geom_histogram(bins=150)
p5 <- ggplot(cov.diff, aes(cov.04v13)) + geom_histogram(bins=150)
p6 <- ggplot(cov.diff, aes(cov.21v23)) + geom_histogram(bins=150)
hist <- p1 + p2 + p3 + p4 + p5 + p6
