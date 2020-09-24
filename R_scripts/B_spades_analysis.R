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

# import files
pviburni.freeze <- read_delim("p.viburni.freeze.v0.softmasked.fa.fai", 
                              "\t", escape_double = FALSE, col_names = FALSE, 
                              trim_ws = TRUE)
B_strict.candidates <- read_csv("B.strict.candidates", 
                                col_names = FALSE)
B_strict.candidates$B.candidate <- "Y"

colnames(pviburni.freeze)[1] <- "scaffold"
colnames(pviburni.freeze)[2] <- "len"

B_strict.candidates$B.candidate <- "Y"
colnames(B_strict.candidates)[1] <- "scaffold"

spades.nucmer.04 <- read_csv("spades/04.spades.v.freeze.v0.scaffolds.list", col_names = FALSE)
spades.nucmer.13 <- read_csv("spades/13.spades.v.freeze.v0.scaffolds.list", col_names = FALSE)
spades.nucmer.21 <- read_csv("spades/21.spades.v.freeze.v0.scaffolds.list", col_names = FALSE)
spades.nucmer.23 <- read_csv("spades/23.spades.v.freeze.v0.scaffolds.list", col_names = FALSE)

# create

pviburni.scaffold.table <- left_join(pviburni.freeze[c(1,2)],B_strict.candidates,by="scaffold")

colnames(spades.nucmer.04)[1] <- "scaffold"
colnames(spades.nucmer.13)[1] <- "scaffold"
colnames(spades.nucmer.21)[1] <- "scaffold"
colnames(spades.nucmer.23)[1] <- "scaffold"
spades.nucmer.04$PV.04 <- "Y"
spades.nucmer.13$PV.13 <- "Y"
spades.nucmer.21$PV.21 <- "Y"
spades.nucmer.23$PV.23 <- "Y"

pviburni.scaffold.table <- left_join(pviburni.scaffold.table,spades.nucmer.04,by="scaffold")
pviburni.scaffold.table <- left_join(pviburni.scaffold.table,spades.nucmer.13,by="scaffold")
pviburni.scaffold.table <- left_join(pviburni.scaffold.table,spades.nucmer.21,by="scaffold")
pviburni.scaffold.table <- left_join(pviburni.scaffold.table,spades.nucmer.23,by="scaffold")

# how many scaffolds are in the B+ lines and not in the B- lines?

pviburni.scaffold.in.Bplus <- pviburni.scaffold.table[!is.na(pviburni.scaffold.table$PV.13) & !is.na(pviburni.scaffold.table$PV.04)
                                                      & is.na(pviburni.scaffold.table$PV.21) & is.na(pviburni.scaffold.table$PV.23)
                                                      & !is.na(pviburni.scaffold.table$B.candidate),]

count(pviburni.scaffold.in.Bplus$B.candidate)
sum(pviburni.scaffold.in.Bplus$len)
