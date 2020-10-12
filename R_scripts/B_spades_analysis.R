rm(list=ls())
ls()
library(tidyverse)
library(gridExtra)
library(patchwork)
library(lattice)
library(grid)
library(gridExtra)
library(reshape2)

A.put.0 <- read.delim("~/Desktop/Aput_kmers_vs_freeze.vo.number.mapped", header=FALSE)
B.put.0 <- read.delim("~/Desktop/Bput_kmers_vs_freeze.vo.number.mapped", header=FALSE)

A.put <- A.put.0[c(1,3)]
B.put <- B.put.0[c(1,3)]

colnames(A.put)[1] <- "scaffold"
colnames(A.put)[2] <- "A"

colnames(B.put)[1] <- "scaffold"
colnames(B.put)[2] <- "B"

AB <- left_join(A.put,B.put,by=c("scaffold"))
View(AB)
