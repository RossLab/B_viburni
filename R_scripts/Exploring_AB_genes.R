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

# genes on the B chromosomes

genes.AB <- merge(freeze.v0.genes.anno,scaffolds.final.assignment[c("seq","length","b.status.final","cov.04v13","b.status","b.status.asn","b.status.kmer")],by="seq")

genes.A <-  genes.AB[genes.AB$b.status.final == "A",]
genes.B1 <- genes.AB[genes.AB$b.status.final == "B1",]
genes.B2 <- genes.AB[genes.AB$b.status.final == "B2",]
genes.B3 <- genes.AB[genes.AB$b.status.final == "B3",]
genes.B4 <- genes.AB[genes.AB$b.status.final == "B4",]

genes.in.Bs.anno <- rbind(genes.B1[genes.B1$anno == "Y",], genes.B2[genes.B2$anno == "Y",], genes.B3[genes.B3$anno == "Y",], genes.B4[genes.B4$anno == "Y",])[c(2,9,1,3,4,5,6)]
#write.table(genes.in.Bs.anno, file = "output/genes.in.Bs.anno.csv",row.names = F,sep = ",")

### examine genes that are overexpressed in B males compared to non-B males, and both B/non-B females
de.over.B.males.genes <- read_delim("output/diff_expr/over.Bmales.vs.all.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE)
de.over.B.males.genes.anno <- left_join(de.over.B.males.genes,genes.AB,by="gene")
table(de.over.B.males.genes.anno$b.status.final)
#write.csv(de.over.B.males.genes.anno,"output/diff_expr/over.Bmales.vs.all.anno.csv")

de.under.B.males.genes <- as.data.frame(c("g4126","g8486","g14550"))
colnames(de.under.B.males.genes)[1] <- "gene"
de.under.B.males.genes.anno <- left_join(de.under.B.males.genes,genes.AB,by="gene")
#write.csv(de.under.B.males.genes.anno,"output/diff_expr/under.Bmales.vs.all.anno.csv")
