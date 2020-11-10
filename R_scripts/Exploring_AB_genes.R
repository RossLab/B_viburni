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

##### Get things ready

# setwd and import data
setwd("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/data")

# master anno
freeze.v0.genes.anno <- read_delim("output/freeze.v0.genes.anno.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE)

# B chromosome assignment
scaffolds.final.assignment <- read_delim("output/scaffolds.final.assignment.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE)

# expression data from RSEM
X04F_1 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/04F_1.genes.results")
X04F_2 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/04F_2.genes.results")
X04F_3 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/04F_3.genes.results")
X04M_1 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/04M_1.genes.results")
X04M_2 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/04M_2.genes.results")
X04M_3 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/04M_3.genes.results")
X13F_1 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/13F_1.genes.results")
X13F_2 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/13F_2.genes.results")
X13F_3 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/13F_3.genes.results")
X13M_1 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/13M_1.genes.results")
X13M_2 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/13M_2.genes.results")
X13M_3 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/13M_3.genes.results")
X13M_4 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/13M_4.genes.results")
X15F_1 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/15F_1.genes.results")
X15F_2 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/15F_2.genes.results")
X15F_3 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/15F_3.genes.results")
X15M_1 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/15M_1.genes.results")
X15M_2 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/15M_2.genes.results")
X15M_3 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/15M_3.genes.results")
X21F_1 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/21F_1.genes.results")
X21F_2 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/21F_2.genes.results")
X21F_3 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/21F_3.genes.results")
X21M_1 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/21M_1.genes.results")
X21M_2 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/21M_2.genes.results")
X21M_3 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/21M_3.genes.results")
X21M_4 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/21M_4.genes.results")

colnames(X04F_1)[6] <- "X04F_1"
colnames(X04F_2)[6] <- "X04F_2"
colnames(X04F_3)[6] <- "X04F_3"
colnames(X04M_1)[6] <- "X04M_1"
colnames(X04M_2)[6] <- "X04M_2"
colnames(X04M_3)[6] <- "X04M_3"
colnames(X13F_1)[6] <- "X13F_1"
colnames(X13F_2)[6] <- "X13F_2"
colnames(X13F_3)[6] <- "X13F_3"
colnames(X13M_1)[6] <- "X13M_1"
colnames(X13M_2)[6] <- "X13M_2"
colnames(X13M_3)[6] <- "X13M_3"
colnames(X13M_4)[6] <- "X13M_4"
colnames(X15F_1)[6] <- "X15F_1"
colnames(X15F_2)[6] <- "X15F_2"
colnames(X15F_3)[6] <- "X15F_3"
colnames(X15M_1)[6] <- "X15M_1"
colnames(X15M_2)[6] <- "X15M_2"
colnames(X15M_3)[6] <- "X15M_3"
colnames(X21F_1)[6] <- "X21F_1"
colnames(X21F_2)[6] <- "X21F_2"
colnames(X21F_3)[6] <- "X21F_3"
colnames(X21M_1)[6] <- "X21M_1"
colnames(X21M_2)[6] <- "X21M_2"
colnames(X21M_3)[6] <- "X21M_3"
colnames(X21M_4)[6] <- "X21M_4"

rsem.tpm <- X04F_1[c(1,6)]
rsem.tpm <- merge(rsem.tpm, X04F_2[c(1,6)], by="gene_id")
rsem.tpm <- merge(rsem.tpm, X04F_3[c(1,6)], by="gene_id")
rsem.tpm <- merge(rsem.tpm, X04M_1[c(1,6)], by="gene_id")
rsem.tpm <- merge(rsem.tpm, X04M_2[c(1,6)], by="gene_id")
rsem.tpm <- merge(rsem.tpm, X04M_3[c(1,6)], by="gene_id")
rsem.tpm <- merge(rsem.tpm, X13F_1[c(1,6)], by="gene_id")
rsem.tpm <- merge(rsem.tpm, X13F_2[c(1,6)], by="gene_id")
rsem.tpm <- merge(rsem.tpm, X13F_3[c(1,6)], by="gene_id")
rsem.tpm <- merge(rsem.tpm, X13M_1[c(1,6)], by="gene_id")
rsem.tpm <- merge(rsem.tpm, X13M_2[c(1,6)], by="gene_id")
rsem.tpm <- merge(rsem.tpm, X13M_3[c(1,6)], by="gene_id")
rsem.tpm <- merge(rsem.tpm, X13M_4[c(1,6)], by="gene_id")
rsem.tpm <- merge(rsem.tpm, X15F_1[c(1,6)], by="gene_id")
rsem.tpm <- merge(rsem.tpm, X15F_2[c(1,6)], by="gene_id")
rsem.tpm <- merge(rsem.tpm, X15F_3[c(1,6)], by="gene_id")
rsem.tpm <- merge(rsem.tpm, X15M_1[c(1,6)], by="gene_id")
rsem.tpm <- merge(rsem.tpm, X15M_2[c(1,6)], by="gene_id")
rsem.tpm <- merge(rsem.tpm, X15M_3[c(1,6)], by="gene_id")
rsem.tpm <- merge(rsem.tpm, X21F_1[c(1,6)], by="gene_id")
rsem.tpm <- merge(rsem.tpm, X21F_2[c(1,6)], by="gene_id")
rsem.tpm <- merge(rsem.tpm, X21F_3[c(1,6)], by="gene_id")
rsem.tpm <- merge(rsem.tpm, X21M_1[c(1,6)], by="gene_id")
rsem.tpm <- merge(rsem.tpm, X21M_2[c(1,6)], by="gene_id")
rsem.tpm <- merge(rsem.tpm, X21M_3[c(1,6)], by="gene_id")
rsem.tpm <- merge(rsem.tpm, X21M_4[c(1,6)], by="gene_id")
head(X21M_2)

# differentially expressed genes -- from Andrés' rerun of Isabelle's script (received 03.11.20, rerun completed 09.11.20)

dt_df <- read_delim("output/diff_expr/dt_df.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE) # complete list of contrasts
de.over.B.males.genes <- read_delim("output/diff_expr/over.Bmales.vs.all.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE)
de.under.B.males.genes <- as.data.frame(c("g4126","g8486","g14550"))

de.B.males.vs.nonB.males <- read_delim("output/diff_expr/B.males.vs.nonB.males.de.treat.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE)
de.B.males.vs.B.females <- read_delim("output/diff_expr/B.males.vs.B.females.de.treat.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE)
de.B.males.vs.nonB.females <- read_delim("output/diff_expr/B.males.vs.nonB.females.de.treat.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE)
de.B.females.vs.nonB.females <- read_delim("output/diff_expr/B.females.vs.nonB.females.de.treat.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE)
de.nonB.males.vs.nonB.females <- read_delim("output/diff_expr/nonB.males.vs.nonB.females.de.treat.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE)

###

##### Genes located on B candidate scaffolds

genes.AB <- merge(freeze.v0.genes.anno,scaffolds.final.assignment[c("seq","length","b.status.final","cov.04v13","b.status","b.status.asn","b.status.kmer")],by="seq")
genes.A <-  genes.AB[genes.AB$b.status.final == "A",]
genes.B1 <- genes.AB[genes.AB$b.status.final == "B1",]
genes.B2 <- genes.AB[genes.AB$b.status.final == "B2",]
genes.B3 <- genes.AB[genes.AB$b.status.final == "B3",]
genes.B4 <- genes.AB[genes.AB$b.status.final == "B4",]

genes.in.Bs.anno <- rbind(genes.B1[genes.B1$anno == "Y",], genes.B2[genes.B2$anno == "Y",], genes.B3[genes.B3$anno == "Y",], genes.B4[genes.B4$anno == "Y",])[c(2,9,1,3,4,5,6)]
#write.table(genes.in.Bs.anno, file = "output/genes.in.Bs.anno.csv",row.names = F,sep = ",")

###

##### Examine genes that are overexpressed in B males compared to all other samples (non-B males and both B/non-B females)

de.over.B.males.genes <- read_delim("output/diff_expr/over.Bmales.vs.all.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE)
de.over.B.males.genes.anno <- left_join(de.over.B.males.genes,genes.AB,by="gene")
table(de.over.B.males.genes.anno$b.status.final)
#write.csv(de.over.B.males.genes.anno,"output/diff_expr/over.Bmales.vs.all.anno.csv")

de.under.B.males.genes <- as.data.frame(c("g4126","g8486","g14550"))
colnames(de.under.B.males.genes)[1] <- "gene"
de.under.B.males.genes.anno <- left_join(de.under.B.males.genes,genes.AB,by="gene")
#write.csv(de.under.B.males.genes.anno,"output/diff_expr/under.Bmales.vs.all.anno.csv")

###

##### Examine differentially expressed genes between groups

de.B.males.vs.nonB.males <- de.B.males.vs.nonB.males[c(-1)]
de.B.males.vs.B.females <- de.B.males.vs.B.females[c(-1)]
de.B.males.vs.nonB.females <- de.B.males.vs.nonB.females[c(-1)]
de.B.females.vs.nonB.females <- de.B.females.vs.nonB.females[c(-1)]
de.nonB.males.vs.nonB.females <- de.nonB.males.vs.nonB.females[c(-1)]

colnames(de.B.males.vs.nonB.males)[1] <- "gene"
colnames(de.B.males.vs.B.females)[1] <- "gene"
colnames(de.B.males.vs.nonB.females)[1] <- "gene"
colnames(de.B.females.vs.nonB.females)[1] <- "gene"
colnames(de.nonB.males.vs.nonB.females)[1] <- "gene"

de.B.males.vs.nonB.males.anno <- left_join(de.B.males.vs.nonB.males,genes.AB,by="gene")
de.B.males.vs.B.females.anno <- left_join(de.B.males.vs.B.females,genes.AB,by="gene")
de.B.males.vs.nonB.female.anno <- left_join(de.B.males.vs.nonB.females,genes.AB,by="gene")
de.B.females.vs.nonB.females.anno <- left_join(de.B.females.vs.nonB.females,genes.AB,by="gene")
de.nonB.males.vs.nonB.females.anno <- left_join(de.nonB.males.vs.nonB.females,genes.AB,by="gene")

## B males v non B males

de.over.B.males.vs.nonB.males.anno <- de.B.males.vs.nonB.males.anno[de.B.males.vs.nonB.males.anno$logFC > 1,]
count(de.over.B.males.vs.nonB.males.anno$b.status.final)
de.over.B.males.vs.nonB.males.anno.yes <- de.over.B.males.vs.nonB.males.anno[de.over.B.males.vs.nonB.males.anno$anno == "Y",]
count(de.over.B.males.vs.nonB.males.anno.yes$b.status.final)
#write.csv(de.over.B.males.vs.nonB.males.anno,"output/diff_expr/de.over.B.males.vs.nonB.males.anno.csv")

de.under.B.males.vs.nonB.males.anno <- de.B.males.vs.nonB.males.anno[de.B.males.vs.nonB.males.anno$logFC < 1,]
count(de.under.B.males.vs.nonB.males.anno$b.status.final)
de.under.B.males.vs.nonB.males.anno.yes <- de.under.B.males.vs.nonB.males.anno[de.under.B.males.vs.nonB.males.anno$anno == "Y",]
count(de.under.B.males.vs.nonB.males.anno.yes$b.status.final)
#write.csv(de.under.B.males.vs.nonB.males.anno,"output/diff_expr/de.under.B.males.vs.nonB.males.anno.csv")

## B females vs non B females

de.over.B.females.vs.nonB.females.anno <- de.B.females.vs.nonB.females.anno[de.B.females.vs.nonB.females.anno$logFC > 1,]
count(de.over.B.females.vs.nonB.females.anno$b.status.final)
de.over.B.females.vs.nonB.females.anno.yes <- de.over.B.females.vs.nonB.females.anno[de.over.B.females.vs.nonB.females.anno$anno == "Y",]
count(de.over.B.females.vs.nonB.females.anno.yes$b.status.final)
#write.csv(de.over.B.females.vs.nonB.females.anno,"output/diff_expr/de.over.B.females.vs.nonB.females.anno.csv")

de.under.B.females.vs.nonB.females.anno <- de.B.females.vs.nonB.females.anno[de.B.females.vs.nonB.females.anno$logFC < 1,]
count(de.under.B.females.vs.nonB.females.anno$b.status.final)
de.under.B.females.vs.nonB.females.anno.yes <- de.under.B.females.vs.nonB.females.anno[de.under.B.females.vs.nonB.females.anno$anno == "Y",]
count(de.under.B.females.vs.nonB.females.anno.yes$b.status.final)
#write.csv(de.under.B.females.vs.nonB.females.anno,"output/diff_expr/de.under.B.females.vs.nonB.females.anno.csv")

## B females vs non B females

de.B.males.vs.nonB.female.anno <- left_join(de.B.males.vs.nonB.females,genes.AB,by="gene")
de.B.females.vs.nonB.females.anno <- left_join(de.B.females.vs.nonB.females,genes.AB,by="gene")
de.nonB.males.vs.nonB.females.anno <- left_join(de.nonB.males.vs.nonB.females,genes.AB,by="gene")

## remaining comparisons

de.B.males.vs.B.females.anno <- left_join(de.B.males.vs.B.females,genes.AB,by="gene")
de.B.males.vs.nonB.females.anno <- left_join(de.B.males.vs.nonB.females,genes.AB,by="gene")
de.nonB.males.vs.nonB.females.anno <- left_join(de.nonB.males.vs.nonB.females,genes.AB,by="gene")

#write.csv(de.B.males.vs.B.females.anno,"output/diff_expr/de.B.males.vs.B.females.anno.csv")
#write.csv(de.B.males.vs.nonB.females.anno,"output/diff_expr/de.B.males.vs.nonB.females.anno.csv")
#write.csv(de.nonB.males.vs.nonB.females.anno,"output/diff_expr/de.nonB.males.vs.nonB.females.anno.csv")

# extract putative B genes that are differentially expressed between B-carrying males and females
de.B.males.vs.B.females.anno.b <- de.B.males.vs.B.females.anno[de.B.males.vs.B.females.anno$b.status.final != "A",]
#write.csv(de.B.males.vs.B.females.anno.b,"output/diff_expr/de.B.males.vs.B.females.anno.b.csv")

###

##### Re-evaluating genes in B scaffolds
genes.B1.anno <- genes.AB[genes.AB$b.status.final == "B1",]
genes.B1.anno.dt1 <- left_join(genes.B1.anno,dt_df[c(-1)],by="gene")
genes.B1.anno.dt1.expr <- left_join(genes.B1.anno.dt1, rsem.counts, by = "gene")

