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
table(scaffolds.final.assignment$b.status.final)
# expression data from RSEM

#this is how I obtained TPM values from RSEM
#X04F_1 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/04F_1.genes.results")
#X04F_2 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/04F_2.genes.results")
#X04F_3 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/04F_3.genes.results")
#X04M_1 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/04M_1.genes.results")
#X04M_2 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/04M_2.genes.results")
#X04M_3 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/04M_3.genes.results")
#X13F_1 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/13F_1.genes.results")
#X13F_2 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/13F_2.genes.results")
#X13F_3 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/13F_3.genes.results")
#X13M_1 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/13M_1.genes.results")
#X13M_2 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/13M_2.genes.results")
#X13M_3 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/13M_3.genes.results")
#X13M_4 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/13M_4.genes.results")
#X15F_1 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/15F_1.genes.results")
#X15F_2 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/15F_2.genes.results")
#X15F_3 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/15F_3.genes.results")
#X15M_1 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/15M_1.genes.results")
#X15M_2 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/15M_2.genes.results")
#X15M_3 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/15M_3.genes.results")
#X21F_1 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/21F_1.genes.results")
#X21F_2 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/21F_2.genes.results")
#X21F_3 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/21F_3.genes.results")
#X21M_1 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/21M_1.genes.results")
#X21M_2 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/21M_2.genes.results")
#X21M_3 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/21M_3.genes.results")
#X21M_4 <- read.delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/R_scripts/input/21M_4.genes.results")
#
#colnames(X04F_1)[6] <- "X04F_1"
#colnames(X04F_2)[6] <- "X04F_2"
#colnames(X04F_3)[6] <- "X04F_3"
#colnames(X04M_1)[6] <- "X04M_1"
#colnames(X04M_2)[6] <- "X04M_2"
#colnames(X04M_3)[6] <- "X04M_3"
#colnames(X13F_1)[6] <- "X13F_1"
#colnames(X13F_2)[6] <- "X13F_2"
#colnames(X13F_3)[6] <- "X13F_3"
#colnames(X13M_1)[6] <- "X13M_1"
#colnames(X13M_2)[6] <- "X13M_2"
#colnames(X13M_3)[6] <- "X13M_3"
#colnames(X13M_4)[6] <- "X13M_4"
#colnames(X15F_1)[6] <- "X15F_1"
#colnames(X15F_2)[6] <- "X15F_2"
#colnames(X15F_3)[6] <- "X15F_3"
#colnames(X15M_1)[6] <- "X15M_1"
#colnames(X15M_2)[6] <- "X15M_2"
#colnames(X15M_3)[6] <- "X15M_3"
#colnames(X21F_1)[6] <- "X21F_1"
#colnames(X21F_2)[6] <- "X21F_2"
#colnames(X21F_3)[6] <- "X21F_3"
#colnames(X21M_1)[6] <- "X21M_1"
#colnames(X21M_2)[6] <- "X21M_2"
#colnames(X21M_3)[6] <- "X21M_3"
#colnames(X21M_4)[6] <- "X21M_4"
#
#rsem.tpm <- X04F_1[c(1,6)]
#rsem.tpm <- merge(rsem.tpm, X04F_2[c(1,6)], by="gene_id")
#rsem.tpm <- merge(rsem.tpm, X04F_3[c(1,6)], by="gene_id")
#rsem.tpm <- merge(rsem.tpm, X04M_1[c(1,6)], by="gene_id")
#rsem.tpm <- merge(rsem.tpm, X04M_2[c(1,6)], by="gene_id")
#rsem.tpm <- merge(rsem.tpm, X04M_3[c(1,6)], by="gene_id")
#rsem.tpm <- merge(rsem.tpm, X13F_1[c(1,6)], by="gene_id")
#rsem.tpm <- merge(rsem.tpm, X13F_2[c(1,6)], by="gene_id")
#rsem.tpm <- merge(rsem.tpm, X13F_3[c(1,6)], by="gene_id")
#rsem.tpm <- merge(rsem.tpm, X13M_1[c(1,6)], by="gene_id")
#rsem.tpm <- merge(rsem.tpm, X13M_2[c(1,6)], by="gene_id")
#rsem.tpm <- merge(rsem.tpm, X13M_3[c(1,6)], by="gene_id")
#rsem.tpm <- merge(rsem.tpm, X13M_4[c(1,6)], by="gene_id")
#rsem.tpm <- merge(rsem.tpm, X15F_1[c(1,6)], by="gene_id")
#rsem.tpm <- merge(rsem.tpm, X15F_2[c(1,6)], by="gene_id")
#rsem.tpm <- merge(rsem.tpm, X15F_3[c(1,6)], by="gene_id")
#rsem.tpm <- merge(rsem.tpm, X15M_1[c(1,6)], by="gene_id")
#rsem.tpm <- merge(rsem.tpm, X15M_2[c(1,6)], by="gene_id")
#rsem.tpm <- merge(rsem.tpm, X15M_3[c(1,6)], by="gene_id")
#rsem.tpm <- merge(rsem.tpm, X21F_1[c(1,6)], by="gene_id")
#rsem.tpm <- merge(rsem.tpm, X21F_2[c(1,6)], by="gene_id")
#rsem.tpm <- merge(rsem.tpm, X21F_3[c(1,6)], by="gene_id")
#rsem.tpm <- merge(rsem.tpm, X21M_1[c(1,6)], by="gene_id")
#rsem.tpm <- merge(rsem.tpm, X21M_2[c(1,6)], by="gene_id")
#rsem.tpm <- merge(rsem.tpm, X21M_3[c(1,6)], by="gene_id")
#rsem.tpm <- merge(rsem.tpm, X21M_4[c(1,6)], by="gene_id")
#write.csv(rsem.tpm,"output/diff_expr/rsem.tpm.csv", row.names=FALSE)

# import the file
rsem.tpm <- read_delim("output/diff_expr/rsem.tpm.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE)
rsem.tpm <- rsem.tpm[c(-1)]

# differentially expressed genes -- from Andrés' rerun of Isabelle's script (received 03.11.20, rerun completed 09.11.20)

dt_df <- read_delim("output/diff_expr/dt_df.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE) # complete list of contrasts
de.over.B.males.genes <- read_delim("output/diff_expr/over.Bmales.vs.all.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE)
de.under.B.males.genes <- as.data.frame(c("g4126","g8486","g14550")) # i didn't bother extracting the file, it's just three genes

de.B.males.vs.nonB.males <- read_delim("output/diff_expr/B.males.vs.nonB.males.de.treat.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE)
de.B.males.vs.B.females <- read_delim("output/diff_expr/B.males.vs.B.females.de.treat.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE)
de.B.males.vs.nonB.females <- read_delim("output/diff_expr/B.males.vs.nonB.females.de.treat.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE)
de.B.females.vs.nonB.females <- read_delim("output/diff_expr/B.females.vs.nonB.females.de.treat.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE)
de.nonB.males.vs.nonB.females <- read_delim("output/diff_expr/nonB.males.vs.nonB.females.de.treat.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE)

de.over.B.genes<- read_delim("output/diff_expr/over.B.vs.nonB.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE)
de.under.B.genes<- read_delim("output/diff_expr/under.B.vs.nonB.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE)
###

##### Genes located on B candidate scaffolds

genes.AB <- merge(freeze.v0.genes.anno,scaffolds.final.assignment[c("seq","length","b.status.final","cov.04v13","b.status","b.status.asn","b.status.kmer")],by="seq")
genes.A <-  genes.AB[genes.AB$b.status.final == "A",]
genes.B1 <- genes.AB[genes.AB$b.status.final == "B1",]
genes.B2 <- genes.AB[genes.AB$b.status.final == "B2",]
genes.B3 <- genes.AB[genes.AB$b.status.final == "B3",]
table(genes.AB$b.status.final)
genes.in.Bs.anno <- rbind(genes.B1[genes.B1$anno == "Y",], genes.B2[genes.B2$anno == "Y",], genes.B3[genes.B3$anno == "Y",])[c(2,9,1,3,4,5,6)]
table(genes.in.Bs.anno$b.status.final)
#write.table(genes.in.Bs.anno, file = "output/genes.in.Bs.anno.csv",row.names = F,sep = ",")

###

##### Examine genes that are overexpressed in B males compared to all other samples (non-B males and both B/non-B females)

de.over.B.males.genes.anno <- left_join(de.over.B.males.genes,genes.AB,by="gene")
table(de.over.B.males.genes.anno$b.status.final)
#write.csv(de.over.B.males.genes.anno,"output/diff_expr/over.Bmales.vs.all.anno.csv")

colnames(de.under.B.males.genes)[1] <- "gene"
de.under.B.males.genes.anno <- left_join(de.under.B.males.genes,genes.AB,by="gene")
#write.csv(de.under.B.males.genes.anno,"output/diff_expr/under.Bmales.vs.all.anno.csv")

###

##### Examine genes that are over/under expressed in B compared to non B lines
colnames(de.over.B.genes)[2] <- "gene"
de.over.B.genes.anno <- left_join(de.over.B.genes,genes.AB,by="gene")
table(de.over.B.genes.anno$b.status.final)
write.csv(de.over.B.genes.anno,"output/diff_expr/over.B.vs.noB.csv")

colnames(de.under.B.genes)[2] <- "gene"
de.under.B.genes.anno <- left_join(de.under.B.genes,genes.AB,by="gene")
table(de.under.B.genes.anno$b.status.final)
write.csv(de.under.B.genes.anno,"output/diff_expr/under.B.vs.noB.csv")

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
table(de.over.B.males.vs.nonB.males.anno$b.status.final)
de.over.B.males.vs.nonB.males.anno.yes <- de.over.B.males.vs.nonB.males.anno[de.over.B.males.vs.nonB.males.anno$anno == "Y",]
table(de.over.B.males.vs.nonB.males.anno.yes$b.status.final)
write.csv(de.over.B.males.vs.nonB.males.anno,"output/diff_expr/de.over.B.males.vs.nonB.males.anno.csv")

de.under.B.males.vs.nonB.males.anno <- de.B.males.vs.nonB.males.anno[de.B.males.vs.nonB.males.anno$logFC < 1,]
table(de.under.B.males.vs.nonB.males.anno$b.status.final)
de.under.B.males.vs.nonB.males.anno.yes <- de.under.B.males.vs.nonB.males.anno[de.under.B.males.vs.nonB.males.anno$anno == "Y",]
table(de.under.B.males.vs.nonB.males.anno.yes$b.status.final)
write.csv(de.under.B.males.vs.nonB.males.anno,"output/diff_expr/de.under.B.males.vs.nonB.males.anno.csv")

## B females vs non B females

de.over.B.females.vs.nonB.females.anno <- de.B.females.vs.nonB.females.anno[de.B.females.vs.nonB.females.anno$logFC > 1,]
table(de.over.B.females.vs.nonB.females.anno$b.status.final)
de.over.B.females.vs.nonB.females.anno.yes <- de.over.B.females.vs.nonB.females.anno[de.over.B.females.vs.nonB.females.anno$anno == "Y",]
table(de.over.B.females.vs.nonB.females.anno.yes$b.status.final)
write.csv(de.over.B.females.vs.nonB.females.anno,"output/diff_expr/de.over.B.females.vs.nonB.females.anno.csv")

de.under.B.females.vs.nonB.females.anno <- de.B.females.vs.nonB.females.anno[de.B.females.vs.nonB.females.anno$logFC < 1,]
table(de.under.B.females.vs.nonB.females.anno$b.status.final)
de.under.B.females.vs.nonB.females.anno.yes <- de.under.B.females.vs.nonB.females.anno[de.under.B.females.vs.nonB.females.anno$anno == "Y",]
table(de.under.B.females.vs.nonB.females.anno.yes$b.status.final)
write.csv(de.under.B.females.vs.nonB.females.anno,"output/diff_expr/de.under.B.females.vs.nonB.females.anno.csv")

## B females vs non B females

de.B.males.vs.nonB.female.anno <- left_join(de.B.males.vs.nonB.females,genes.AB,by="gene")
de.B.females.vs.nonB.females.anno <- left_join(de.B.females.vs.nonB.females,genes.AB,by="gene")
de.nonB.males.vs.nonB.females.anno <- left_join(de.nonB.males.vs.nonB.females,genes.AB,by="gene")

## remaining comparisons

de.B.males.vs.B.females.anno <- left_join(de.B.males.vs.B.females,genes.AB,by="gene")
de.B.males.vs.nonB.females.anno <- left_join(de.B.males.vs.nonB.females,genes.AB,by="gene")
de.nonB.males.vs.nonB.females.anno <- left_join(de.nonB.males.vs.nonB.females,genes.AB,by="gene")

write.csv(de.B.males.vs.B.females.anno,"output/diff_expr/de.B.males.vs.B.females.anno.csv")
write.csv(de.B.males.vs.nonB.females.anno,"output/diff_expr/de.B.males.vs.nonB.females.anno.csv")
write.csv(de.nonB.males.vs.nonB.females.anno,"output/diff_expr/de.nonB.males.vs.nonB.females.anno.csv")






# extract putative B genes that are differentially expressed between B-carrying males and females
de.B.males.vs.B.females.anno.b <- de.B.males.vs.B.females.anno[de.B.males.vs.B.females.anno$b.status.final != "A",]
write.csv(de.B.males.vs.B.females.anno.b,"output/diff_expr/de.B.males.vs.B.females.anno.b.csv")






###

##### Re-evaluating genes in B scaffolds

# merge B assignment with expression data
rsem.tpm.B.males <- subset(rsem.tpm[c(1,5,6,7,11,12,13,14)])
rsem.tpm.B.females <- subset(rsem.tpm[c(1,2,3,4,8,9,10)])
rsem.tpm.nonB.males <- subset(rsem.tpm[c(1,18,19,20,24,25,26,27)])
rsem.tpm.nonB.females <- subset(rsem.tpm[c(1,15,16,17,21,22,23)])

# B vs non B groups
rsem.tpm.B <- subset(rsem.tpm[c(1,5,6,7,11,12,13,14,2,3,4,8,9,10)])
rsem.tpm.nonB <- subset(rsem.tpm[c(1,18,19,20,24,25,26,27,15,16,17,21,22,23)])



gene <- rsem.tpm[1]
B.males.tpm <- rowMeans(rsem.tpm.B.males[-1])
B.females.tpm <- rowMeans(rsem.tpm.B.females[-1])
nonB.males.tpm <- rowMeans(rsem.tpm.nonB.males[-1])
nonB.females.tpm <- rowMeans(rsem.tpm.nonB.females[-1])

B.tpm <- rowMeans(rsem.tpm.B[-1])
B.tpm <- rowMeans(rsem.tpm.nonB[-1])



rsem.tpm.avg <- data.frame(gene,B.males.tpm,B.females.tpm,nonB.males.tpm,nonB.females.tpm)
colnames(rsem.tpm.avg)[1] <- "gene"
genes.AB.tpm <- merge(rsem.tpm.avg, genes.AB, by ="gene")

#for b vs non b
rsem.tpm.avg2 <- data.frame(gene,B.tpm,B.tpm)
colnames(rsem.tpm.avg2)[1] <- "gene"
genes.AB.tpm2 <- merge(rsem.tpm.avg2, genes.AB, by ="gene")




library(patchwork)

genes.AB.tpm


b.males.tpm <- ggplot(genes.AB.tpm, aes(b.status.final, log10(B.males.tpm+1e-3),fill=b.status.final)) + 
  geom_jitter(data=genes.AB.tpm, aes(b.status.final, log10(B.males.tpm+1e-3), group=b.status.final),
              size=0.5, width = 0.4, alpha=0.5, show.legend=FALSE) +
  geom_boxplot(alpha=0.75,outlier.shape = NA,notch=TRUE,lwd=0.6) +
  labs(title="B males",x="", y ="log10(Average TPM + 1e-4)") +
  scale_fill_manual(values=c("gray85","royalblue4", "deepskyblue", "cadetblue", "lavenderblush4")) +
  theme_bw()
#b.males.tpm

b.females.tpm <- ggplot(genes.AB.tpm, aes(b.status.final, log10(B.females.tpm+1e-3),fill=b.status.final)) + 
  geom_jitter(data=genes.AB.tpm, aes(b.status.final, log10(B.females.tpm+1e-3), group=b.status.final),
              size=0.5, width = 0.4, alpha=0.5, show.legend=FALSE) +
  geom_boxplot(alpha=0.75,outlier.shape = NA,notch=TRUE,lwd=0.6) +
  labs(title="B females",x="", y ="log10(Average TPM + 1e-4)") +
  scale_fill_manual(values=c("gray85","royalblue4", "deepskyblue", "cadetblue", "lavenderblush4")) +
  theme_bw()

nonb.males.tpm <- ggplot(genes.AB.tpm, aes(b.status.final, log10(nonB.males.tpm+1e-3),fill=b.status.final)) + 
  geom_jitter(data=genes.AB.tpm, aes(b.status.final, log10(nonB.males.tpm+1e-3), group=b.status.final),
              size=0.5, width = 0.4, alpha=0.5, show.legend=FALSE) +
  geom_boxplot(alpha=0.75,outlier.shape = NA,notch=TRUE,lwd=0.6) +
  labs(title="Non B males",x="", y ="log10(Average TPM + 1e-4)") +
  scale_fill_manual(values=c("gray85","royalblue4", "deepskyblue", "cadetblue", "lavenderblush4")) +
  theme_bw()

nonb.females.tpm <- ggplot(genes.AB.tpm, aes(b.status.final, log10(nonB.females.tpm+1e-3),fill=b.status.final)) + 
  geom_jitter(data=genes.AB.tpm, aes(b.status.final, log10(nonB.females.tpm+1e-3), group=b.status.final),
              size=0.5, width = 0.4, alpha=0.5, show.legend=FALSE) +
  geom_boxplot(alpha=0.75,outlier.shape = NA,notch=TRUE,lwd=0.6) +
  labs(title="Non B females",x="", y ="log10(Average TPM + 1e-4)",guide="Scaffolds") +
  scale_fill_manual(values=c("gray85","royalblue4", "deepskyblue", "cadetblue", "lavenderblush4")) +
  theme_bw()

#b vs non b
b.tpm <- ggplot(genes.AB.tpm2, aes(b.status.final, log10(B.tpm+1e-3),fill=b.status.final)) + 
  geom_jitter(data=genes.AB.tpm, aes(b.status.final, log10(B.tpm+1e-3), group=b.status.final),
              size=0.5, width = 0.4, alpha=0.5, show.legend=FALSE) +
  geom_boxplot(alpha=0.75,outlier.shape = NA,notch=TRUE,lwd=0.6) +
  labs(title="B+ lines",x="", y ="log10(Average TPM + 1e-4)") +
  scale_fill_manual(values=c("gray85","royalblue4", "deepskyblue", "cadetblue", "lavenderblush4")) +
  theme_bw()


nonb.tpm <- ggplot(genes.AB.tpm2, aes(b.status.final, log10(nonB.tpm+1e-3),fill=b.status.final)) + 
  geom_jitter(data=genes.AB.tpm2, aes(b.status.final, log10(nonB.tpm+1e-3), group=b.status.final),
              size=0.5, width = 0.4, alpha=0.5, show.legend=FALSE) +
  geom_boxplot(alpha=0.75,outlier.shape = NA,notch=TRUE,lwd=0.6) +
  labs(title="B- lines",x="", y ="log10(Average TPM + 1e-4)") +
  scale_fill_manual(values=c("gray85","royalblue4", "deepskyblue", "cadetblue", "lavenderblush4")) +
  theme_bw()


library(patchwork)
jpeg("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/misc/gene.tpm.by.b.status.jpeg",
     width = 4200, height = 3800, units = 'px', res = 300)
b.males.tpm + b.females.tpm + nonb.males.tpm + nonb.females.tpm + b.tpm + nonb.tpm
dev.off()

# examine classes of genes

genes.B1.tpm <- genes.AB.tpm[genes.AB.tpm$b.status.final == "B1",]
nrow(genes.B1.tpm)
genes.B1.tpm.dt <- left_join(genes.B1.tpm[1:14], dt_df[-c(1)], by ="gene")
write.csv(genes.B1.tpm.dt,"output/diff_expr/genes.B1.tpm.dt.csv")

genes.B2.tpm <- genes.AB.tpm[genes.AB.tpm$b.status.final == "B2",]
nrow(genes.B2.tpm)
genes.B2.tpm.dt <- left_join(genes.B2.tpm[1:14], dt_df[-c(1)], by ="gene")
write.csv(genes.B2.tpm.dt,"output/diff_expr/genes.B2.tpm.dt.csv")

genes.B3.tpm <- genes.AB.tpm[genes.AB.tpm$b.status.final == "B3",]
nrow(genes.B3.tpm)
genes.B3.tpm.dt <- left_join(genes.B3.tpm[1:14], dt_df[-c(1)], by ="gene")
write.csv(genes.B3.tpm.dt,"output/diff_expr/genes.B3.tpm.dt.csv")

genes.A.tpm <- genes.AB.tpm[genes.AB.tpm$b.status.final == "A",]
nrow(genes.A.tpm.dt)
genes.A.tpm.dt <- left_join(genes.A.tpm[1:14], dt_df[-c(1)], by ="gene")
write.csv(genes.A.tpm.dt,"output/diff_expr/genes.A.tpm.dt.csv")
head(genes.A.tpm.dt)

