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

# master anno
freeze.v0.genes.anno <- read_delim("output/freeze.v0.genes.anno.complete.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE)

# B chromosome assignment
# read.table("output/genes.by.scaffolds.tsv", sep = "\t", header = T)
scaffolds.final.assignment <- read_delim("output/scaffolds.final.assignment.tsv","\t", escape_double = FALSE, col_names = T,trim_ws = TRUE)
scaffolds.final.assignment[scaffolds.final.assignment$seq %in% c("scaffold_360", "scaffold_957"), 'b.status.final'] <- 'B-A'

table(scaffolds.final.assignment$b.status.final)

# import the file
rsem.tpm <- read_delim("output/B_diff_expr/alternative_logFC1/rsem.tpm.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE)
rsem.tpm <- rsem.tpm[c(-1)]
rsem.tpm

##### Genes located on B candidate scaffolds

genes.AB <- merge(freeze.v0.genes.anno,scaffolds.final.assignment[c("seq","length","b.status.final","cov.b.status","b.status.asn","b.status.kmer")],by="seq")
genes.AB$anno <- ifelse(!is.na(genes.AB$blast) | !is.na(genes.AB$diamond) | !is.na(genes.AB$pfam_acc) | !is.na(genes.AB$ipr_acc) | !is.na(genes.AB$GO),"Y","N")
genes.AB$b.status.final <- factor(genes.AB$b.status.final, levels = c('A', 'B-A', 'B', 'Bc'))

genes.A <-  genes.AB[genes.AB$b.status.final == "A",]
genes.B1 <- genes.AB[genes.AB$b.status.final == "B" | genes.AB$b.status.final == "B-A",]
genes.B2 <- genes.AB[genes.AB$b.status.final == "Bc",]

# which or these are annotated?
genes.AB$anno <- ifelse(!is.na(genes.AB$blast) | !is.na(genes.AB$diamond) | !is.na(genes.AB$pfam_acc) | !is.na(genes.AB$ipr_acc) | !is.na(genes.AB$GO),"Y","N")
genes.in.Bs.anno <- rbind(genes.B1[genes.B1$anno == "Y",], genes.B2[genes.B2$anno == "Y",])[c(2,1,3,4,5,6,7,8,9,10,13)]
table(genes.in.Bs.anno$b.status.final)

# table of annotated genes
write.table(genes.in.Bs.anno, file = "output/genes.in.Bs.anno.csv",row.names = F,sep = ",")

sum(genes.B1$gene_len)
length(unique(genes.B1$seq))

# which or these are expressed?
tpm <- rsem.tpm

tpm$tpm04F <- (tpm$X04F_1+tpm$X04F_2+tpm$X04F_3)/3
tpm$tpm04M <- (tpm$X04M_1+tpm$X04M_2+tpm$X04M_3)/3
tpm$tpm13F <- (tpm$X13F_1+tpm$X13F_2+tpm$X13F_3)/3
tpm$tpm13M <- (tpm$X13M_1+tpm$X13M_2+tpm$X13M_3+tpm$X13M_4)/4
tpm$tpm15F <- (tpm$X15F_1+tpm$X15F_2+tpm$X15F_3)/3
tpm$tpm15M <- (tpm$X15M_1+tpm$X15M_2+tpm$X15M_3)/3
tpm$tpm21F <- (tpm$X21F_1+tpm$X21F_2+tpm$X21F_3)/3
tpm$tpm21M <- (tpm$X21M_1+tpm$X21M_2+tpm$X21M_3+tpm$X21M_4)/4

tpm <- tpm[c("gene_id","tpm04F","tpm04M","tpm13F","tpm13M","tpm15F","tpm15M","tpm21F","tpm21M")]
colnames(tpm)[1] <- "gene"
genes.B1.tpm <- merge(genes.B1,tpm,by="gene")
genes.B2.tpm <- merge(genes.B2,tpm,by="gene")

genes.B.tpm <- rbind(genes.B1.tpm,genes.B2.tpm)
#write.table(genes.B.tpm, file = "output/genes.in.Bs.tpm.csv",row.names = F,sep = ",")

genes.B1.tpm$expr <- "no"
genes.B1.tpm$expr <- ifelse((genes.B1.tpm$tpm04F > 1 | genes.B1.tpm$tpm04M > 1 | genes.B1.tpm$tpm13F > 1 | genes.B1.tpm$tpm13M > 1 | genes.B1.tpm$tpm15F > 1 | genes.B1.tpm$tpm15M > 1 | genes.B1.tpm$tpm21F > 1 | genes.B1.tpm$tpm21M > 1), "yes_1", genes.B1.tpm$expr)
genes.B1.tpm$expr <- ifelse((genes.B1.tpm$tpm04F > 10 | genes.B1.tpm$tpm04M > 10 | genes.B1.tpm$tpm13F > 10 | genes.B1.tpm$tpm13M > 10 | genes.B1.tpm$tpm15F > 10 | genes.B1.tpm$tpm15M > 10 | genes.B1.tpm$tpm21F > 10 | genes.B1.tpm$tpm21M > 10), "yes_10", genes.B1.tpm$expr)
table(genes.B1.tpm$expr, genes.B1.tpm$b.status.final)
genes.B1.tpm[(genes.B1.tpm$expr == "yes_10"),]

# differences in expression rates across genes

rsem.tpm.B.males <- subset(rsem.tpm[c(1,5,6,7,11,12,13,14)])
rsem.tpm.B.females <- subset(rsem.tpm[c(1,2,3,4,8,9,10)])
rsem.tpm.nonB.males <- subset(rsem.tpm[c(1,18,19,20,24,25,26,27)])
rsem.tpm.nonB.females <- subset(rsem.tpm[c(1,15,16,17,21,22,23)])

# B vs non B groups
rsem.tpm.B <- subset(rsem.tpm[c(1,5,6,7,11,12,13,14,2,3,4,8,9,10)])
rsem.tpm.nonB <- subset(rsem.tpm[c(1,18,19,20,24,25,26,27,15,16,17,21,22,23)])

# differences in expression rates across B+ lines and B- lines

gene <- rsem.tpm[1]
B.males.tpm <- rowMeans(rsem.tpm.B.males[-1])
B.females.tpm <- rowMeans(rsem.tpm.B.females[-1])
nonB.males.tpm <- rowMeans(rsem.tpm.nonB.males[-1])
nonB.females.tpm <- rowMeans(rsem.tpm.nonB.females[-1])

B.tpm <- rowMeans(rsem.tpm.B[-1])
nonB.tpm <- rowMeans(rsem.tpm.nonB[-1])

rsem.tpm.avg <- data.frame(gene,B.males.tpm,B.females.tpm,nonB.males.tpm,nonB.females.tpm)
colnames(rsem.tpm.avg)[1] <- "gene"
genes.AB.tpm <- merge(rsem.tpm.avg, genes.AB, by ="gene")

#for b vs non b
rsem.tpm.avg2 <- data.frame(gene,B.tpm,nonB.tpm)
colnames(rsem.tpm.avg2)[1] <- "gene"
genes.AB.tpm2 <- merge(rsem.tpm.avg2, genes.AB, by ="gene")

pal <- c("gray85","lightpink2", "royalblue4", "deepskyblue")

b.males.tpm <- ggplot(genes.AB.tpm, aes(b.status.final, log10(B.males.tpm+1e-3), fill=b.status.final)) +
  geom_jitter(data=genes.AB.tpm, aes(b.status.final, log10(B.males.tpm+1e-3), group=b.status.final),
              size=0.5, width = 0.4, alpha=0.5, show.legend=FALSE) +
  geom_boxplot(alpha=0.75,outlier.shape = NA,notch=TRUE,lwd=0.6) +
  labs(title="B males",x="", y ="log10(Average TPM + 1e-4)") +
  scale_fill_manual(values=pal) +
  theme_bw()

b.females.tpm <- ggplot(genes.AB.tpm, aes(b.status.final, log10(B.females.tpm+1e-3),fill=b.status.final)) +
  geom_jitter(data=genes.AB.tpm, aes(b.status.final, log10(B.females.tpm+1e-3), group=b.status.final),
              size=0.5, width = 0.4, alpha=0.5, show.legend=FALSE) +
  geom_boxplot(alpha=0.75,outlier.shape = NA,notch=TRUE,lwd=0.6) +
  labs(title="B females",x="", y ="log10(Average TPM + 1e-4)") +
  scale_fill_manual(values=pal) +
  theme_bw()

nonb.males.tpm <- ggplot(genes.AB.tpm, aes(b.status.final, log10(nonB.males.tpm+1e-3),fill=b.status.final)) +
  geom_jitter(data=genes.AB.tpm, aes(b.status.final, log10(nonB.males.tpm+1e-3), group=b.status.final),
              size=0.5, width = 0.4, alpha=0.5, show.legend=FALSE) +
  geom_boxplot(alpha=0.75,outlier.shape = NA,notch=TRUE,lwd=0.6) +
  labs(title="Non B males",x="", y ="log10(Average TPM + 1e-4)") +
  scale_fill_manual(values=pal) +
  theme_bw()

nonb.females.tpm <- ggplot(genes.AB.tpm, aes(b.status.final, log10(nonB.females.tpm+1e-3),fill=b.status.final)) +
  geom_jitter(data=genes.AB.tpm, aes(b.status.final, log10(nonB.females.tpm+1e-3), group=b.status.final),
              size=0.5, width = 0.4, alpha=0.5, show.legend=FALSE) +
  geom_boxplot(alpha=0.75,outlier.shape = NA,notch=TRUE,lwd=0.6) +
  labs(title="Non B females",x="", y ="",guide="Scaffolds") +
  scale_fill_manual(values=pal) +
  theme_bw()

# are expression levels different across B groups?

# library("PMCMR")
genes.AB.tpm2.noA <- genes.AB.tpm2[genes.AB.tpm2$b.status.final != "A",]

genes.AB.tpm2.noA$b.status.final <- as.factor(genes.AB.tpm2.noA$b.status.final)
kruskal.test(B.tpm~b.status.final, data=genes.AB.tpm2.noA)
# posthoc.kruskal.nemenyi.test(B.tpm~b.status.final, data=genes.AB.tpm2.noA,dist="Chisquare")
kruskal.test(nonB.tpm~b.status.final, data=genes.AB.tpm2.noA)
# posthoc.kruskal.nemenyi.test(nonB.tpm~b.status.final, data=genes.AB.tpm2.noA,dist="Chisquare")

###########
## B lines expression boxplot
###########
b.tpm.plot <- ggplot(data = genes.AB.tpm2, aes(b.status.final, log10(B.tpm+1e-3),fill=b.status.final)) +
  geom_jitter(data = genes.AB.tpm2, aes(b.status.final, log10(B.tpm+1e-3), group=b.status.final),
              size=0.5, width = 0.4, alpha=0.5, show.legend=FALSE) +
  geom_boxplot(alpha=0.75,outlier.shape = NA,notch=TRUE,lwd=0.6) +
  labs(title="Gene expression levels in B+ samples",x="", y ="Mean TPM (log10)", fill = "Location") +
  scale_fill_manual(values=pal) + theme_bw()
  # geom_text(x=2.5,y=2.40, label="1.2e-12",size=3) +
  # geom_text(x=3.0,y=2.80, label="5.9e-06",size=3) +
  # geom_text(x=3.5,y=3.30, label="0.15" ,size=3) +
  # geom_segment(aes(x=2,xend=3,y=2.30,yend=2.30),size=0.1) +
  # geom_segment(aes(x=2,xend=4,y=2.70,yend=2.70),size=0.1) +
  # geom_segment(aes(x=3,xend=4,y=3.20,yend=3.20),size=0.1)

###########
## NO B lines expression boxplot
###########
nonb.tpm.plot <- ggplot(data = genes.AB.tpm2, aes(b.status.final, log10(nonB.tpm+1e-3),fill=b.status.final)) +
  geom_jitter(data = genes.AB.tpm2, aes(b.status.final, log10(nonB.tpm+1e-3), group=b.status.final),
              size=0.5, width = 0.4, alpha=0.5, show.legend=FALSE) +
  geom_boxplot(alpha=0.75,outlier.shape = NA,notch=TRUE,lwd=0.6) +
  labs(title="Gene expression levels in B- samples",x="", y ="Mean TPM (log10)", fill = "Location") +
  scale_fill_manual(values=pal) + theme_bw()

tiff("manuscript/figures_revision/supplementary_fig3_expression_B.tiff", width = 8, height = 8, units = 'in', res = 300)
    b.tpm.plot
dev.off()

tiff("manuscript/figures_revision/supplementary_fig3_expression_noB.tiff", width = 8, height = 8, units = 'in', res = 300)
    nonb.tpm.plot
dev.off()

# nonb.tpm.plot <- ggplot(genes.AB.tpm2, aes(b.status.final, log10(nonB.tpm+1e-3),fill=factor(b.status.final, levels = c('A', 'B-A', 'Bc', 'B'))) +
#   geom_jitter(data=genes.AB.tpm2, aes(b.status.final, log10(nonB.tpm+1e-3), group=factor(b.status.final, levels = c('A', 'B-A', 'Bc', 'B')),
#               size=0.5, width = 0.4, alpha=0.5, show.legend=FALSE) +
#   geom_boxplot(alpha=0.75,outlier.shape = NA,notch=TRUE,lwd=0.6) +
#   labs(title="Gene expression levels in B- samples",x="", y ="Mean TPM (log10)", fill = "Location") +
#   scale_fill_manual(values=c("gray85", "lightpink2", "royalblue4", "deepskyblue", "lightpink2", "lavenderblush4")) +
#   geom_text(x=2.5,y=2.40, label="< 2e-16",size=3) +
#   geom_text(x=3.0,y=2.80, label="5.3e-10",size=3) +
#   geom_text(x=3.5,y=3.30, label="0.067" ,size=3) +
#   geom_segment(aes(x=2,xend=3,y=2.30,yend=2.30),size=0.1) +
#   geom_segment(aes(x=2,xend=4,y=2.70,yend=2.70),size=0.1) +
#   geom_segment(aes(x=3,xend=4,y=3.20,yend=3.20),size=0.1) +
#   theme_bw()

# B_scf_with_genes <- names(table(genes.AB.tpm2[genes.AB.tpm2$b.status.final == 'B', 'seq']))
# getTspSummary <- function(scf){ summary(log10(genes.AB.tpm2[genes.AB.tpm2$seq == scf, 'nonB.tpm'] + 1e-3)) }
# B_per_scf_exp <- t(sapply(B_scf_with_genes, getTspSummary))
# hist(log10(genes.AB.tpm2[genes.AB.tpm2$seq == 'scaffold_360', 'nonB.tpm'] + 1e-3))
# scaffold_360 looks FUNKY!

b.tpm.plot + nonb.tpm.plot

# differentially expressed genes -- from Andres' rerun of Isabelle's script (received 03.11.20, rerun completed 09.11.20)
# Kamil revisions 14. 04. 2023

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
de.B.vs.nonB <- read_delim("output/diff_expr/B.vs.nonB.de.treat.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE)

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
write.csv(de.over.B.genes.anno,"output/diff_expr/over.B.vs.noB.anno.csv")

colnames(de.under.B.genes)[2] <- "gene"
de.under.B.genes.anno <- left_join(de.under.B.genes,genes.AB,by="gene")


table(de.under.B.genes.anno$b.status.final)
write.csv(de.under.B.genes.anno,"output/diff_expr/under.B.vs.noB.anno.csv")

##### Examine differentially expressed genes between groups

de.B.males.vs.nonB.males <- de.B.males.vs.nonB.males[c(-1)]
de.B.males.vs.B.females <- de.B.males.vs.B.females[c(-1)]
de.B.males.vs.nonB.females <- de.B.males.vs.nonB.females[c(-1)]
de.B.females.vs.nonB.females <- de.B.females.vs.nonB.females[c(-1)]
de.nonB.males.vs.nonB.females <- de.nonB.males.vs.nonB.females[c(-1)]
de.B.vs.nonB <- de.B.vs.nonB[c(-1)]


colnames(de.B.males.vs.nonB.males)[1] <- "gene"
colnames(de.B.males.vs.B.females)[1] <- "gene"
colnames(de.B.males.vs.nonB.females)[1] <- "gene"
colnames(de.B.females.vs.nonB.females)[1] <- "gene"
colnames(de.nonB.males.vs.nonB.females)[1] <- "gene"
colnames(de.B.vs.nonB)[1] <- "gene"


de.B.males.vs.nonB.males.anno <- left_join(de.B.males.vs.nonB.males,genes.AB,by="gene")
de.B.males.vs.B.females.anno <- left_join(de.B.males.vs.B.females,genes.AB,by="gene")
de.B.males.vs.nonB.female.anno <- left_join(de.B.males.vs.nonB.females,genes.AB,by="gene")
de.B.females.vs.nonB.females.anno <- left_join(de.B.females.vs.nonB.females,genes.AB,by="gene")
de.nonB.males.vs.nonB.females.anno <- left_join(de.nonB.males.vs.nonB.females,genes.AB,by="gene")
de.B.vs.nonB.anno <- left_join(de.B.vs.nonB,genes.AB,by="gene")

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

## B lines vs non B lines
de.over.B.vs.nonB.anno <- de.B.vs.nonB.anno[de.B.vs.nonB.anno$logFC > 1,]
table(de.over.B.vs.nonB.anno$b.status.final)
de.over.B.vs.nonB.anno.yes <- de.over.B.vs.nonB.anno[de.over.B.vs.nonB.anno$anno == "Y",]
table(de.over.B.vs.nonB.anno.yes$b.status.final)
write.csv(de.over.B.vs.nonB.anno,"output/diff_expr/de.over.B.vs.nonB.anno.csv")

de.under.B.vs.nonB.anno <- de.B.vs.nonB.anno[de.B.vs.nonB.anno$logFC < 1,]
table(de.under.B.vs.nonB.anno$b.status.final)
de.under.B.vs.nonB.anno.yes <- de.under.B.vs.nonB.anno[de.under.B.vs.nonB.anno$anno == "Y",]
table(de.under.B.vs.nonB.anno.yes$b.status.final)
write.csv(de.under.B.vs.nonB.anno,"output/diff_expr/de.under.B.vs.nonB.anno.csv")




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
  scale_fill_manual(values=c("gray85","royalblue4", "deepskyblue", "lightpink2", "lavenderblush4")) +
  theme_bw()
#b.males.tpm

b.females.tpm <- ggplot(genes.AB.tpm, aes(b.status.final, log10(B.females.tpm+1e-3),fill=b.status.final)) +
  geom_jitter(data=genes.AB.tpm, aes(b.status.final, log10(B.females.tpm+1e-3), group=b.status.final),
              size=0.5, width = 0.4, alpha=0.5, show.legend=FALSE) +
  geom_boxplot(alpha=0.75,outlier.shape = NA,notch=TRUE,lwd=0.6) +
  labs(title="B females",x="", y ="log10(Average TPM + 1e-4)") +
  scale_fill_manual(values=c("gray85","royalblue4", "deepskyblue", "lightpink2", "lavenderblush4")) +
  theme_bw()

nonb.males.tpm <- ggplot(genes.AB.tpm, aes(b.status.final, log10(nonB.males.tpm+1e-3),fill=b.status.final)) +
  geom_jitter(data=genes.AB.tpm, aes(b.status.final, log10(nonB.males.tpm+1e-3), group=b.status.final),
              size=0.5, width = 0.4, alpha=0.5, show.legend=FALSE) +
  geom_boxplot(alpha=0.75,outlier.shape = NA,notch=TRUE,lwd=0.6) +
  labs(title="Non B males",x="", y ="log10(Average TPM + 1e-4)") +
  scale_fill_manual(values=c("gray85","royalblue4", "deepskyblue", "lightpink2", "lavenderblush4")) +
  theme_bw()

nonb.females.tpm <- ggplot(genes.AB.tpm, aes(b.status.final, log10(nonB.females.tpm+1e-3),fill=b.status.final)) +
  geom_jitter(data=genes.AB.tpm, aes(b.status.final, log10(nonB.females.tpm+1e-3), group=b.status.final),
              size=0.5, width = 0.4, alpha=0.5, show.legend=FALSE) +
  geom_boxplot(alpha=0.75,outlier.shape = NA,notch=TRUE,lwd=0.6) +
  labs(title="Non B females",x="", y ="log10(Average TPM + 1e-4)",guide="Scaffolds") +
  scale_fill_manual(values=c("gray85","royalblue4", "deepskyblue", "lightpink2", "lavenderblush4")) +
  theme_bw()

#b vs non b
b.tpm <- ggplot(genes.AB.tpm2, aes(b.status.final, log10(B.tpm+1e-3),fill=b.status.final)) +
  geom_jitter(data=genes.AB.tpm, aes(b.status.final, log10(B.tpm+1e-3), group=b.status.final),
              size=0.5, width = 0.4, alpha=0.5, show.legend=FALSE) +
  geom_boxplot(alpha=0.75,outlier.shape = NA,notch=TRUE,lwd=0.6) +
  labs(title="B+ lines",x="", y ="log10(Average TPM + 1e-4)") +
  scale_fill_manual(values=c("gray85","royalblue4", "deepskyblue", "lightpink2", "lavenderblush4")) +
  theme_bw()


nonb.tpm <- ggplot(genes.AB.tpm2, aes(b.status.final, log10(nonB.tpm+1e-3),fill=b.status.final)) +
  geom_jitter(data=genes.AB.tpm2, aes(b.status.final, log10(nonB.tpm+1e-3), group=b.status.final),
              size=0.5, width = 0.4, alpha=0.5, show.legend=FALSE) +
  geom_boxplot(alpha=0.75,outlier.shape = NA,notch=TRUE,lwd=0.6) +
  labs(title="B- lines",x="", y ="log10(Average TPM + 1e-4)") +
  scale_fill_manual(values=c("gray85","royalblue4", "deepskyblue", "lightpink2", "lavenderblush4")) +
  theme_bw()


library(patchwork)
#jpeg("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/misc/gene.tpm.by.b.status.jpeg",
    # width = 4200, height = 3800, units = 'px', res = 300)
jpeg("misc/gene.tpm.by.b.status.jpeg",
     width = 4200, height = 3800, units = 'px', res = 300)
b.males.tpm + b.females.tpm + nonb.males.tpm + nonb.females.tpm
dev.off()
jpeg("misc/gene.tpm.by.b.statusnew.jpeg",
     width = 5400, height = 4200, units = 'px', res = 300)
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
