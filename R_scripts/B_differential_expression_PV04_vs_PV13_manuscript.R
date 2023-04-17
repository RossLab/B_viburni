#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("edgeR")
rm(list=ls())
library(edgeR)
library(methods)
library(dplyr)
library(tidyverse)
library(limma)
#BiocManager::install("Glimma")
library(Glimma)
library(gplots)
#BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
library(RColorBrewer)
#BiocManager::install("GOstats")
library(GOstats)
library(GSEABase)
#BiocManager::install("treemap")
library(treemap)
library(patchwork)

B_col = "royalblue4"
B_c_col = "deepskyblue"

# count file from all samples
# need a dataframe containing all gene info per gene id

freeze.v0.genes.anno <- read_delim("output/freeze.v0.genes.anno.complete.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE) # master anno
genes.by.scaffold <- read_delim("output/genes.by.scaffolds.tsv","\t", escape_double = FALSE, col_names = T,trim_ws = TRUE)
rsem.counts <-read.delim("R_scripts/RSEM_digi.counts.matrix",header=TRUE) #matrix generated from rsem


a <- rsem.counts
head(a,2)

colnames(a) <- substr(colnames(a),start=1,stop=6) #removing ".genes.results" in colnames
colnames(a)

#first column X is the gene id, this needs to be removed but has the order of gene id has to match when merging with gene info
head(a)
a1 <- a[,-1]
head(a1)
rownames(a1) <- a[,1]
head(a1)
a <- a1
head(a)

# we want just PV04 vs PV13 contrast, so I will remove the rest
a <- a[, 1:13]
replicate <- as.factor(c("X04F","X04F","X04F","X04M","X04M","X04M","X13F","X13F","X13F","X13M","X13M","X13M","X13M"))
x <- DGEList(counts=round(a), genes=rownames(a), group = replicate)
head(x)

unfilt_logcounts <- cpm(x, log=TRUE)

# rowSums(x$counts==0) # for each sample where the count is 0, then sums the samples that have counts = 0 for a gene
table(rowSums(x$counts==0) == 13) #number of genes that have all samples with 0 counts
proportion_0_count <- 2871 * 100 / nrow(x)
proportion_0_count

# compare filtering options. By default, the function keeps genes with about 10 read counts or more in a minimum number of samples, where the number of samples is chosen according to the minimum group sample size. The actual filtering uses CPM values rather than counts in order to avoid giving preference to samples with large library sizes.
# reduce min.count to 5 to account for the low expression of B genes
keep.exprs.group <- filterByExpr(x, group=x$samples$group,min.count=5)
# keep.exprs.group[keep.exprs.group == FALSE]
x1 <- x[keep.exprs.group, keep.lib.sizes=FALSE]
dim(x1)
x <- x1

x <- calcNormFactors(x, method = "TMM")
logcounts <- cpm(x, log=TRUE)

x$samples$norm.factors

# define groups by sex and B presence or absence
group1=c("PV04_F","PV04_F","PV04_F","PV04_M","PV04_M","PV04_M",
         "PV13_F","PV13_F","PV13_F","PV13_M","PV13_M","PV13_M","PV13_M")

# design: Model Matrix by group defined above
design1 <- model.matrix(~0 + group1)
colnames(design1)
rownames(design1) = rownames(x$samples[1:13,])

# removing heteroscedascity from count data: voom plots
v1 <- voom(x, design1,plot = TRUE)
v1
# limma lm fit
fit1 <- lmFit(v1)

# comparison 1: all transcripts that are just B male
# contrast matrix: called fit.cont1
colnames(design1)

# I only want transcripts differentially expressed in male B samples.
cont.matrix1 <- makeContrasts(F04.vs.F13 = group1PV04_F - group1PV13_F, M04.vs.M13 = group1PV04_M - group1PV13_M, levels=design1)
cont.matrix1

fit.cont1 <- contrasts.fit(fit1, cont.matrix1)
fit.cont1 <- eBayes(fit.cont1)
summary(decideTests(fit.cont1))

## Examine the number of DE genes
tfit <- treat(fit.cont1, lfc=0.58)
dt <- decideTests(tfit)
summary(dt)

F04.vs.F13 <- topTreat(tfit, coef=1, n=Inf)
M04.vs.M13 <- topTreat(tfit, coef=2, n=Inf)

colnames(F04.vs.F13)[1] <- "gene"
colnames(M04.vs.M13)[1] <- "gene"
F04.vs.F13$de <- "NS"
M04.vs.M13$de <- "NS"

F04.vs.F13$de <- ifelse(F04.vs.F13$adj.P.Val < 0.05 & F04.vs.F13$logFC < 0, "B-", F04.vs.F13$de)
F04.vs.F13$de <- ifelse(F04.vs.F13$adj.P.Val < 0.05 & F04.vs.F13$logFC > 0, "B+", F04.vs.F13$de)

M04.vs.M13$de <- ifelse(M04.vs.M13$adj.P.Val < 0.05 & M04.vs.M13$logFC < 0, "B-", M04.vs.M13$de)
M04.vs.M13$de <- ifelse(M04.vs.M13$adj.P.Val < 0.05 & M04.vs.M13$logFC > 0, "B+", M04.vs.M13$de)

table(F04.vs.F13$de)
#
#    B-    B+    NS
#   131   217 17196
table(M04.vs.M13$de)
#
#    B-    B+    NS
#   193   120 17231

# merge with annotation
F04.vs.F13.anno <- left_join(F04.vs.F13, genes.by.scaffold, by="gene")
M04.vs.M13.anno <- left_join(M04.vs.M13, genes.by.scaffold, by="gene")
F04.vs.F13.anno <- left_join(F04.vs.F13.anno, freeze.v0.genes.anno, by="gene")
M04.vs.M13.anno <- left_join(M04.vs.M13.anno, freeze.v0.genes.anno, by="gene")
F04.vs.F13.anno.de <- F04.vs.F13.anno[F04.vs.F13.anno$de != "NS",]
M04.vs.M13.anno.de <- M04.vs.M13.anno[M04.vs.M13.anno$de != "NS",]
#write.csv(F04.vs.F13.anno.de[F04.vs.F13.anno.de$anno == "Y",], file="output/B_diff_expr/MB.vs.MnoB.de.annotated.csv") #export results
#write.csv(M04.vs.M13.anno.de[M04.vs.M13.anno.de$anno == "Y",], file="output/B_diff_expr/FB.vs.FnoB.de.annotated.csv") #export results

F04.vs.F13.plot <- F04.vs.F13.anno
F04.vs.F13.plot$de <- factor(F04.vs.F13.plot$de, levels = c("NS", "B-", "B+"))
M04.vs.M13.plot <- M04.vs.M13.anno
M04.vs.M13.plot$de <- factor(M04.vs.M13.plot$de, levels = c("NS", "B-", "B+"))

F04.vs.F13.plot.b1.de <- F04.vs.F13.plot[F04.vs.F13.plot$b.status.final == "B" & F04.vs.F13.plot$de != "NS",]
F04.vs.F13.plot.b1.ns <- F04.vs.F13.plot[F04.vs.F13.plot$b.status.final == "B" & F04.vs.F13.plot$de == "NS",]
F04.vs.F13.plot.Bc.de <- F04.vs.F13.plot[F04.vs.F13.plot$b.status.final == "Bc" & F04.vs.F13.plot$de != "NS",]
F04.vs.F13.plot.Bc.ns <- F04.vs.F13.plot[F04.vs.F13.plot$b.status.final == "Bc" & F04.vs.F13.plot$de == "NS",]
F04.vs.F13.plot.BA.de <- F04.vs.F13.plot[F04.vs.F13.plot$b.status.final == "B-A" & F04.vs.F13.plot$de != "NS",]
F04.vs.F13.plot.BA.ns <- F04.vs.F13.plot[F04.vs.F13.plot$b.status.final == "B-A" & F04.vs.F13.plot$de == "NS",]
F04.vs.F13.plot <- F04.vs.F13.plot[F04.vs.F13.plot$b.status.final != "B",]
F04.vs.F13.plot.de <- F04.vs.F13.plot[F04.vs.F13.plot$de != "NS",]
F04.vs.F13.plot.ns <- F04.vs.F13.plot[F04.vs.F13.plot$de == "NS",]

M04.vs.M13.plot.b1.de <- M04.vs.M13.plot[M04.vs.M13.plot$b.status.final == "B" & M04.vs.M13.plot$de != "NS",]
M04.vs.M13.plot.b1.ns <- M04.vs.M13.plot[M04.vs.M13.plot$b.status.final == "B" & M04.vs.M13.plot$de == "NS",]
M04.vs.M13.plot.Bc.de <- M04.vs.M13.plot[M04.vs.M13.plot$b.status.final == "Bc" & M04.vs.M13.plot$de != "NS",]
M04.vs.M13.plot.Bc.ns <- M04.vs.M13.plot[M04.vs.M13.plot$b.status.final == "Bc" & M04.vs.M13.plot$de == "NS",]
M04.vs.M13.plot.BA.de <- M04.vs.M13.plot[M04.vs.M13.plot$b.status.final == "B-A" & M04.vs.M13.plot$de != "NS",]
M04.vs.M13.plot.BA.ns <- M04.vs.M13.plot[M04.vs.M13.plot$b.status.final == "B-A" & M04.vs.M13.plot$de == "NS",]
M04.vs.M13.plot <- M04.vs.M13.plot[M04.vs.M13.plot$b.status.final != "B",]
M04.vs.M13.plot.de <- M04.vs.M13.plot[M04.vs.M13.plot$de != "NS",]
M04.vs.M13.plot.ns <- M04.vs.M13.plot[M04.vs.M13.plot$de == "NS",]

p1 <- ggplot() +
  geom_vline(xintercept=(0)) +
  geom_point(data=F04.vs.F13.plot.ns, aes(x = logFC, y = AveExpr),fill="#c2b8b8",colour="#c2b8b8", size=1.5,alpha=0.5) +
  geom_point(data=F04.vs.F13.plot.de, aes(x = logFC, y = AveExpr,shape=de),fill="#c2b8b8",colour="gray60") +
  geom_point(data=F04.vs.F13.plot.BA.ns, aes(x = logFC, y = AveExpr),fill="lightpink2",colour="lightpink2") +
  geom_point(data=F04.vs.F13.plot.BA.de, aes(x = logFC, y = AveExpr,shape=de,colour=de),fill="lightpink2",colour="lightpink2") +
	geom_point(data=F04.vs.F13.plot.Bc.ns, aes(x = logFC, y = AveExpr),fill="deepskyblue",colour="deepskyblue") +
	geom_point(data=F04.vs.F13.plot.Bc.de, aes(x = logFC, y = AveExpr,shape=de,colour=de),fill="deepskyblue",colour="deepskyblue") +
	geom_point(data=F04.vs.F13.plot.b1.ns, aes(x = logFC, y = AveExpr),fill="royalblue4",colour="royalblue4") +
	geom_point(data=F04.vs.F13.plot.b1.de, aes(x = logFC, y = AveExpr,shape=de,colour=de),fill="royalblue4",colour="royalblue4") +
  scale_shape_manual(name="Expression",values=c(25,24)) +
  labs(x = "Expression log-ratio (PV04 v PV13)", y = "Average log-expression (TPM)", title = "B+ Females") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (20)),
        axis.title = element_text(family = "Helvetica", size = (15)),
        axis.text = element_text(family = "Helvetica", size = (13)),
        legend.position="none") + guides(scale = 'none')

p2 <- ggplot() +
  geom_vline(xintercept=(0)) +
  geom_point(data=M04.vs.M13.plot.ns, aes(x = logFC, y = AveExpr),fill="#b8bac2",colour="#b8bac2", size=1.5,alpha=0.5) +
  geom_point(data=M04.vs.M13.plot.de, aes(x = logFC, y = AveExpr,shape=de),fill="#b8bac2",colour="gray60") +
	geom_point(data=M04.vs.M13.plot.BA.ns, aes(x = logFC, y = AveExpr),fill="lightpink2",colour="lightpink2") +
	geom_point(data=M04.vs.M13.plot.BA.de, aes(x = logFC, y = AveExpr,shape=de,colour=de),fill="lightpink2",colour="lightpink2") +
  geom_point(data=M04.vs.M13.plot.Bc.ns, aes(x = logFC, y = AveExpr),fill="deepskyblue",colour="deepskyblue") +
  geom_point(data=M04.vs.M13.plot.Bc.de, aes(x = logFC, y = AveExpr,shape=de,colour=de),fill="deepskyblue",colour="deepskyblue") +
  geom_point(data=M04.vs.M13.plot.b1.ns, aes(x = logFC, y = AveExpr),fill="royalblue4",colour="royalblue4") +
  geom_point(data=M04.vs.M13.plot.b1.de, aes(x = logFC, y = AveExpr,shape=de,colour=de),fill="royalblue4",colour="royalblue4") +
  scale_shape_manual(name="Expression",values=c(25,24)) +
  labs(x = "Expression log-ratio (PV04 v PV13)", y = "Average log-expression (TPM)", title = "B+ males") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (20)),
        axis.title = element_text(family = "Helvetica", size = (15)),
        axis.text = element_text(family = "Helvetica", size = (13)),
        legend.position="none") + guides(scale = 'none')

tiff("manuscript/figures_revision/fig6a_expression_PV04_vs_PV13_femlaes.tiff", width = 8, height = 8, units = 'in', res = 300)
    p1
dev.off()

tiff("manuscript/figures_revision/fig6a_expression_PV04_vs_PV13_males.tiff", width = 8, height = 8, units = 'in', res = 300)
    p2
dev.off()

######################### expression barplots

X13F_cpm <- rowMeans(logcounts[, c("X13F_1", "X13F_2", "X13F_3")])
X13M_cpm <- rowMeans(logcounts[, c("X13M_1", "X13M_2", "X13M_3", "X13M_4")])
X04F_cpm <- rowMeans(logcounts[, c("X04F_1","X04F_2","X04F_3")])
X04M_cpm <- rowMeans(logcounts[, c("X04M_1", "X04M_2", "X04M_3")])

cpms <- data.frame(gene = names(X13F_cpm), X13F_cpm = X13F_cpm, X13M_cpm = X13M_cpm, X04F_cpm = X04F_cpm, X04M_cpm = X04M_cpm)
cpms <- left_join(cpms, genes.by.scaffold[, c('gene', 'b.status.final')], by="gene")


unfilt_logcounts_df <- data.frame(gene = rownames(unfilt_logcounts), cmp_means = rowMeans(unfilt_logcounts), cpm_max = apply(unfilt_logcounts, 1, max))
unfilt_logcounts_df <- left_join(unfilt_logcounts_df, genes.by.scaffold[, c('gene', 'b.status.final')], by="gene")

unfilt_logcounts_df[unfilt_logcounts_df$b.status.final == "B" & unfilt_logcounts_df$cpm_max > 0, ]
# g13953, g20085, g2644

# unfilt_logcounts_df[unfilt_logcounts_df$gene %in% c('g13953', 'g1208', 'g5582', 'g9061', 'g9062'), ]

BA_col = "lightpink2"
B_col = "royalblue4"
Bc_col = "deepskyblue"

pdf('output/cpm_expression_PV13_PV04_per_sex.pdf', width = 14, height = 6)

plot(NULL, xlim = c(0, 4), ylim = c(-5, 14), xlab = '', ylab = 'Expression (log2 counts per million)')
boxplot(cpms[cpms$b.status.final == 'A', 'X13F_cpm'], at = 0.15, add = T, boxwex = 0.3)
boxplot(cpms[cpms$b.status.final == 'A', 'X04F_cpm'], at = 0.35, add = T, boxwex = 0.3)
boxplot(cpms[cpms$b.status.final == 'A', 'X13M_cpm'], at = 0.65, add = T, boxwex = 0.3)
boxplot(cpms[cpms$b.status.final == 'A', 'X04M_cpm'], at = 0.85, add = T, boxwex = 0.3)

boxplot(cpms[cpms$b.status.final == 'B-A', 'X13F_cpm'], at = 1.15, add = T, col = BA_col, boxwex = 0.3)
boxplot(cpms[cpms$b.status.final == 'B-A', 'X04F_cpm'], at = 1.35, add = T, col = BA_col, boxwex = 0.3)
boxplot(cpms[cpms$b.status.final == 'B-A', 'X13M_cpm'], at = 1.65, add = T, col = BA_col, boxwex = 0.3)
boxplot(cpms[cpms$b.status.final == 'B-A', 'X04M_cpm'], at = 1.85, add = T, col = BA_col, boxwex = 0.3)

boxplot(cpms[cpms$b.status.final == 'B', 'X13F_cpm'], at = 2.15, add = T, col = B_col, boxwex = 0.3)
boxplot(cpms[cpms$b.status.final == 'B', 'X04F_cpm'], at = 2.35, add = T, col = B_col, boxwex = 0.3)
boxplot(cpms[cpms$b.status.final == 'B', 'X13M_cpm'], at = 2.65, add = T, col = B_col, boxwex = 0.3)
boxplot(cpms[cpms$b.status.final == 'B', 'X04M_cpm'], at = 2.85, add = T, col = B_col, boxwex = 0.3)

boxplot(cpms[cpms$b.status.final == 'Bc', 'X13F_cpm'], at = 3.15, add = T, col = Bc_col, boxwex = 0.3)
boxplot(cpms[cpms$b.status.final == 'Bc', 'X04F_cpm'], at = 3.35, add = T, col = Bc_col, boxwex = 0.3)
boxplot(cpms[cpms$b.status.final == 'Bc', 'X13M_cpm'], at = 3.65, add = T, col = Bc_col, boxwex = 0.3)
boxplot(cpms[cpms$b.status.final == 'Bc', 'X04M_cpm'], at = 3.85, add = T, col = Bc_col, boxwex = 0.3)
#
dev.off()
