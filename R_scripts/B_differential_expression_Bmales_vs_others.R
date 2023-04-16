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

replicate <- as.factor(c("X04F","X04F","X04F","X04M","X04M","X04M","X13F","X13F","X13F","X13M","X13M","X13M","X13M","X15F","X15F","X15F","X15M","X15M","X15M","X21F","X21F","X21F","X21M","X21M","X21M","X21M"))
x <- DGEList(counts=round(a), genes=rownames(a), group = replicate)
head(x)

unfilt_logcounts <- cpm(x, log=TRUE)

# rowSums(x$counts==0) # for each sample where the count is 0, then sums the samples that have counts = 0 for a gene
zero_count <- sum(rowSums(x$counts==0) == length(replicate)) #number of genes that have all samples with 0 counts
proportion_0_count <- zero_count * 100 / nrow(x)


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

# define groups by being B-males or not
group1 <- rep("noBM", length(replicate))
group1[replicate %in% c("X04M", "X13M")] <- "BM"

# design: Model Matrix by group defined above
design1 <- model.matrix(~0 + group1)
colnames(design1)
rownames(design1) = rownames(x$samples)

# removing heteroscedascity from count data: voom plots
v1 <- voom(x, design1,plot = TRUE)
v1
# limma lm fit
fit1 <- lmFit(v1)

# comparison 1: all transcripts that are just B male
# contrast matrix: called fit.cont1
colnames(design1)

# I only want transcripts differentially expressed in male B samples.
cont.matrix1 <- makeContrasts(BM.vs.all = group1BM - group1noBM, levels=design1)
cont.matrix1

fit.cont1 <- contrasts.fit(fit1, cont.matrix1)
fit.cont1 <- eBayes(fit.cont1)
summary(decideTests(fit.cont1))

## Examine the number of DE genes
tfit <- treat(fit.cont1, lfc=0.58)
dt <- decideTests(tfit)
summary(dt)

BM.vs.all <- topTreat(tfit, coef=1, n=Inf)

colnames(BM.vs.all)[1] <- "gene"
BM.vs.all$de <- "NS"

BM.vs.all$de <- ifelse(BM.vs.all$adj.P.Val < 0.05 & BM.vs.all$logFC < 0, "BM down", BM.vs.all$de)
BM.vs.all$de <- ifelse(BM.vs.all$adj.P.Val < 0.05 & BM.vs.all$logFC > 0, "BM up", BM.vs.all$de)


table(BM.vs.all$de)
#
#BM down   BM up      NS
#     28     875   17163

# merge with annotation
BM.vs.all.anno <- left_join(BM.vs.all, genes.by.scaffold, by="gene")
BM.vs.all.anno <- left_join(BM.vs.all.anno, freeze.v0.genes.anno, by="gene")
BM.vs.all.anno.de <- BM.vs.all.anno[BM.vs.all.anno$de != "NS",]
#write.csv(BM.vs.all.anno.de[BM.vs.all.anno.de$anno == "Y",], file="output/B_diff_expr/MB.vs.MnoB.de.annotated.csv") #export results
#write.csv(M04.vs.M13.anno.de[M04.vs.M13.anno.de$anno == "Y",], file="output/B_diff_expr/FB.vs.FnoB.de.annotated.csv") #export results

BM.vs.all.plot <- BM.vs.all.anno
BM.vs.all.plot$de <- factor(BM.vs.all.plot$de)

BM.vs.all.plot.b1.de <- BM.vs.all.plot[BM.vs.all.plot$b.status.final == "B" & BM.vs.all.plot$de != "NS",]
BM.vs.all.plot.b1.ns <- BM.vs.all.plot[BM.vs.all.plot$b.status.final == "B" & BM.vs.all.plot$de == "NS",]
BM.vs.all.plot.Bc.de <- BM.vs.all.plot[BM.vs.all.plot$b.status.final == "Bc" & BM.vs.all.plot$de != "NS",]
BM.vs.all.plot.Bc.ns <- BM.vs.all.plot[BM.vs.all.plot$b.status.final == "Bc" & BM.vs.all.plot$de == "NS",]
BM.vs.all.plot.BA.de <- BM.vs.all.plot[BM.vs.all.plot$b.status.final == "B-A" & BM.vs.all.plot$de != "NS",]
BM.vs.all.plot.BA.ns <- BM.vs.all.plot[BM.vs.all.plot$b.status.final == "B-A" & BM.vs.all.plot$de == "NS",]
BM.vs.all.plot <- BM.vs.all.plot[BM.vs.all.plot$b.status.final != "B",]
BM.vs.all.plot.de <- BM.vs.all.plot[BM.vs.all.plot$de != "NS",]
BM.vs.all.plot.ns <- BM.vs.all.plot[BM.vs.all.plot$de == "NS",]

p1 <- ggplot() +
  geom_vline(xintercept=(0)) +
  geom_point(data=BM.vs.all.plot.ns, aes(x = logFC, y = AveExpr),fill="#c2b8b8",colour="#c2b8b8", size=1.5,alpha=0.5) +
  geom_point(data=BM.vs.all.plot.de, aes(x = logFC, y = AveExpr,shape=de),fill="#c2b8b8",colour="gray60") +
  geom_point(data=BM.vs.all.plot.BA.ns, aes(x = logFC, y = AveExpr),fill="lightpink2",colour="lightpink2") +
  geom_point(data=BM.vs.all.plot.BA.de, aes(x = logFC, y = AveExpr,shape=de,colour=de),fill="lightpink2",colour="lightpink2") +
	geom_point(data=BM.vs.all.plot.Bc.ns, aes(x = logFC, y = AveExpr),fill="deepskyblue",colour="deepskyblue") +
	geom_point(data=BM.vs.all.plot.Bc.de, aes(x = logFC, y = AveExpr,shape=de,colour=de),fill="deepskyblue",colour="deepskyblue") +
	geom_point(data=BM.vs.all.plot.b1.ns, aes(x = logFC, y = AveExpr),fill="royalblue4",colour="royalblue4") +
	geom_point(data=BM.vs.all.plot.b1.de, aes(x = logFC, y = AveExpr,shape=de,colour=de),fill="royalblue4",colour="royalblue4") +
  scale_shape_manual(name="Expression",values=c(25,24)) +
  labs(x = "Expression log-ratio (B males v all)", y = "Average log-expression (TPM)", title = "B+ Males vs all") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (20)),
        axis.title = element_text(family = "Helvetica", size = (15)),
        axis.text = element_text(family = "Helvetica", size = (13)),
        legend.position="none") + guides(scale = 'none')


tiff("output/B+_males_vs_all.tiff", width = 8, height = 8, units = 'in', res = 300)
    p1
dev.off()

BM.vs.all.anno <- BM.vs.all.anno[order(BM.vs.all.anno$logFC, decreasing = T),]
BM_upregulated_Autosomal <- BM.vs.all.anno[BM.vs.all.anno$anno == 'Y' & BM.vs.all.anno$de == 'BM up' & BM.vs.all.anno$b.status.final == 'A', ]
BM_upregulated_Blinked <- BM.vs.all.anno[BM.vs.all.anno$anno == 'Y' & BM.vs.all.anno$de == 'BM up' & BM.vs.all.anno$b.status.final != 'A', ]
BM_downregulated <- BM.vs.all.anno[BM.vs.all.anno$anno == 'Y' & BM.vs.all.anno$de == 'BM down', ]

write.table(BM_upregulated_Autosomal, 'output/BM_upregulated_Autosomal.tsv', sep = '\t', quote = F, row.names = F, col.names = T)
write.table(BM_upregulated_Blinked, 'output/BM_upregulated_Blinked.tsv', sep = '\t', quote = F, row.names = F, col.names = T)
write.table(BM_downregulated, 'output/BM_downregulated.tsv', sep = '\t', quote = F, row.names = F, col.names = T)
