#!/usr/bin/env Rscript

library(edgeR)
library(dplyr)
library(tidyverse)
library(readr)

B_col = "royalblue4"
B_c_col = "deepskyblue"

# genetotranscript <- read.delim("annotation/p.viburni.freeze.v0.braker.transcripts.to.genes.txt", header=FALSE) #mapping genes to transcripts
# genetotranscript <- genetotranscript[order(genetotranscript$V1),]
rsem.counts <- read.delim("R_scripts/RSEM_digi.counts.matrix", header=TRUE) #matrix generated from rsem
sampleinfo <- read.csv("R_scripts/sampleinfoPviburniB.csv") #sample group info
freeze.v0.genes.anno <- read_delim("output/freeze.v0.genes.anno.complete.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE) # master anno
# # scaffolds.preprint.assignment <- read_delim("output/scaffolds.preprint.assignment.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE)
# scaffolds.final.assignment <- read_delim('output/scaffolds.final.assignment.tsv',"\t", escape_double = FALSE, col_names = T,trim_ws = TRUE)

# Assignments
genes.by.scaffold <- read.table("output/genes.by.scaffolds.tsv", sep = "\t", header = T)
# genes.by.scaffold$b.status.final <- factor(genes.by.scaffold$b.status.final, levels = c('A', 'B-A', 'B', 'Bc'))
# gene           seq b.status.final
# 1     g1 scaffold_1617              A
# 2    g10   scaffold_93              A
# 3   g100  scaffold_360            B-A


# count file from all samples
# need a dataframe containing all gene info per gene id
colnames(rsem.counts) <- substr(colnames(rsem.counts),start=1,stop=6) #removing ".genes.results" in colnames

#first column X is the gene id, this needs to be removed but has the order of gene id has to match when merging with gene info
gene_names <- rsem.counts[, 1]
rsem.counts <- rsem.counts[, -1]
rownames(rsem.counts) <- gene_names

#adding replicate information
replicate <- as.factor(substr(colnames(rsem.counts), 1, 4))

#create the DGE object including all count matrix, replicate info and annotations when available
dge_list <- DGEList(counts=round(rsem.counts), genes=gene_names, group=replicate)

# find all genes with NO reads mapped from any of the replicates
table(rowSums(dge_list$counts==0) == 26) #number of genes that have all samples with 0 counts
proportion_0_count <- 2393 * 100 / nrow(dge_list)

# compare filtering options. By default, the function keeps genes with about 10 read counts or more in a minimum number of samples, where the number of samples is chosen according to the minimum group sample size. The actual filtering uses CPM values rather than counts in order to avoid giving preference to samples with large library sizes.
# reduce min.count to 5 to account for the low expression of B genes
keep.exprs.group <- filterByExpr(dge_list, group=dge_list$samples$group, min.count = 5)
dge_list <- dge_list[keep.exprs.group, keep.lib.sizes=FALSE]

# normalize distribution (TMM normalization)
dge_list <- calcNormFactors(dge_list, method = "TMM")

## for supplementary material

tiff('manuscript/figures_revision/supplmentary_figure_S4_expression.tiff', width = 8, height = 8, units = 'in', res = 150)
	par(mfrow=c(2,2))
	# QC
	dge_list$samples$lib.size
	barplot(dge_list$samples$lib.size,names=colnames(dge_list),las=2)
	title("Barplot of library sizes")
	logcounts <- cpm(dge_list, log=TRUE)
	boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
	abline(h=median(logcounts),col="blue")
	title("Boxplots of logCPMs (normalized)")

	# MDS plot to check for grouping
	# checking grouping by sex, and by B presence

	sampleinfo2 <- sampleinfo %>% mutate(Sex = as.factor(Sex))
	levels(sampleinfo2$Sex)
	col.sex <- c("red","blue")[sampleinfo2$Sex]
	data.frame(sampleinfo2$Sex,col.sex)

	# plot by sex
	plotMDS(dge_list, col=col.sex, main="Sex", cex=0.8)
	legend("top",fill=c("red","blue"),legend=levels(sampleinfo2$Sex))

	# plot by B presence
	sampleinfo2 <- sampleinfo %>% mutate(Bpresence = as.factor(Bpresence))
	levels(sampleinfo2$Bpresence)
	col.B <- c("orange","purple")[sampleinfo2$Bpresence]
	data.frame(sampleinfo2$Bpresence, col.B)
	plotMDS(dge_list, col=col.B, main="B presence",cex=0.8)
	legend("top", fill=c("orange","purple"), legend=levels(sampleinfo2$Bpresence))
dev.off()

## Hierarchical clustering
var_genes <- apply(logcounts, 1, var)
select_var <- names(sort(var_genes, decreasing=TRUE))[1:1000]
highly_variable_lcpm <- logcounts[select_var,]

## Model design

# The model design and fit is based on sex by B presence groups. Model will be fit based on this group, then we can compare pairs of groups by creating different contrast matrices.

# define groups by sex and B presence or absence (corresponding lines in comments)
group1 <- c("FB","FB","FB", # 04
            "MB","MB","MB",
            "FB","FB","FB", # 13
            "MB","MB","MB","MB",
            "FnoB","FnoB","FnoB", # 15
            "MnoB","MnoB","MnoB",
            "FnoB","FnoB","FnoB", # 21
            "MnoB","MnoB","MnoB","MnoB")

# design: Model Matrix by group defined above
design1 <- model.matrix(~0 + group1)
colnames(design1)
rownames(design1) <- rownames(dge_list$samples)

# removing heteroscedascity from count data: voom plots
v1 <- voom(dge_list, design1, plot = F)

# limma lm fit
fit1 <- lmFit(v1)


# I only want transcripts differentially expressed in male B samples.
cont.matrix1 <- makeContrasts(MB.vs.MnoB = group1MB - group1MnoB, MB.vs.FB = group1MB - group1FB, MB.vs.FnoB = group1MB - group1FnoB, FB.vs.FnoB = group1FB - group1FnoB, levels=design1)

fit.cont1 <- contrasts.fit(fit1, cont.matrix1)
fit.cont1 <- eBayes(fit.cont1)

# compare mean-variance trend
v1 <- voom(dge_list, design1, plot = F)

## Examine the number of DE genes
tfit <- treat(fit.cont1, lfc=0.58)
dt <- decideTests(tfit)
summary(dt)
#write.fit(tfit, dt, file="output/B_diff_expr/results.txt")

# these are all the genes that are over expressed and underexpressed, respectively, between B males and all the others
de.over.B.males <- which(dt[,1]==1 & dt[,2]==1 &dt[,3]==1)
de.under.B.males <- which(dt[,1]==-1 & dt[,2]==-1 &dt[,3]==-1)
length(de.over.B.males) # B genes overexpressed in B males compared to the rest of the groups
length(de.under.B.males) # B genes underexpressed in B males compared to the rest of the groups

# examine DE genes in M and F
MB.vs.MnoB.results <- topTreat(tfit, coef=1, n=Inf)
FB.vs.FnoB.results <- topTreat(tfit, coef=4, n=Inf)

m.tfit <- MB.vs.MnoB.results
f.tfit <- FB.vs.FnoB.results

colnames(m.tfit)[1] <- "gene"
colnames(f.tfit)[1] <- "gene"
m.tfit$de <- "NS"
f.tfit$de <- "NS"

m.tfit$de <- ifelse(m.tfit$adj.P.Val < 0.05 & m.tfit$logFC < 0, "B-", m.tfit$de)
m.tfit$de <- ifelse(m.tfit$adj.P.Val < 0.05 & m.tfit$logFC > 0, "B+", m.tfit$de)

f.tfit$de <- ifelse(f.tfit$adj.P.Val < 0.05 & f.tfit$logFC < 0, "B-", f.tfit$de)
f.tfit$de <- ifelse(f.tfit$adj.P.Val < 0.05 & f.tfit$logFC > 0, "B+", f.tfit$de)

table(m.tfit$de)
table(f.tfit$de)

# merge with annotation
m.tfit.anno <- left_join(m.tfit, genes.by.scaffold, by="gene")
f.tfit.anno <- left_join(f.tfit, genes.by.scaffold, by="gene")
m.tfit.anno <- left_join(m.tfit.anno, freeze.v0.genes.anno, by="gene")
f.tfit.anno <- left_join(f.tfit.anno, freeze.v0.genes.anno, by="gene")
m.tfit.anno.de <- m.tfit.anno[m.tfit.anno$de != "NS",]
f.tfit.anno.de <- f.tfit.anno[f.tfit.anno$de != "NS",]
#write.csv(m.tfit.anno.de[m.tfit.anno.de$anno == "Y",], file="output/B_diff_expr/MB.vs.MnoB.de.annotated.csv") #export results
#write.csv(f.tfit.anno.de[f.tfit.anno.de$anno == "Y",], file="output/B_diff_expr/FB.vs.FnoB.de.annotated.csv") #export results

m.tfit.plot <- m.tfit.anno
m.tfit.plot$de <- factor(m.tfit.plot$de, levels = c("NS", "B-", "B+"))
f.tfit.plot <- f.tfit.anno
f.tfit.plot$de <- factor(f.tfit.plot$de, levels = c("NS", "B-", "B+"))

m.tfit.plot.b1.de <- m.tfit.plot[m.tfit.plot$b.status.final == "B" & m.tfit.plot$de != "NS",]
m.tfit.plot.b1.ns <- m.tfit.plot[m.tfit.plot$b.status.final == "B" & m.tfit.plot$de == "NS",]
m.tfit.plot.Bc.de <- m.tfit.plot[m.tfit.plot$b.status.final == "Bc" & m.tfit.plot$de != "NS",]
m.tfit.plot.Bc.ns <- m.tfit.plot[m.tfit.plot$b.status.final == "Bc" & m.tfit.plot$de == "NS",]
m.tfit.plot.BA.de <- m.tfit.plot[m.tfit.plot$b.status.final == "B-A" & m.tfit.plot$de != "NS",]
m.tfit.plot.BA.ns <- m.tfit.plot[m.tfit.plot$b.status.final == "B-A" & m.tfit.plot$de == "NS",]
m.tfit.plot <- m.tfit.plot[m.tfit.plot$b.status.final != "B",]
m.tfit.plot.de <- m.tfit.plot[m.tfit.plot$de != "NS",]
m.tfit.plot.ns <- m.tfit.plot[m.tfit.plot$de == "NS",]

f.tfit.plot.b1.de <- f.tfit.plot[f.tfit.plot$b.status.final == "B" & f.tfit.plot$de != "NS",]
f.tfit.plot.b1.ns <- f.tfit.plot[f.tfit.plot$b.status.final == "B" & f.tfit.plot$de == "NS",]
f.tfit.plot.BA.de <- f.tfit.plot[f.tfit.plot$b.status.final == "B-A" & f.tfit.plot$de != "NS",]
f.tfit.plot.BA.ns <- f.tfit.plot[f.tfit.plot$b.status.final == "B-A" & f.tfit.plot$de == "NS",]
f.tfit.plot <- f.tfit.plot[f.tfit.plot$b.status.final != "B",]
f.tfit.plot.de <- f.tfit.plot[f.tfit.plot$de != "NS",]
f.tfit.plot.ns <- f.tfit.plot[f.tfit.plot$de == "NS",]

p1 <- ggplot() +
  geom_vline(xintercept=(0)) +
  geom_point(data=m.tfit.plot.ns, aes(x = logFC, y = AveExpr),fill="#b8bac2",colour="#b8bac2", size=1.5,alpha=0.5) +
  geom_point(data=m.tfit.plot.de, aes(x = logFC, y = AveExpr,shape=de),fill="#b8bac2",colour="gray60") +
  geom_point(data=m.tfit.plot.BA.ns, aes(x = logFC, y = AveExpr),fill="lightpink2",colour="lightpink2") +
  geom_point(data=m.tfit.plot.BA.de, aes(x = logFC, y = AveExpr,shape=de,colour=de),fill="lightpink2",colour="lightpink2") +
	# geom_point(data=m.tfit.plot.Bc.ns, aes(x = logFC, y = AveExpr),fill="deepskyblue",colour="deepskyblue") +
	# geom_point(data=m.tfit.plot.Bc.de, aes(x = logFC, y = AveExpr,shape=de,colour=de),fill="deepskyblue",colour="deepskyblue") +
	geom_point(data=m.tfit.plot.b1.ns, aes(x = logFC, y = AveExpr),fill="royalblue4",colour="royalblue4") +
	geom_point(data=m.tfit.plot.b1.de, aes(x = logFC, y = AveExpr,shape=de,colour=de),fill="royalblue4",colour="royalblue4") +
  scale_shape_manual(name="Expression",values=c(25,24)) +
  labs(x = "Expression log-ratio (B+ v B-)", y = "Average log-expression (TPM)", title = "Males") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (20)),
        axis.title = element_text(family = "Helvetica", size = (15)),
        axis.text = element_text(family = "Helvetica", size = (13)),
        legend.position="none") + guides(scale = 'none')

p2 <- ggplot() +
  geom_vline(xintercept=(0)) +
  geom_point(data=f.tfit.plot.ns, aes(x = logFC, y = AveExpr),fill="#c2b8b8",colour="#c2b8b8", size=1.5,alpha=0.5) +
  geom_point(data=f.tfit.plot.de, aes(x = logFC, y = AveExpr,shape=de),fill="#c2b8b8",colour="gray60") +
	geom_point(data=f.tfit.plot.BA.ns, aes(x = logFC, y = AveExpr),fill="lightpink2",colour="lightpink2") +
	geom_point(data=f.tfit.plot.BA.de, aes(x = logFC, y = AveExpr,shape=de,colour=de),fill="lightpink2",colour="lightpink2") +
  geom_point(data=f.tfit.plot.b1.ns, aes(x = logFC, y = AveExpr),fill="royalblue4",colour="royalblue4") +
  geom_point(data=f.tfit.plot.b1.de, aes(x = logFC, y = AveExpr,shape=de,colour=de),fill="royalblue4",colour="royalblue4") +
  scale_shape_manual(name="Expression",values=c(25,24)) +
  labs(x = "Expression log-ratio (B+ v B-)", y = "Average log-expression (TPM)", title = "Females") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (20)),
        axis.title = element_text(family = "Helvetica", size = (15)),
        axis.text = element_text(family = "Helvetica", size = (13)),
        legend.position="none") + guides(scale = 'none')

# get gene counts

with(m.tfit.anno.de, table(b.status.final, anno, de))
with(f.tfit.anno.de, table(b.status.final, anno, de))

# extracting and examined DE that are differentially expressed in B males compared to the other groups

over.bm <- which(dt[,1]==1 & dt[,2]==1 &dt[,3]==1)
under.bm <- which(dt[,1]==-1 & dt[,2]==-1 &dt[,3]==-1)
length(under.bm)
# list of genes
fit.cont1$genes[over.bm,]
fit.cont1$genes[under.bm,]

# export overexpressed genes in B+ males vs all
over.bm.genes <- data.frame(fit.cont1$genes[over.bm,])
colnames(over.bm.genes) <- "gene"
over.bm.genes <- left_join(over.bm.genes, genes.by.scaffold, by="gene")
over.bm.genes <- left_join(over.bm.genes, freeze.v0.genes.anno, by="gene")
over.bm.genes

# export underexpressed genes in B+ males vs all
under.bm.genes <- data.frame(fit.cont1$genes[under.bm,])
colnames(under.bm.genes) <- "gene"
under.bm.genes <- left_join(under.bm.genes, genes.by.scaffold, by="gene")
under.bm.genes <- left_join(under.bm.genes, freeze.v0.genes.anno, by="gene")
under.bm.genes

#write.csv(over.bm.genes,file="output/B_diff_expr/over.Bmales.vs.all.csv")
#write.csv(under.bm.genes,file="output/B_diff_expr/under.Bmales.vs.all.csv")
with(under.bm.genes, table(b.status.final,anno))

counts.a <- data.frame(with(over.bm.genes, table(b.status.final))) # table of overexpressed genes vs B status
with(over.bm.genes, table(b.status.final[anno == 'Y']))  # table of overexpressed genes with annotation vs B status
counts.b <- data.frame(with(under.bm.genes, table(b.status.final))) # the smae for underexpressed
with(under.bm.genes, table(b.status.final[anno == 'Y']))
counts.a$dir <- "O"
counts.b$dir <- "U"
counts <- rbind(counts.a,counts.b)
# data.frame()

Status <- as.factor(c('A','B','Bc','A','B','Bc'))
# Status <- as.factor(c('A','B1','B2/B3','A','B1','B2/B3'))
#   b.status.final Freq dir
# 1              A   85   O
# 2              B    1   O
# 3             Bc    2   O
# 4              A    9   U
Total <- as.integer(c(85,1,2,9,0,0))
Anno <- as.factor(c("85 (39)","1 (1)","2 (1)","9 (6)","0 (0)","0 (0)"))
Dir <- as.factor(c("Overexpressed","Overexpressed","Overexpressed","Underexpressed","Underexpressed","Underexpressed"))
counts <- data.frame(Status, Total, Anno, Dir)
counts

p4 <- ggplot(counts,aes(Status, Total, fill=Dir)) + ylim(0,90) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(label=Anno),position = position_dodge(0.9),vjust=-1) +
  scale_fill_manual(name="",values=c("gray20","gray60")) +
  labs(x = "Location on scaffold", y = "Number of genes", title = "DE in B+ males") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (20)),
        axis.title = element_text(family = "Helvetica", size = (15)),
        axis.text = element_text(family = "Helvetica", size = (13)),
        legend.position=c(0.75,0.8), legend.title = element_blank(), legend.text = element_text(family = "Helvetica", size = (15)))

tiff("manuscript/figures_revision/fig6a_expression.tiff", width = 8, height = 8, units = 'in', res = 300)
    p1
dev.off()

tiff("manuscript/figures_revision/fig6b_expression.tiff", width = 8, height = 8, units = 'in', res = 300)
    p2
dev.off()

tiff("manuscript/figures_revision/fig4c_expression.tiff", width = 8, height = 8, units = 'in', res = 300)
    vennDiagram(dt[,1:3], circle.col=c("purple", "green","orange"),include=c("up","down"), mar = c(0,0,0,0))
dev.off()

tiff("manuscript/figures_revision/fig6d_expression.tiff", width = 8, height = 8, units = 'in', res = 300)
    p4
dev.off()
