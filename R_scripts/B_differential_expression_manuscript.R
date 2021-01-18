rm(list=ls())
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

# get files

setwd("E:/agdel/Documents/projects_rosslab/B_viburni")

genetotranscript<-read.delim("annotation/p.viburni.freeze.v0.braker.transcripts.to.genes.txt", header=FALSE) #mapping genes to transcripts
genetotranscript <-genetotranscript[order(genetotranscript$V1),]

sampleinfo<-read.csv("R_scripts/sampleinfoPviburniB.csv") #sample group info
annotation<-read.csv("output/freeze.v0.genes.anno.complete.csv") # annotation

# count file from all samples
# need a dataframe containing all gene info per gene id
rsem.counts <-read.delim("E:/agdel/Documents/projects_rosslab/B_viburni/R_scripts/RSEM_digi.counts.matrix",header=TRUE) #matrix generated from rsem
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
a<-a1
head(a)

#adding replicate information
replicate = as.factor(c("X04F","X04F","X04F","X04M","X04M","X04M","X13F","X13F","X13F","X13M","X13M","X13M","X13M","X15F","X15F","X15F","X15M","X15M","X15M","X21F","X21F","X21F","X21M","X21M","X21M","X21M"))

#create the DGE object including all count matrix, replicate info and annotations when available
x<-DGEList(counts=round(a),genes=rownames(a), group = replicate)
head(x)

# filter out low count

rowSums(x$counts==0) # for each sample where the count is 0, then sums the samples that have counts = 0 for a gene
table(rowSums(x$counts==0) == 26) #number of genes that have all samples with 0 counts
proportion_0_count= 2393*100/nrow(x)
proportion_0_count

dim(x)

# compare filtering options. By default, the function keeps genes with about 10 read counts or more in a minimum number of samples, where the number of samples is chosen according to the minimum group sample size. The actual filtering uses CPM values rather than counts in order to avoid giving preference to samples with large library sizes.
# reduce min.count to 5 to account for the low expression of B genes
keep.exprs.group <- filterByExpr(x, group=x$samples$group,min.count=5)
keep.exprs.group[keep.exprs.group == FALSE]
x1 <- x[keep.exprs.group, keep.lib.sizes=FALSE]
dim(x1)
x <- x1

# plot  individual samples pre normalization
for (i in 1:ncol(x)) {
	plotMD(cpm(x, log=TRUE), column=i, xlab = "Average log-expression", ylab = "Expression log-ratio (this sample vs others)", main = colnames(x)[i])
	abline(h=0, col="red", lty=2, lwd=2)
}

# normalize distribution (TMM normalization)
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

# plot  individual samples post normalization
for (i in 1:ncol(x)) {
	plotMD(cpm(x, log=TRUE), column=i, xlab = "Average log-expression", ylab = "Expression log-ratio (this sample vs others)", main = colnames(x)[i])
	abline(h=0, col="red", lty=2, lwd=2)
}


## for supplementary material

par(mfrow=c(2,2))
# QC
x$samples$lib.size
barplot(x$samples$lib.size,names=colnames(x),las=2)
title("Barplot of library sizes")
logcounts <- cpm(x,log=TRUE)
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
plotMDS(x,col=col.sex, main="Sex",cex=0.8)
legend("top",fill=c("red","blue"),legend=levels(sampleinfo2$Sex))

# plot by B presence
sampleinfo2 <- sampleinfo %>% mutate(Bpresence = as.factor(Bpresence))
levels(sampleinfo2$Bpresence)
col.B <- c("orange","purple")[sampleinfo2$Bpresence]
data.frame(sampleinfo2$Bpresence,col.B)
plotMDS(x,col=col.B, main="B presence",cex=0.8)
legend("top",fill=c("orange","purple"),legend=levels(sampleinfo2$Bpresence))

## Hierarchical clustering
var_genes <- apply(logcounts, 1, var)
head(var_genes)
select_var <- names(sort(var_genes, decreasing=TRUE))[1:1000]
head(select_var)
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)

# plot the heatmap by sex
library(RColorBrewer)
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 1000 most variable genes across samples",ColSideColors=col.sex,scale="row")

# plot the heatmap by B presence
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 1000 most variable genes across samples",ColSideColors=col.B,scale="row")

## Model design

# The model design and fit is based on sex by B presence groups. Model will be fit based on this group, then we can compare pairs of groups by creating different contrast matrices.

# define groups by sex and B presence or absence
group1=c("FB","FB","FB","MB","MB","MB","FB","FB","FB","MB","MB","MB","MB","FnoB","FnoB","FnoB","MnoB","MnoB","MnoB","FnoB","FnoB","FnoB","MnoB","MnoB","MnoB","MnoB")

# design: Model Matrix by group defined above
design1 <- model.matrix(~0 + group1)
colnames(design1)
rownames(design1) = rownames(x$samples)

# removing heteroscedascity from count data: voom plots
v1 <- voom(x,design1,plot = TRUE)
v1
# limma lm fit
fit1 <- lmFit(v1)

# comparison 1: all transcripts that are just B male
# contrast matrix: called fit.cont1
colnames(design1)

# I only want transcripts differentially expressed in male B samples.
cont.matrix1 <- makeContrasts(MB.vs.MnoB = group1MB - group1MnoB, MB.vs.FB = group1MB - group1FB, MB.vs.FnoB = group1MB - group1FnoB, FB.vs.FnoB = group1FB - group1FnoB, levels=design1)
cont.matrix1

fit.cont1 <- contrasts.fit(fit1, cont.matrix1)
fit.cont1 <- eBayes(fit.cont1)
summary(decideTests(fit.cont1))

# compare mean-variance trend
v1 <- voom(x,design1,plot = TRUE)
plotSA(fit.cont1, main="Final model: Mean-variance trend")

## Examine the number of DE genes
tfit <- treat(fit.cont1, lfc=0.58) 
dt <- decideTests(tfit)
summary(dt)
#write.fit(tfit, dt, file="output/B_diff_expr/results.txt")

par(mfrow=c(2,2))
#Venn Diagram for B vs no B in males
vennDiagram(dt[,1], circle.col=c("turquoise", "salmon","orange"),include=c("up","down"))
#Venn Diagram for B vs no B in males
vennDiagram(dt[,2], circle.col=c("turquoise", "salmon","orange"),include=c("up","down"))
#Venn Diagram for B vs no B in males
vennDiagram(dt[,3], circle.col=c("turquoise", "salmon","orange"),include=c("up","down"))
#Venn Diagram for B vs no B in females
vennDiagram(dt[,4], circle.col=c("turquoise", "salmon","orange"),include=c("up","down"))

#Venn Diagram for Bmale
par(mfrow=c(1,1))
vennDiagram(dt[,1:3], circle.col=c("purple", "green","orange"),include=c("up","down"))

# figs
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

# extracting only the B male specific

# these are all the genes that are DE between B males and all the others
de.B.vs.all <- which(dt[,1]!=0 & dt[,2]!=0 &dt[,3]!=0)
length(de.B.vs.all)

# these are all the genes that are over expressed and underexpressed, respectively, between B males and all the others
de.over.B.males <- which(dt[,1]==1 & dt[,2]==1 &dt[,3]==1)
de.under.B.males <- which(dt[,1]==-1 & dt[,2]==-1 &dt[,3]==-1)
length(de.over.B.males) # B genes overexpressed in B males compared to the rest of the groups
length(de.under.B.males) # B genes underexpressed in B males compared to the rest of the groups

# list of genes from Bmale
fit.cont1$genes[de.over.B.males,]
fit.cont1$genes[de.under.B.males,]

# export overexpressed genes in B+ males vs all
de.over.B.males.genes <- data.frame(fit.cont1$genes[de.over.B.males,])
colnames(de.over.B.males.genes) <- "genes"
#write.csv(de.over.B.males.genes,"output/diff_expr/over.Bmales.vs.all.csv")
head(de.over.B.males.genes)

# export all lists of genes 
maleB.noB <- topTreat(tfit, coef=1,number = summary(dt)[1]+summary(dt)[3])
maleB.femaleB <- topTreat(tfit, coef=2, number = summary(dt)[4]+summary(dt)[6]) 
maleB.femalenoB <- topTreat(tfit, coef=3, number = summary(dt)[7]+summary(dt)[9] ) 
femaleB.femalenoB <- topTreat(tfit, coef=4, number = summary(dt)[10]+summary(dt)[12]) 
malenoB.femalenoB <- topTreat(tfit, coef=5, number = summary(dt)[13]+summary(dt)[15]) 

#write.csv(maleB.noB,"output/diff_expr/B.males.vs.nonB.males.de.treat.csv")
#write.csv(maleB.femaleB,"output/diff_expr/B.males.vs.B.females.de.treat.csv")
#write.csv(maleB.femalenoB,"output/diff_expr/B.males.vs.nonB.females.de.treat.csv")
#write.csv(femaleB.femalenoB,"output/diff_expr/B.females.vs.nonB.females.de.treat.csv")
#write.csv(malenoB.femalenoB,"output/diff_expr/nonB.males.vs.nonB.females.de.treat.csv")

#additional exports: complete list of contrasts

dt_df <- as.data.frame(dt)
dt_df$gene <- row.names(dt_df)
#write.csv(dt_df,"output/diff_expr/dt_df.csv")
