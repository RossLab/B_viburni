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
rsem.counts <-read.delim("R_scripts/RSEM_digi.counts.matrix",header=TRUE) #matrix generated from rsem
sampleinfo<-read.csv("R_scripts/sampleinfoPviburniB.csv") #sample group info
freeze.v0.genes.anno <- read_delim("output/freeze.v0.genes.anno.complete.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE) # master anno
scaffolds.final.assignment <- read_delim("output/scaffolds.final.assignment.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE)
genes.by.scaffold <- read_delim("output/genes.by.scaffolds.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE)


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

m.tfit.plot.b1.de <- m.tfit.plot[m.tfit.plot$b.status.final == "B1" & m.tfit.plot$de != "NS",]
m.tfit.plot.b1.ns <- m.tfit.plot[m.tfit.plot$b.status.final == "B1" & m.tfit.plot$de == "NS",]
m.tfit.plot <- m.tfit.plot[m.tfit.plot$b.status.final != "B1",]
m.tfit.plot.de <- m.tfit.plot[m.tfit.plot$de != "NS",]
m.tfit.plot.ns <- m.tfit.plot[m.tfit.plot$de == "NS",]

f.tfit.plot.b1.de <- f.tfit.plot[f.tfit.plot$b.status.final == "B1" & f.tfit.plot$de != "NS",]
f.tfit.plot.b1.ns <- f.tfit.plot[f.tfit.plot$b.status.final == "B1" & f.tfit.plot$de == "NS",]
f.tfit.plot <- f.tfit.plot[f.tfit.plot$b.status.final != "B1",]
f.tfit.plot.de <- f.tfit.plot[f.tfit.plot$de != "NS",]
f.tfit.plot.ns <- f.tfit.plot[f.tfit.plot$de == "NS",]

p1 <- ggplot() +
  geom_vline(xintercept=(0)) +
  geom_point(data=m.tfit.plot.ns, aes(x = logFC, y = AveExpr),fill="#b8bac2",colour="#b8bac2", size=1.5,alpha=0.5) +
  geom_point(data=m.tfit.plot.de, aes(x = logFC, y = AveExpr,shape=de),fill="#b8bac2",colour="gray60") +
  geom_point(data=m.tfit.plot.b1.ns, aes(x = logFC, y = AveExpr),fill="royalblue4",colour="royalblue4") +
  geom_point(data=m.tfit.plot.b1.de, aes(x = logFC, y = AveExpr,shape=de,colour=de),fill="royalblue4",colour="royalblue4") +
  scale_shape_manual(name="Expression",values=c(25,24)) +
  labs(x = "Expression log-ratio (B+ v B-)", y = "Average log-expression (TPM)", title = "Males") +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) + guides(shape = FALSE) + theme_bw()

p2 <- ggplot() +
  geom_vline(xintercept=(0)) +
  geom_point(data=f.tfit.plot.ns, aes(x = logFC, y = AveExpr),fill="#c2b8b8",colour="#c2b8b8", size=1.5,alpha=0.5) +
  geom_point(data=f.tfit.plot.de, aes(x = logFC, y = AveExpr,shape=de),fill="#c2b8b8",colour="gray60") +
  geom_point(data=f.tfit.plot.b1.ns, aes(x = logFC, y = AveExpr),fill="royalblue4",colour="royalblue4") +
  geom_point(data=f.tfit.plot.b1.de, aes(x = logFC, y = AveExpr,shape=de,colour=de),fill="royalblue4",colour="royalblue4") +
  scale_shape_manual(name="Expression",values=c(25,24)) +
  labs(x = "Expression log-ratio (B+ v B-)", y = "Average log-expression (TPM)", title = "Females") +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) + guides(shape = FALSE) + theme_bw()

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

counts.a <- data.frame(with(over.bm.genes, table(b.status.final)))
counts.b <- data.frame(with(under.bm.genes, table(b.status.final)))
counts.a$dir <- "O"
counts.b$dir <- "U"
counts <- rbind(counts.a,counts.b)
data.frame()

Status <- as.factor(c('A','B1','B2/B3','A','B1','B2/B3'))
Total <- as.integer(c(82,1,5,8,0,1))
Anno <- as.factor(c("82 (38)","1 (1)","5 (2)","8 (5)","0","1 (1)"))
Dir <- as.factor(c("Overexpressed","Overexpressed","Overexpressed","Underexpressed","Underexpressed","Underexpressed"))
counts <- data.frame(Status, Total, Anno, Dir)
counts

p4 <- ggplot(counts,aes(Status, Total, fill=Dir)) + ylim(0,90) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(label=Anno),position = position_dodge(0.9),vjust=-1) +
  scale_fill_manual(name="",values=c("gray20","gray60")) +
  labs(x = "Location on scaffold", y = "Number of genes", title = "DE in B+ males") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (12)), 
        axis.title = element_text(family = "Helvetica", size = (11)),
        axis.text = element_text(family = "Helvetica", size = (10)),
        legend.position=c(0.75,0.8), legend.title = element_blank(), legend.text = element_text(family = "Helvetica", size = (10)))

### Go analyses: how does having a B change your expression profiles if you are a male or a female?

# Make a compatible GO annotation file

Pviburni_genes_with_GO <- read_table2("output/pviburni.gene.GO", 
                                      col_names = FALSE)
colnames(Pviburni_genes_with_GO) <- c("gene","go")
Pviburni_genes_with_GO <- separate_rows(Pviburni_genes_with_GO, go, sep =';')

# Make the expression gene lists

MlogFC_DEgenes <- MB.vs.MnoB.results
FlogFC_DEgenes <- FB.vs.FnoB.results
MB_DEgenes <- MlogFC_DEgenes[MlogFC_DEgenes$adj.P.Val < 0.05,]
FB_DEgenes <- FlogFC_DEgenes[FlogFC_DEgenes$adj.P.Val < 0.05,]

nrow(MB_DEgenes)
nrow(FB_DEgenes)
background <- MlogFC_DEgenes[c(1)] #the background pops are identical so MlogFC_DEgenes[c(1)] = FlogFC_DEgenes[c(1)]
colnames(background)[1] <- "gene"
head(background)
background <- merge(Pviburni_genes_with_GO, background)

# load packages

#BiocManager::install("GOstats")
#BiocManager::install("treemap")
library(GOstats)
library(GSEABase)
library(treemap)

# Read in background GO set and make compatible with GOstats

GO_annotations <- background
GO_annotations[,3] <- paste("IEA")
names(GO_annotations) <- c("genes","GOIds","evi")
GO_annotations[,3] <- paste("IEA")
GO_annotations <- GO_annotations[c(2,3,1)]

# Create necessary objects

GO_frame <- GOFrame(GO_annotations,organism = "Pseudococcus viburni")
goAllFrame <- GOAllFrame(GO_frame)
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
universe <- as.vector(unique(GO_annotations[,3]))
nrow(background)

# Read in genes of interest 

nrow(MB_DEgenes)
m_genes <- as.data.frame(MB_DEgenes$gene)
m_genes <- as.data.frame(na.omit(m_genes[,1]))
m_genes <- as.vector(m_genes[,1])
length(m_genes)

f_genes <- as.data.frame(FB_DEgenes$gene)
f_genes <- as.data.frame(na.omit(f_genes[,1]))
f_genes <- as.vector(f_genes[,1])
length(f_genes)

# Keep only genes with annotated GOs
m_genes <- m_genes[m_genes %in% universe]
f_genes <- f_genes[f_genes %in% universe]
length(m_genes)
length(f_genes)
# background pop: 7375 annotated
# female biased genes: 122
# male biased genes: 119

# run hypergeometric test in M and collect results

Get_GO_params_all <- function(genes_of_i,universe,pvalue_cut){
  onto_terms <- c("BP","CC","MF")
  directions <- c("over","under")
  param_list <- list()
  name_1 <- list()
  for(i in 1:3){
    for(j in 1:2){
      name_1 <- c(name_1,paste(onto_terms[i],directions[j],sep = "_"))
      parameters <- GSEAGOHyperGParams(name="Hygeo params",
                                       geneSetCollection = gsc,
                                       universeGeneIds = universe,
                                       geneIds = m_genes,
                                       ontology = paste(onto_terms[i]),
                                       pvalueCutoff = pvalue_cut,
                                       conditional = T,testDirection = paste(directions[j]))
      param_list <- c(param_list,parameters)
    }
  }
  names(param_list) <- name_1
  return(param_list)
}

param_list <- Get_GO_params_all(genes_of_i = DE_Genes_A,universe = universe,
                                pvalue_cut = 0.01)

Hyper_G_test <- function(param_list){
  Hyper_G_list <- list()
  for(i in 1:length(param_list)){
    res <- hyperGTest(param_list[[i]])
    Hyper_G_list <- c(Hyper_G_list,res)
  }
  names(Hyper_G_list) <- names(param_list)
  return(Hyper_G_list)
}

GO_enrichment.M <- Hyper_G_test(param_list = param_list)

Result.BP.M <- summary(GO_enrichment.M[["BP_over"]])
Result.CC.M <- summary(GO_enrichment.M[["CC_over"]])
Result.MF.M <- summary(GO_enrichment.M[["MF_over"]])

colnames(Result.BP.M)[1] <- "GO"
colnames(Result.CC.M)[1] <- "GO"
colnames(Result.MF.M)[1] <- "GO"

Result.BP.M$Category <- "BP"
Result.CC.M$Category <- "CC"
Result.MF.M$Category <- "MF"
GO.enriched.M <- rbind(Result.BP.M,Result.CC.M,Result.MF.M)

# run hypergeometric test in F and collect results

Get_GO_params_all <- function(genes_of_i,universe,pvalue_cut){
  onto_terms <- c("BP","CC","MF")
  directions <- c("over","under")
  param_list <- list()
  name_1 <- list()
  for(i in 1:3){
    for(j in 1:2){
      name_1 <- c(name_1,paste(onto_terms[i],directions[j],sep = "_"))
      parameters <- GSEAGOHyperGParams(name="Hygeo params",
                                       geneSetCollection = gsc,
                                       universeGeneIds = universe,
                                       geneIds = f_genes,
                                       ontology = paste(onto_terms[i]),
                                       pvalueCutoff = pvalue_cut,
                                       conditional = T,testDirection = paste(directions[j]))
      param_list <- c(param_list,parameters)
    }
  }
  names(param_list) <- name_1
  return(param_list)
}

param_list <- Get_GO_params_all(genes_of_i = DE_Genes_A,universe = universe,
                                pvalue_cut = 0.01)

Hyper_G_test <- function(param_list){
  Hyper_G_list <- list()
  for(i in 1:length(param_list)){
    res <- hyperGTest(param_list[[i]])
    Hyper_G_list <- c(Hyper_G_list,res)
  }
  names(Hyper_G_list) <- names(param_list)
  return(Hyper_G_list)
}

GO_enrichment.F <- Hyper_G_test(param_list = param_list)

Result.BP.F <- summary(GO_enrichment.F[["BP_over"]])
Result.CC.F <- summary(GO_enrichment.F[["CC_over"]])
Result.FF.F <- summary(GO_enrichment.F[["MF_over"]])

colnames(Result.BP.F)[1] <- "GO"
colnames(Result.CC.F)[1] <- "GO"
colnames(Result.FF.F)[1] <- "GO"

Result.BP.F$Category <- "BP"
Result.CC.F$Category <- "CC"
Result.FF.F$Category <- "MF"
GO.enriched.F <- rbind(Result.BP.F,Result.CC.F,Result.FF.F)

#write.csv(GO.enriched.F, file="output/B_diff_expr/GO.enriched.BvsnoB.F.csv", quote = F, row.names = F)
#write.csv(GO.enriched.M, file="output/B_diff_expr/GO.enriched.BvsnoB.M.csv", quote = F, row.names = F)
p1 + p2 + p3.blank + p4
library(patchwork)
tiff("E:/agdel/Documents/projects_rosslab/B_viburni/manuscript/figures/fig4.tiff",
     width = 3000, height = 2500, units = 'px', res = 300)
p1 + p2 + p3.blank + p4
dev.off()

tiff("E:/agdel/Documents/projects_rosslab/B_viburni/manuscript/figures/fig4c.tiff",
     width = 2700, height = 2700, units = 'px', res = 300)
vennDiagram(dt[,1:3], circle.col=c("purple", "green","orange"),include=c("up","down"))
dev.off()

