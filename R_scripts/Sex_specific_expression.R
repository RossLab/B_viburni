rm(list=ls())
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("edgeR")
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

setwd("E:/agdel/Documents/projects_rosslab/B_viburni")


##### Differential expression analysis

rsem.counts <-read.delim("R_scripts/RSEM_digi.counts.matrix",header=TRUE) #matrix generated from rsem
sampleinfo <- read.csv("sampleinfoPviburnisex.csv") #sample group info

a <- rsem.counts
head(a,2)

colnames(a) <- substr(colnames(a),start=1,stop=6) #removing ".genes.results" in colnames
colnames(a)

# gene id as row name

head(a) 
a1 <- a[,-1]
head(a1)
rownames(a1) <- a[,1]
head(a1)
a<-a1
head(a)

# adding replicate information

group <- as.factor(c("F","F","F","M","M","M","F","F","F","M","M","M","M","F","F","F","M","M","M","F","F","F","M","M","M","M"))

# create the DGE object including all count matrix, replicate info and annotations when available

x <- DGEList(counts=round(a),genes=rownames(a), group = group)
x$samples

# filter out low count

rowSums(x$counts==0) # for each sample where the count is 0, then sums the samples that have counts = 0 for a gene
table(rowSums(x$counts==0) == 26) #number of genes that have all samples with 0 counts
proportion_0_count = 21236/2393


dim(x)

# filter with default options

keep.exprs.group <- filterByExpr(x, group= group) 
x1 <- x[keep.exprs.group, keep.lib.sizes=FALSE]
dim(x1)
x<-x1

# normalize distribution (TMM normalization)
x <- calcNormFactors(x, method = "TMM")

x$samples$norm.factors

# plot individual samples post normalization

for (i in 1:ncol(x)) {
	plotMD(cpm(x, log=TRUE), column=i, xlab = "Average log-expression", ylab = "Expression log-ratio (this sample vs others)", main = colnames(x)[i])
	abline(h=0, col="red", lty=2, lwd=2)
}

# QC, MDS, heatmap plots

par(mfrow=c(3,1))

barplot(x$samples$lib.size,names=colnames(x),las=2)
title("Barplot of library sizes")
logcounts <- cpm(x,log=TRUE)

boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (normalized)")

sampleinfo2 <- sampleinfo %>% mutate(Sex = as.factor(Sex))
levels(sampleinfo2$Sex)
col.sex <- c("red","blue")[sampleinfo2$Sex]
data.frame(sampleinfo2$Sex,col.sex)
plotMDS(x,col=col.sex, main="Sex")
legend("top",fill=c("red","blue"),legend=levels(sampleinfo2$Sex))

var_genes <- apply(logcounts, 1, var)
head(var_genes)
select_var <- names(sort(var_genes, decreasing=TRUE))[1:1000]
head(select_var)
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)

library(RColorBrewer)
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 1000 most variable genes across samples",ColSideColors=col.sex,scale="row")

# design model matrix

design <- model.matrix(~0 + group)
colnames(design)
rownames(design) = rownames(x$samples)

# removing heteroscedascity

par(mfrow=c(2,1))
v1 <- voom(x,design,plot = TRUE)
v1

# run model

contr.matrix <- makeContrasts(FvsM = groupF - groupM, levels = design)
vfit <- lmFit(v1, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")
summary(decideTests(efit))

# set logFc to 1

tfit <- treat(efit, lfc = 1) 
dt <- decideTests(tfit)
summary(dt)
#write.fit(tfit, dt, file="output/sex_diff_expr/FvsM_results_tfit.txt") #export results

# results

FvsM.results.all <- topTreat(tfit, coef=1, n=Inf)

FvsM.anno <- FvsM.results.all
colnames(FvsM.anno)[1] <- "gene"
FvsM.anno$de <- "NS"
FvsM.anno$de <- ifelse(FvsM.anno$adj.P.Val < 0.05 & FvsM.anno$logFC < 0,
                                   "MB",FvsM.anno$de)
FvsM.anno$de <- ifelse(FvsM.anno$adj.P.Val < 0.05 & FvsM.anno$logFC > 0,
                                   "FB",FvsM.anno$de)
table(FvsM.anno$de)

# plot

sb.plot <- ggplot(FvsM.anno, aes(x = logFC, y = AveExpr)) +
  geom_vline(xintercept=(0)) +
  geom_point(aes(colour=de, size=de, alpha=de)) +
  scale_colour_manual(name="Expression",values=c("firebrick3","cornflowerblue","azure4")) +
  scale_size_manual(values=c(1, 1, 0.8),guide=FALSE) +
  scale_alpha_manual(values=c(0.8, 0.8, 0.3),guide=FALSE) +
  labs(x = "Expression log-ratio (F v M)", y = "Average log-expression", title = "Sex biased expression") +
  theme_bw()
table(FvsM.anno$de)
#tiff("E:/agdel/Documents/projects_rosslab/B_viburni/manuscript/plots/FIG.sex.biased.expr.tiff",
#    width = 3000, height = 2100, units = 'px', res = 300)
#sb.plot
#dev.off()
#
#jpeg("E:/agdel/Documents/projects_rosslab/B_viburni/misc/FIG.sex.biased.expr.jpg",
#     width = 3000, height = 2100, units = 'px', res = 300)
#sb.plot
#dev.off()

# annotate genes with the master annotation

freeze.v0.genes.anno <- read_delim("output/freeze.v0.genes.anno.complete.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE) # master anno
FvsM.anno <- left_join(FvsM.anno,freeze.v0.genes.anno,by="gene")

#write.csv(FvsM.anno, file="output/sex_diff_expr/FvsM_results_anno.csv") #export results


##### GO analysis

# Make a compatible GO annotation file

Pcitri_genes_with_GO <- read_table2("output/pviburni.gene.GO", 
                                    col_names = FALSE)
colnames(Pcitri_genes_with_GO) <- c("gene","go")
Pcitri_genes_with_GO <- separate_rows(Pcitri_genes_with_GO, go, sep =';')
#write.table(new_annotations, file="P_citri_GO_terms.txt", sep="\t", quote = F,
#            col.names = T, row.names = F)


# Make the expression gene lists

logFC_DEgenes <- FvsM.anno
nrow(FB_DEgenes)

MB_DEgenes <- logFC_DEgenes[logFC_DEgenes$de == "MB",]
FB_DEgenes <- logFC_DEgenes[logFC_DEgenes$de == "FB",]

background <- merge(Pcitri_genes_with_GO, logFC_DEgenes[c(1)])

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

GO_frame <- GOFrame(GO_annotations,organism = "Planococcus citri")
goAllFrame <- GOAllFrame(GO_frame)
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
universe <- as.vector(unique(GO_annotations[,3]))

# Read in genes of interest 

f_genes <- as.data.frame(logFC_DEgenes$gene[logFC_DEgenes$de=="FB"])
f_genes <- as.data.frame(na.omit(f_genes[,1]))
f_genes <- as.vector(f_genes[,1])

m_genes <- as.data.frame(logFC_DEgenes$gene[logFC_DEgenes$de=="MB"])
m_genes <- as.data.frame(na.omit(m_genes[,1]))
m_genes <- as.vector(m_genes[,1])

# Keep only genes with annotated GOs

f_genes <- f_genes[f_genes %in% universe]
m_genes <- m_genes[m_genes %in% universe]
# background pop: 6958 annotated / 15673 total
# female biased genes: 803/2186
# male biased genes: 653/1644

# Set up parameters for hypergeometric test

my_genes <- f_genes
#my_genes <- m_genes

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
                                       geneIds = my_genes,
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

# Hypergeometric test

Hyper_G_test <- function(param_list){
  Hyper_G_list <- list()
  for(i in 1:length(param_list)){
    res <- hyperGTest(param_list[[i]])
    Hyper_G_list <- c(Hyper_G_list,res)
  }
  names(Hyper_G_list) <- names(param_list)
  return(Hyper_G_list)
}

GO_enrichment <- Hyper_G_test(param_list = param_list)

# Extract results

Result.BP.F <- summary(GO_enrichment[["BP_over"]])
Result.CC.F <- summary(GO_enrichment[["CC_over"]])
Result.MF.F <- summary(GO_enrichment[["MF_over"]])

colnames(Result.BP.F)[1] <- "GO"
colnames(Result.CC.F)[1] <- "GO"
colnames(Result.MF.F)[1] <- "GO"

Result.BP.F$Category <- "BP"
Result.CC.F$Category <- "CC"
Result.MF.F$Category <- "MF"
GO.enriched.F <- rbind(Result.BP.F,Result.CC.F,Result.MF.F)

Result.BP.M <- summary(GO_enrichment[["BP_over"]])
Result.CC.M <- summary(GO_enrichment[["CC_over"]])
Result.MF.M <- summary(GO_enrichment[["MF_over"]])

colnames(Result.BP.M)[1] <- "GO"
colnames(Result.CC.M)[1] <- "GO"
colnames(Result.MF.M)[1] <- "GO"

Result.BP.M$Category <- "BP"
Result.CC.M$Category <- "CC"
Result.MF.M$Category <- "MF"
GO.enriched.M <- rbind(Result.BP.M,Result.MF.M)

# export results

#write.table(GO.enriched.F, file="output/sex_diff_expr/GO.enriched.F.tsv", sep = "\t", quote = F, row.names = F)
#write.table(GO.enriched.M, file="output/sex_diff_expr/GO.enriched.M.tsv", sep = "\t", quote = F, row.names = F)