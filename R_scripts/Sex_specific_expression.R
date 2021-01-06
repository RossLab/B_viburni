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

# import counts and sample info

setwd("E:\agdel\Documents\projects_rosslab\B_viburni\R_scripts")

rsem.counts <-read.delim("RSEM_digi.counts.matrix",header=TRUE) #matrix generated from rsem

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

tfit <- treat(efit, lfc=1) 
dt <- decideTests(tfit)
summary(dt)
write.fit(tfit, dt, file="diff_expression_sex_tfit.txt")









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
