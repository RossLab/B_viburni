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

# setwd and import data

setwd("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/data")

# import assembly features and annotation hits

freeze.v0.genes.mapped <- read_delim("annotation/p.viburni.freeze.v0.braker.genes.mapped.txt",
                                     "\t", escape_double = FALSE, trim_ws = TRUE)

freeze.v0.blast.uniprot <- read_delim("annotation/p.viburni.freeze.v0.braker.aa.blast.vs.uniprot.out", 
                                      "\t", escape_double = FALSE, col_names = FALSE, 
                                      trim_ws = TRUE)

freeze.v0.blast.uniprot <- freeze.v0.blast.uniprot[c(1,2)]
colnames(freeze.v0.blast.uniprot)[1] <- "transcript"
colnames(freeze.v0.blast.uniprot)[2] <- "uniprot"
freeze.v0.blast.uniprot.best$gene <- gsub("\\..*","", freeze.v0.blast.uniprot.best$transcript)

freeze.v0.diamond.refprot <- read_delim("annotation/p.viburni.freeze.v0.braker.aa.diamond.vs.refprot.out", 
                                        "\t", escape_double = FALSE, col_names = FALSE, 
                                        trim_ws = TRUE)

freeze.v0.diamond.refprot <- freeze.v0.diamond.refprot[c(1,2)]
colnames(freeze.v0.diamond.refprot)[1] <- "transcript"
colnames(freeze.v0.diamond.refprot)[2] <- "refprot"
freeze.v0.diamond.refprot.best$gene <- gsub("\\..*","", freeze.v0.diamond.refprot.best$transcript)

freeze.v0.blast.interproscan <- read_delim("annotation/p.viburni.freeze.v0.braker.interproscan.tsv", 
                                            "\t", escape_double = FALSE, col_names = FALSE, 
                                            trim_ws = TRUE)
colnames(freeze.v0.blast.interproscan)[1] <- "transcript"
colnames(freeze.v0.blast.interproscan)[6] <- "interproscan"
freeze.v0.blast.interproscan$gene <- gsub("\\..*","", freeze.v0.blast.interproscan$transcript)

freeze.v0.blast.transcriptome <- read_delim("annotation/p.viburni.freeze.v0.braker.aa.blast.vs.transcriptome.out", 
                                            "\t", escape_double = FALSE, col_names = FALSE, 
                                            trim_ws = TRUE)
colnames(freeze.v0.blast.transcriptome)[1] <- "transcript"
colnames(freeze.v0.blast.transcriptome)[2] <- "protein_id"
freeze.v0.blast.transcriptome$gene <- gsub("\\..*","", freeze.v0.blast.transcriptome$transcript)

# make gene annotation file, collapsing by gene and keeping first (best) blast/diamond hits

# blast
freeze.v0.blast.uniprot.best <- freeze.v0.blast.uniprot %>% distinct(transcript, .keep_all = TRUE)
freeze.v0.blast.uniprot.best.dedup <- freeze.v0.blast.uniprot.best[c(3,2)][!duplicated(freeze.v0.blast.uniprot.best[c(3,2)]), ]
freeze.v0.blast.uniprot.best.dedup.nona <- freeze.v0.blast.uniprot.best.dedup[!is.na(freeze.v0.blast.uniprot.best.dedup$uniprot),]
freeze.v0.blast.uniprot.best.dedup.nona.collapsed <- ddply(freeze.v0.blast.uniprot.best.dedup.nona,c("gene"), summarize,
                                            blast = paste(uniprot,collapse=",")) # collapse annotations for genes with >1 descriptor


# diamond
freeze.v0.diamond.refprot.best <- freeze.v0.diamond.refprot %>% distinct(transcript, .keep_all = TRUE)
freeze.v0.diamond.refprot.best.dedup <- freeze.v0.diamond.refprot.best[c(3,2)][!duplicated(freeze.v0.diamond.refprot.best[c(3,2)]), ]
freeze.v0.diamond.refprot.best.dedup.nona <- freeze.v0.diamond.refprot.best.dedup[!is.na(freeze.v0.diamond.refprot.best.dedup$refprot),]
freeze.v0.diamond.refprot.best.dedup.nona.collapsed <- ddply(freeze.v0.diamond.refprot.best.dedup.nona,c("gene"), summarize,
                                                           diamond = paste(refprot,collapse=",")) # collapse annotations for genes with >1 descriptor

# interproscan
freeze.v0.blast.interproscan.dedup <- freeze.v0.blast.interproscan[c(15,6)][!duplicated(freeze.v0.blast.interproscan[c(15,6)]), ]
freeze.v0.blast.interproscan.dedup.nona <- freeze.v0.blast.interproscan.dedup[!is.na(freeze.v0.blast.interproscan.dedup$interproscan),]
freeze.v0.blast.interproscan.dedup.nona.collapsed <- ddply(freeze.v0.blast.interproscan.dedup.nona,c("gene"), summarize,
                                                             interpro = paste(interproscan,collapse=",")) # collapse annotations for genes with >1 descriptor
nrow(freeze.v0.blast.interproscan.dedup.nona.collapsed)

# make and export master annotation
 
freeze.v0.genes <- merge(freeze.v0.genes.mapped, freeze.v0.blast.uniprot.best.dedup.nona.collapsed, by="gene", all=TRUE)
freeze.v0.genes <- merge(freeze.v0.genes, freeze.v0.diamond.refprot.best.dedup.nona.collapsed,by="gene",all=TRUE)
freeze.v0.genes <- merge(freeze.v0.genes, freeze.v0.blast.interproscan.dedup.nona.collapsed,by="gene",all=TRUE)
head(freeze.v0.genes)
(nrow(freeze.v0.genes) - nrow(freeze.v0.genes[is.na(freeze.v0.genes$blast) & is.na(freeze.v0.genes$diamond) & is.na(freeze.v0.genes$interpro),])) / nrow(freeze.v0.genes)

write.table(freeze.v0.genes, file = "output/freeze.v0.genes.anno.csv",row.names = F,sep = ",")
