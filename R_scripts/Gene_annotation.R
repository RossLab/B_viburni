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

setwd("E:/agdel/Documents/projects_rosslab/B_viburni")

# import assembly features and annotation hits

freeze.v0.genes.mapped <- read_delim("annotation/p.viburni.freeze.v0.braker.genes.mapped.txt",
                                     "\t", escape_double = FALSE, trim_ws = TRUE)

freeze.v0.blast.uniprot <- read_delim("annotation/p.viburni.freeze.v0.braker.aa.blast.vs.uniprot.out", 
                                      "\t", escape_double = FALSE, col_names = FALSE, 
                                      trim_ws = TRUE)

freeze.v0.blast.uniprot <- freeze.v0.blast.uniprot[c(1,2)]
colnames(freeze.v0.blast.uniprot)[1] <- "transcript"
colnames(freeze.v0.blast.uniprot)[2] <- "uniprot"
freeze.v0.blast.uniprot$gene <- gsub("\\..*","", freeze.v0.blast.uniprot$transcript)

freeze.v0.diamond.refprot <- read_delim("annotation/p.viburni.freeze.v0.braker.aa.diamond.vs.refprot.out", 
                                        "\t", escape_double = FALSE, col_names = FALSE, 
                                        trim_ws = TRUE)

freeze.v0.diamond.refprot <- freeze.v0.diamond.refprot[c(1,2)]
colnames(freeze.v0.diamond.refprot)[1] <- "transcript"
colnames(freeze.v0.diamond.refprot)[2] <- "refprot"
freeze.v0.diamond.refprot$gene <- gsub("\\..*","", freeze.v0.diamond.refprot$transcript)

freeze.v0.interproscan <- read_delim("annotation/p.viburni.freeze.v0.braker.interproscan.tsv", 
                                            "\t", escape_double = FALSE, col_names = FALSE, 
                                            trim_ws = TRUE)

colnames(freeze.v0.interproscan)[1] <- "transcript"
colnames(freeze.v0.interproscan)[5] <- "pfam_acc"
colnames(freeze.v0.interproscan)[6] <- "pfam_descr"
colnames(freeze.v0.interproscan)[12] <- "ipr_acc"
colnames(freeze.v0.interproscan)[13] <- "ipr_descr"
colnames(freeze.v0.interproscan)[14] <- "go"
freeze.v0.interproscan$gene <- gsub("\\..*","", freeze.v0.interproscan$transcript)

# make gene annotation file, collapsing by gene and keeping first (best) blast/diamond hits

# blast
freeze.v0.blast.uniprot.best <- freeze.v0.blast.uniprot %>% distinct(transcript, .keep_all = TRUE) # will keep the first (best) match
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

freeze.v0.interproscan <- freeze.v0.interproscan[freeze.v0.interproscan$X4 != "SignalP_EUK",] #remove entries annotated with SignalP, not relevant

# we have to separate pfam, ipr and go terms to merge later

freeze.v0.pfam <- freeze.v0.interproscan[c(15,5,6)]
freeze.v0.pfam.dedup <- freeze.v0.pfam[!duplicated(freeze.v0.pfam), ]
freeze.v0.pfam.dedup.nona <- freeze.v0.pfam.dedup[!is.na(freeze.v0.pfam.dedup$pfam_descr),] # remove entries with NA in description -- there should be none though
freeze.v0.pfam.dedup.nona.collapsed <- ddply(freeze.v0.pfam.dedup.nona,c("gene"), summarize,
                                        pfam_acc = paste(pfam_acc,collapse=";"),
                                        pfam_descr = paste(pfam_descr,collapse="; ")) # collapse annotations for genes with >1 descriptor
nrow(freeze.v0.pfam.collapsed)

freeze.v0.ipr <- freeze.v0.interproscan[c(15,12,13)]
freeze.v0.ipr.dedup <- freeze.v0.ipr[!duplicated(freeze.v0.ipr), ]
freeze.v0.ipr.dedup.nona <- freeze.v0.ipr.dedup[!is.na(freeze.v0.ipr.dedup$ipr_acc),] # remove entries with NA in description
freeze.v0.ipr.dedup.nona.collapsed <- ddply(freeze.v0.ipr.dedup.nona,c("gene"), summarize,
                                        ipr_acc = paste(ipr_acc,collapse=";"),
                                        ipr_descr = paste(ipr_descr,collapse="; ")) # collapse annotations for genes with >1 descriptor

freeze.v0.go <- freeze.v0.interproscan[c(15,14)]
freeze.v0.go.dedup <- freeze.v0.go[!duplicated(freeze.v0.go), ]
nrow(freeze.v0.go)
freeze.v0.go.dedup.nona <- freeze.v0.go.dedup[!is.na(freeze.v0.go.dedup$go),] # remove transcripts without GO anno
nrow(freeze.v0.go.dedup)
freeze.v0.go.split <- freeze.v0.go.dedup.nona %>% 
  mutate(GO = strsplit(as.character(go), "|", fixed = TRUE)) %>% 
  unnest(GO) # split multiple go terms in individual rows
freeze.v0.go.split <- freeze.v0.go.split[c(1,3)]
freeze.v0.go.split.dedup <- unique(freeze.v0.go.split) # remove duplicate entries again
freeze.v0.go.split.dedup.nona <- freeze.v0.go.split.dedup[!is.na(freeze.v0.go.split.dedup$GO),] # remove entries with NA -- there should be none though

freeze.v0.go.split.dedup.nona.collapsed <- ddply(freeze.v0.go.split.dedup.nona,c("gene"), summarize,
                          GO=paste(GO,collapse=";"))
nrow(freeze.v0.go.split.dedup.nona.collapsed)

# make and export master annotation

freeze.v0.genes <- merge(freeze.v0.genes.mapped, freeze.v0.blast.uniprot.best.dedup.nona.collapsed, by="gene", all=TRUE)
freeze.v0.genes <- merge(freeze.v0.genes, freeze.v0.diamond.refprot.best.dedup.nona.collapsed,by="gene",all=TRUE)
freeze.v0.genes <- merge(freeze.v0.genes, freeze.v0.pfam.dedup.nona.collapsed,by="gene",all=TRUE)
freeze.v0.genes <- merge(freeze.v0.genes, freeze.v0.ipr.dedup.nona.collapsed,by="gene",all=TRUE)
freeze.v0.genes <- merge(freeze.v0.genes, freeze.v0.go.split.dedup.nona.collapsed,by="gene",all=TRUE)


freeze.v0.genes$anno <- ifelse((!is.na(freeze.v0.genes$blast) | !is.na(freeze.v0.genes$diamond) | !is.na(freeze.v0.genes$pfam_acc) |
                                  !is.na(freeze.v0.genes$ipr_acc) | !is.na(freeze.v0.genes$GO)),"Y","N")

write.table(freeze.v0.genes, file = "output/freeze.v0.genes.anno.complete.csv",row.names = F,sep = ",")
