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

# B chromosome assignment

scaffolds.final.assignment <- read_delim("output/scaffolds.final.assignment.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE)

# differentially expressed genes -- from Isabelle (received 03.11.20)

fb.vs.fnob <- read_delim("output/rsem_gene_femaleB.femalenoB.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE)
mb.vs.fb <- read_delim("output/rsem_gene_maleB.femaleB.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE)
mb.vs.fnob <- read_delim("output/rsem_gene_maleB.femalenoB.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE)
mb.vs.mnob <- read_delim("output/rsem_gene_maleB.noB_annot.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE)
mnob.vs.fnb <- read_delim("output/rsem_gene_malenoB.femalenoB.csv",",", escape_double = FALSE, col_names = T,trim_ws = TRUE)

# import assembly features and annotation hits

freeze.v0.genes.mapped <- read_delim("annotation/p.viburni.freeze.v0.braker.genes.mapped.txt", 
                                      "\t", escape_double = FALSE, col_names = FALSE, 
                                      trim_ws = TRUE)

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
colnames(freeze.v0.blast.interproscan)[6] <- "interpro"
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
freeze.v0.blast.uniprot.best$gene <- gsub("\\..*","", freeze.v0.blast.uniprot.best$transcript)
freeze.v0.blast.uniprot.best.dedup <- freeze.v0.blast.uniprot.best[c(3,2)][!duplicated(freeze.v0.blast.uniprot.best[c(3,2)]), ]
freeze.v0.blast.uniprot.best.dedup.nona <- freeze.v0.blast.uniprot.best.dedup[!is.na(freeze.v0.blast.uniprot.best.dedup$uniprot),]
freeze.v0.blast.uniprot.best.dedup.nona.collapsed <- ddply(freeze.v0.blast.uniprot.best.dedup.nona,c("gene"), summarize,
                                            blast = paste(uniprot,collapse=",")) # collapse annotations for genes with >1 descriptor
nrow(freeze.v0.blast.uniprot.best.dedup.nona)

# diamond
freeze.v0.diamond.refprot.best <- freeze.v0.diamond.refprot %>% distinct(transcript, .keep_all = TRUE)
freeze.v0.diamond.refprot.best$gene <- gsub("\\..*","", freeze.v0.diamond.refprot.best$transcript)
freeze.v0.diamond.refprot.best.dedup <- freeze.v0.diamond.refprot.best[c(3,2)][!duplicated(freeze.v0.diamond.refprot.best[c(3,2)]), ]
freeze.v0.diamond.refprot.best.dedup.nona <- freeze.v0.diamond.refprot.best.dedup[!is.na(freeze.v0.diamond.refprot.best.dedup$refprot),]
freeze.v0.diamond.refprot.best.dedup.nona.collapsed <- ddply(freeze.v0.diamond.refprot.best.dedup.nona,c("gene"), summarize,
                                                           diamond = paste(refprot,collapse=",")) # collapse annotations for genes with >1 descriptor



freeze.v0.genes <- merge(freeze.v0.blast.uniprot.best,freeze.v0.diamond.refprot.best,by="transcript",all=TRUE)


head(freeze.v0.genes)
head(scaffolds.final.assignment)
colnames(genes.in.B.strict.candidates)[1] <- "gene"
colnames(genes.in.B.strict.candidates)[2] <- "length"
colnames(genes.in.B.strict.candidates)[3] <- "seq"

genes.in.B.strict.anno <- merge(genes.in.B.strict.candidates,freeze.v0.genes,by="gene")
nrow(genes.in.B.strict.anno)

# transcriptome

freeze.v0.blast.transcriptome <- read_delim("p.viburni.freeze.v0.braker.aa.blast.vs.transcriptome.out", 
                                            "\t", escape_double = FALSE, col_names = FALSE, 
                                            trim_ws = TRUE)
colnames(freeze.v0.blast.transcriptome)[1] <- "transcript"
colnames(freeze.v0.blast.transcriptome)[2] <- "protein_id"
freeze.v0.blast.transcriptome$gene <- gsub("\\..*","", freeze.v0.blast.transcriptome$transcript)

genes.in.B.strict.transcriptome <- merge(genes.in.B.strict.candidates,freeze.v0.blast.transcriptome,by="gene")
nrow(genes.in.B.strict.transcriptome)

# import differentially expressed transcripts

overexpressed.B.transcripts <- read_csv("~/Documents/genomics/B_viburni_ross_lab/data/preliminary_rna_seq/sleuth_v1/overexpressed.B.transcripts.csv")
overexpressed.B.transcripts.cat <- overexpressed.B.transcripts[c("transcript_id","info","sprot_Top_BLASTX_hit","sprot_Top_BLASTP_hit","Pfam","eggnog","Kegg","gene_ontology_BLASTX","gene_ontology_BLASTP","gene_ontology_Pfam")]

genes.in.B.strict.transcriptome$transcript_id <- gsub("\\..*","", genes.in.B.strict.transcriptome$protein_id)
genes.in.B.strict.transcriptome.expr <- merge(genes.in.B.strict.transcriptome, overexpressed.B.transcripts.cat, by = "transcript_id")

count(unique(genes.in.B.strict.transcriptome$gene))
count(unique(genes.in.B.strict.transcriptome.expr$gene))

table(genes.in.B.strict.transcriptome.expr$info)