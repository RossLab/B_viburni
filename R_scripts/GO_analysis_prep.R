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

##### setwd and import data

setwd("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/data")
freeze.v0.blast.interproscan <- read_delim("annotation/p.viburni.freeze.v0.braker.interproscan.tsv", 
                                           "\t", escape_double = FALSE, col_names = FALSE, 
                                           trim_ws = TRUE) # interproscan file with GO terms
dt_df <- read_csv("~/Documents/genomics/B_viburni_ross_lab/output/diff_expr/dt_df.csv") # differential expression analysis

##### Prepare GO association file

freeze.v0.blast.interproscan.go <- freeze.v0.blast.interproscan[c(1,14)]
colnames(freeze.v0.blast.interproscan.go)[1] <- "transcript"
colnames(freeze.v0.blast.interproscan.go)[2] <- "GO"
nrow(freeze.v0.blast.interproscan.go)
freeze.v0.blast.interproscan.go <- freeze.v0.blast.interproscan.go[!is.na(freeze.v0.blast.interproscan.go$GO),] # remove transcripts without GO anno
nrow(freeze.v0.blast.interproscan.go)
freeze.v0.blast.interproscan.go <- unique(freeze.v0.blast.interproscan.go) # remove duplicate entries

freeze.v0.blast.interproscan.go.split <- freeze.v0.blast.interproscan.go %>% 
  mutate(GO = strsplit(as.character(GO), "|", fixed = TRUE)) %>% 
  unnest(GO) # split multiple go terms in individual rows
freeze.v0.blast.interproscan.go.split$gene <- gsub("\\..*","", freeze.v0.blast.interproscan.go.split$transcript)
pviburni.gene.GO.split <- freeze.v0.blast.interproscan.go.split[c(3,2)]
pviburni.gene.GO.split <- unique(pviburni.gene.GO.split) # remove duplicate entries again

pviburni.gene.GO <- ddply(pviburni.gene.GO.split,c("gene"), summarize,
                          go=paste(GO,collapse=";"))
nrow(pviburni.gene.GO)

##### Prepare GO files

background.pop.go <- dt_df[c("gene")][dt_df$gene %in% pviburni.gene.GO$gene,]
DE.males.genes.go <- dt_df[c("gene")][dt_df$BmalevnoBmale != "0",]
DE.males.genes.go <- DE.males.genes.go[DE.males.genes.go$gene %in% pviburni.gene.GO$gene,]
DE.females.genes.go <- dt_df[c("gene")][dt_df$BfemalevsnoBfemale != "0",]
DE.females.genes.go <- DE.females.genes.go[DE.females.genes.go$gene %in% pviburni.gene.GO$gene,]

write.table(background.pop.go, file = "output/diff_expr/go/background.pop.go",col.names= FALSE,row.names=FALSE, quote = FALSE)
write.table(DE.males.genes.go, file = "output/diff_expr/go/DE.males.genes.go",col.names= FALSE,row.names=FALSE, quote = FALSE)
write.table(DE.females.genes.go, file = "output/diff_expr/go/DE.females.genes.go",col.names= FALSE,row.names=FALSE, quote = FALSE)
write.table(pviburni.gene.GO, file = "output/diff_expr/go/pviburni.gene.GO",col.names= FALSE,row.names=FALSE, quote = FALSE,sep="\t")
