rm(list=ls())
ls()
library(tidyverse)
library(gridExtra)
library(patchwork)
library(lattice)
library(grid)
library(gridExtra)
library(reshape2)

setwd("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/data/coverage_analysis")

# import number of reads mapped to scaffold

PV04.reads.mapped <- read_delim("PV_18-04.primary.reads.mapped.no.mismatches.count", 
                                "\t", escape_double = FALSE, col_names = FALSE, 
                                trim_ws = TRUE)
PV13.reads.mapped <- read_delim("PV_18-13.primary.reads.mapped.no.mismatches.count", 
                                "\t", escape_double = FALSE, col_names = FALSE, 
                                trim_ws = TRUE)
PV21.reads.mapped <- read_delim("PV_18-21.primary.reads.mapped.no.mismatches.count", 
                                "\t", escape_double = FALSE, col_names = FALSE, 
                                trim_ws = TRUE)
PV23.reads.mapped <- read_delim("PV_18-23.primary.reads.mapped.no.mismatches.count", 
                                "\t", escape_double = FALSE, col_names = FALSE, 
                                trim_ws = TRUE)

colnames(PV04.reads.mapped) <- c("seq","length","PV04.mapped","unmapped")
colnames(PV13.reads.mapped) <- c("seq","length","PV13.mapped","unmapped")
colnames(PV21.reads.mapped) <- c("seq","length","PV21.mapped","unmapped")
colnames(PV23.reads.mapped) <- c("seq","length","PV23.mapped","unmapped")

# explore differences in coverage. Let's always do +1 to avoid 0s. We can plot this in all pairs of samples. 

cov.13v21 <- log2((PV13.reads.mapped$PV13.mapped+1)/(PV21.reads.mapped$PV21.mapped+1))
cov.13v23 <- log2((PV13.reads.mapped$PV13.mapped+1)/(PV23.reads.mapped$PV23.mapped+1))
cov.04v21 <- log2((PV04.reads.mapped$PV04.mapped+1)/(PV21.reads.mapped$PV21.mapped+1))
cov.04v23 <- log2((PV04.reads.mapped$PV04.mapped+1)/(PV23.reads.mapped$PV23.mapped+1))
cov.04v13 <- log2((PV04.reads.mapped$PV04.mapped+1)/(PV13.reads.mapped$PV13.mapped+1))
cov.21v23 <- log2((PV21.reads.mapped$PV21.mapped+1)/(PV23.reads.mapped$PV23.mapped+1))

cov.diff <- data.frame(cov.13v21,cov.13v23,cov.04v21,cov.04v23,cov.04v13,cov.21v23)

p1 <- ggplot(cov.diff, aes(cov.13v21)) + geom_histogram(bins=150)
p2 <- ggplot(cov.diff, aes(cov.13v23)) + geom_histogram(bins=150)
p3 <- ggplot(cov.diff, aes(cov.04v21)) + geom_histogram(bins=150)
p4 <- ggplot(cov.diff, aes(cov.04v23)) + geom_histogram(bins=150)
p5 <- ggplot(cov.diff, aes(cov.04v13)) + geom_histogram(bins=150)
p6 <- ggplot(cov.diff, aes(cov.21v23)) + geom_histogram(bins=150)

hist1 <- p1 + p2 + p3 + p4
hist2 <- p5 + p6

# looks promising. Let's create a datafile

seq <- PV04.reads.mapped[ ,1]
length <- PV04.reads.mapped[ ,2]
PV04 <- PV04.reads.mapped[ ,3]
PV13 <- PV13.reads.mapped[ ,3]
PV21 <- PV21.reads.mapped[ ,3]
PV23 <- PV23.reads.mapped[ ,3]
reads.all.lines0 <- data.frame(seq,length,PV04,PV13,PV21,PV23)
reads.all.lines <-reads.all.lines0[1:(nrow(reads.all.lines0)-1),] # remove last line (not a scaffold)

# normalise by total number of reads

sum(reads.all.lines[, 'PV04.mapped']) # 180975436
sum(reads.all.lines[, 'PV13.mapped']) # 186868354
sum(reads.all.lines[, 'PV21.mapped']) # 88933548
sum(reads.all.lines[, 'PV23.mapped']) # 95633941

# normalisation factor

cov.13v21_norm <- round((PV13.reads.mapped$PV13.mapped+1)/(PV21.reads.mapped$PV21.mapped+1),2)
cov.13v23_norm <- round((PV13.reads.mapped$PV13.mapped+1)/(PV23.reads.mapped$PV23.mapped+1),2)
cov.04v21_norm <- round((PV04.reads.mapped$PV04.mapped+1)/(PV21.reads.mapped$PV21.mapped+1),2)
cov.04v23_norm <- round((PV04.reads.mapped$PV04.mapped+1)/(PV23.reads.mapped$PV23.mapped+1),2)
cov.04v13_norm <- round((PV04.reads.mapped$PV04.mapped+1)/(PV13.reads.mapped$PV13.mapped+1),2)
cov.21v23_norm <- round((PV21.reads.mapped$PV21.mapped+1)/(PV23.reads.mapped$PV23.mapped+1),2)
median(cov.13v21_norm)
median(cov.13v23_norm)
median(cov.04v21_norm)
median(cov.04v23_norm)
median(cov.04v13_norm)
median(cov.21v23_norm)

norm.04 <- sum(reads.all.lines[, 'PV21.mapped']) / sum(reads.all.lines[, 'PV04.mapped']) # 0.49
norm.13 <- sum(reads.all.lines[, 'PV21.mapped']) / sum(reads.all.lines[, 'PV13.mapped']) # 0.48
norm.23 <- sum(reads.all.lines[, 'PV21.mapped']) / sum(reads.all.lines[, 'PV23.mapped']) # 0.93
norm.21 <- 1

# coverage differences (with normalised read counts)

reads.all.lines$cov.13v21 <- log2(((reads.all.lines$PV13.mapped + 1)/(reads.all.lines$PV21.mapped + 1))/median(cov.13v21_norm))
reads.all.lines$cov.13v23 <- log2(((reads.all.lines$PV13.mapped + 1)/(reads.all.lines$PV23.mapped + 1))/median(cov.13v23_norm))
reads.all.lines$cov.04v21 <- log2(((reads.all.lines$PV04.mapped + 1)/(reads.all.lines$PV21.mapped + 1))/median(cov.04v21_norm))
reads.all.lines$cov.04v23 <- log2(((reads.all.lines$PV04.mapped + 1)/(reads.all.lines$PV23.mapped + 1))/median(cov.04v23_norm))
reads.all.lines$cov.04v13 <- log2(((reads.all.lines$PV04.mapped + 1)/(reads.all.lines$PV13.mapped + 1))/median(cov.04v13_norm))
reads.all.lines$cov.21v23 <- log2(((reads.all.lines$PV21.mapped + 1)/(reads.all.lines$PV23.mapped + 1))/median(cov.21v23_norm))

# plot again

p1 <- ggplot(reads.all.lines, aes(cov.13v21)) + geom_histogram(bins=150) + theme_bw() + labs(title="PV13 v PV21",x="log2(norm read count + 1) ratio", y="Scaffold count")
p2 <- ggplot(reads.all.lines, aes(cov.13v23)) + geom_histogram(bins=150) + theme_bw() + labs(title="PV13 v PV23",x="log2(norm read count + 1) ratio", y="Scaffold count")
p3 <- ggplot(reads.all.lines, aes(cov.04v21)) + geom_histogram(bins=150) + theme_bw() + labs(title="PV04 v PV21",x="log2(norm read count + 1) ratio", y="Scaffold count")
p4 <- ggplot(reads.all.lines, aes(cov.04v23)) + geom_histogram(bins=150) + theme_bw() + labs(title="PV04 v PV23",x="log2(norm read count + 1) ratio", y="Scaffold count")
p5 <- ggplot(reads.all.lines, aes(cov.04v13)) + geom_histogram(bins=150) + theme_bw() + labs(title="PV04 v PV13",x="log2(norm read count + 1) ratio", y="Scaffold count")
p6 <- ggplot(reads.all.lines, aes(cov.21v23)) + geom_histogram(bins=150) + theme_bw() + labs(title="PV21 v PV23",x="log2(norm read count + 1) ratio", y="Scaffold count")

hist3 <- p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 2)

# let's look at B candidates

reads.all.lines$b.status <- ifelse((reads.all.lines$cov.13v21 >= 0.58 & reads.all.lines$cov.13v23 >= 0.58 & reads.all.lines$cov.04v21 >= 0.58 & reads.all.lines$cov.04v23 >= 0.58), "B.loose", "A")
reads.all.lines$b.status <- ifelse((reads.all.lines$cov.13v21 >= 2 & reads.all.lines$cov.13v23 >= 2 & reads.all.lines$cov.04v21 >= 2 & reads.all.lines$cov.04v23 >= 2), "B.strict", reads.all.lines$b.status)

table(reads.all.lines$b.status)
sum(reads.all.lines[reads.all.lines$b.status == "B.strict",]$length)
sum(reads.all.lines[reads.all.lines$b.status != "A",]$length)

# let's look at B candidates

reads.B.lines <- reads.all.lines[c(1,2,3,4,13)]

norm2.13 <- sum(reads.all.lines[, 'PV04.mapped']) / sum(reads.all.lines[, 'PV13.mapped']) # 0.97
norm2.04 <- 1

reads.B.lines$PV13.read.cov <- reads.B.lines$PV13.mapped*norm2.13/reads.B.lines$length
reads.B.lines$PV04.read.cov <- reads.B.lines$PV04.mapped*norm2.04/reads.B.lines$length

p1 <- ggplot(reads.B.lines, aes(log10(PV13.read.cov+1e-4),log10(PV04.read.cov+1e-4))) + geom_point(aes(colour=b.status),size=1)  + scale_color_manual(values=c("azure3", "darkgreen", "deeppink3")) +
  labs(title="log10(norm read cov + 1e-4)", y="PV04", x = "PV13") + theme_bw()

reads.B.lines.cov <- reads.B.lines[c(1,5,6,7)]
colnames(reads.B.lines.cov)[3] <- "PV13"
colnames(reads.B.lines.cov)[4] <- "PV04"

reads.B.lines.long <- melt(reads.B.lines.cov, id.vars=c("seq","b.status"))
colnames(reads.B.lines.long)[3] <- "B.line"
colnames(reads.B.lines.long)[4] <- "read.cov"

p2 <- ggplot(reads.B.lines.long, aes(B.line, log10(read.cov+1e-4),fill=b.status)) +
  geom_boxplot(alpha=0.75,outlier.shape = NA,notch=TRUE,lwd=0.6) + #ylim(-1.5,2) +
  scale_fill_manual(breaks = c("A","B.loose","B.strict"), values = c("azure3", "darkgreen", "deeppink3")) + 
  theme_bw() 

p.depth <- p1 + p2

# coverage differences between PV04 and PV13:

aggregate((reads.B.lines$PV04.read.cov+1e-4)/(reads.B.lines$PV13.read.cov+1e-4)~b.status, FUN=mean, data = reads.B.lines)
aggregate((reads.B.lines$PV04.read.cov+1e-4)/(reads.B.lines$PV13.read.cov+1e-4)~b.status, FUN=sd, data = reads.B.lines)

# inspect B strict set

B.strict <- reads.all.lines[reads.all.lines$b.status == "B.strict",]
ggplot(B.strict, aes(length)) + geom_bar() + scale_x_binned(n.breaks = 20, limits = c(1,200000)) + labs(x="Length", y="Scaffold count") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# import hits

freeze.v0.blast.uniprot <- read_delim("p.viburni.freeze.v0.braker.aa.blast.vs.uniprot.out", 
                                      "\t", escape_double = FALSE, col_names = FALSE, 
                                      trim_ws = TRUE)

freeze.v0.blast.uniprot <- freeze.v0.blast.uniprot[c(1,2)]
colnames(freeze.v0.blast.uniprot)[1] <- "transcript"
colnames(freeze.v0.blast.uniprot)[2] <- "uniprot"

freeze.v0.diamond.refprot <- read_delim("p.viburni.freeze.v0.braker.aa.diamond.vs.refprot.out", 
                                        "\t", escape_double = FALSE, col_names = FALSE, 
                                        trim_ws = TRUE)

freeze.v0.diamond.refprot <- freeze.v0.diamond.refprot[c(1,2)]
colnames(freeze.v0.diamond.refprot)[1] <- "transcript"
colnames(freeze.v0.diamond.refprot)[2] <- "refprot"

# keep first (best) hit only
freeze.v0.blast.uniprot.best <- freeze.v0.blast.uniprot %>% distinct(transcript, .keep_all = TRUE)
freeze.v0.diamond.refprot.best <- freeze.v0.diamond.refprot %>% distinct(transcript, .keep_all = TRUE)

# merge with candidates
freeze.v0.genes <- merge(freeze.v0.blast.uniprot.best,freeze.v0.diamond.refprot.best,by="transcript",all=TRUE)
freeze.v0.genes <- merge(freeze.v0.genes,freeze.v0.blast.transcriptome,by="transcript",all=TRUE)
freeze.v0.genes$gene <- gsub("\\..*","", freeze.v0.genes$transcript)

genes.in.B.strict.candidates <- read_delim("genes.in.B.strict.candidates.txt", 
                                           "\t", escape_double = FALSE, col_names = FALSE, 
                                           trim_ws = TRUE)

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

# export lists
# write.csv(B.strict,"/Users/agarcia/Documents/genomics/B_viburni_ross_lab/data/coverage_analysis/B.strict.primary.reads.mapped.no.mismatches.csv", row.names = FALSE)
# write.csv(genes.in.B.strict.anno,"/Users/agarcia/Documents/genomics/B_viburni_ross_lab/data/coverage_analysis/genes.in.B.strict.anno.csv", row.names = FALSE)

# get alignments of SPAdes assemblies to the Pacbio reference

# import files: many-to-many
pviburni.freeze <- read_delim("p.viburni.freeze.v0.softmasked.fa.fai", 
                              "\t", escape_double = FALSE, col_names = FALSE, 
                              trim_ws = TRUE)
B_strict.candidates <- read_csv("B.strict.candidates", 
                                col_names = FALSE)
B_strict.candidates$B.candidate <- "Y"

colnames(pviburni.freeze)[1] <- "scaffold"
colnames(pviburni.freeze)[2] <- "len"

B_strict.candidates$B.candidate <- "Y"
colnames(B_strict.candidates)[1] <- "scaffold"

spades.nucmer.04 <- read_csv("spades/04.spades.v.freeze.v0.dnadiff.mcoords.list", col_names = FALSE)
spades.nucmer.13 <- read_csv("spades/13.spades.v.freeze.v0.dnadiff.mcoords.list", col_names = FALSE)
spades.nucmer.21 <- read_csv("spades/21.spades.v.freeze.v0.dnadiff.mcoords.list", col_names = FALSE)
spades.nucmer.23 <- read_csv("spades/23.spades.v.freeze.v0.dnadiff.mcoords.list", col_names = FALSE)

# create

pviburni.scaffold.table <- left_join(pviburni.freeze[c(1,2)],B_strict.candidates,by="scaffold")

colnames(spades.nucmer.04)[1] <- "scaffold"
colnames(spades.nucmer.13)[1] <- "scaffold"
colnames(spades.nucmer.21)[1] <- "scaffold"
colnames(spades.nucmer.23)[1] <- "scaffold"
spades.nucmer.04$PV.04 <- "Y"
spades.nucmer.13$PV.13 <- "Y"
spades.nucmer.21$PV.21 <- "Y"
spades.nucmer.23$PV.23 <- "Y"

pviburni.scaffold.table <- left_join(pviburni.scaffold.table,spades.nucmer.04,by="scaffold")
pviburni.scaffold.table <- left_join(pviburni.scaffold.table,spades.nucmer.13,by="scaffold")
pviburni.scaffold.table <- left_join(pviburni.scaffold.table,spades.nucmer.21,by="scaffold")
pviburni.scaffold.table <- left_join(pviburni.scaffold.table,spades.nucmer.23,by="scaffold")

# how many scaffolds are in the B+ lines and not in the B- lines?

pviburni.scaffold.in.Bplus.m <- pviburni.scaffold.table[!is.na(pviburni.scaffold.table$PV.13) & !is.na(pviburni.scaffold.table$PV.04)
                                                        & is.na(pviburni.scaffold.table$PV.21) & is.na(pviburni.scaffold.table$PV.23)
                                                        & !is.na(pviburni.scaffold.table$B.candidate),]

count(pviburni.scaffold.in.Bplus.m$B.candidate)
sum(pviburni.scaffold.in.Bplus.m$len)

# reimport files: 1-to-1

spades.nucmer.04 <- read_csv("spades/04.spades.v.freeze.v0.dnadiff.1coords.list", col_names = FALSE)
spades.nucmer.13 <- read_csv("spades/13.spades.v.freeze.v0.dnadiff.1coords.list", col_names = FALSE)
spades.nucmer.21 <- read_csv("spades/21.spades.v.freeze.v0.dnadiff.1coords.list", col_names = FALSE)
spades.nucmer.23 <- read_csv("spades/23.spades.v.freeze.v0.dnadiff.1coords.list", col_names = FALSE)

# create

pviburni.scaffold.table <- left_join(pviburni.freeze[c(1,2)],B_strict.candidates,by="scaffold")

colnames(spades.nucmer.04)[1] <- "scaffold"
colnames(spades.nucmer.13)[1] <- "scaffold"
colnames(spades.nucmer.21)[1] <- "scaffold"
colnames(spades.nucmer.23)[1] <- "scaffold"
spades.nucmer.04$PV.04 <- "Y"
spades.nucmer.13$PV.13 <- "Y"
spades.nucmer.21$PV.21 <- "Y"
spades.nucmer.23$PV.23 <- "Y"

nrow(spades.nucmer.04)
nrow(spades.nucmer.13)
nrow(spades.nucmer.21)
nrow(spades.nucmer.23)

pviburni.scaffold.table <- left_join(pviburni.scaffold.table,spades.nucmer.04,by="scaffold")
pviburni.scaffold.table <- left_join(pviburni.scaffold.table,spades.nucmer.13,by="scaffold")
pviburni.scaffold.table <- left_join(pviburni.scaffold.table,spades.nucmer.21,by="scaffold")
pviburni.scaffold.table <- left_join(pviburni.scaffold.table,spades.nucmer.23,by="scaffold")

# how many scaffolds are in the B+ lines and not in the B- lines?

pviburni.B.scaffolds <- pviburni.scaffold.table[!is.na(pviburni.scaffold.table$PV.13) & !is.na(pviburni.scaffold.table$PV.04)
                                                & is.na(pviburni.scaffold.table$PV.21) & is.na(pviburni.scaffold.table$PV.23),]
count(pviburni.B.scaffolds$B.candidate)

# let's plot this

reads.B.lines.spades <- reads.B.lines
reads.B.lines.spades$b.status <- ifelse(reads.B.lines.spades$seq %in% pviburni.B.scaffolds$scaffold & reads.B.lines.spades$b.status == "B.strict", "B.strict.plus.assembly", reads.B.lines.spades$b.status)
reads.B.lines.spades$b.status <- ifelse(reads.B.lines.spades$seq %in% pviburni.B.scaffolds$scaffold & reads.B.lines.spades$b.status == "B.loose", "B.loose.plus.assembly", reads.B.lines.spades$b.status)
reads.B.lines.spades$b.status <- ifelse(reads.B.lines.spades$seq %in% pviburni.B.scaffolds$scaffold & reads.B.lines.spades$b.status == "A", "B.assembly", reads.B.lines.spades$b.status)
count(reads.B.lines.spades$b.status)
ddply(reads.B.lines.spades,c("b.status"),summarise, N = length(seq), size = sum(length)/1000000) # get counts

reads.B.lines.spades$b.status <- factor(reads.B.lines.spades$b.status, levels = c("B.strict.plus.assembly","B.strict","B.loose.plus.assembly","B.loose","B.assembly","A"))

p1 <- ggplot(reads.B.lines.spades, aes(log10(PV13.read.cov+1e-4),log10(PV04.read.cov+1e-4))) + geom_point(aes(colour=b.status),size=1) +
  scale_color_manual(values=c("royalblue4", "dodgerblue", "green4", "green1", "lavenderblush4", "lavenderblush1")) +
  labs(title="log10(norm read cov + 1e-4)", y="PV04", x = "PV13") + theme_bw()

reads.B.lines.spades.cov <- reads.B.lines.spades[c(1,5,6,7)]
colnames(reads.B.lines.spades.cov)[3] <- "PV13"
colnames(reads.B.lines.spades.cov)[4] <- "PV04"

reads.B.lines.spades.long <- melt(reads.B.lines.spades.cov, id.vars=c("seq","b.status"))
colnames(reads.B.lines.spades.long)[3] <- "B.line"
colnames(reads.B.lines.spades.long)[4] <- "read.cov"

p2 <- ggplot(reads.B.lines.spades.long, aes(B.line, log10(read.cov+1e-4),fill=b.status)) +
  geom_boxplot(alpha=0.75,outlier.shape = NA,notch=FALSE,lwd=0.6) + ylim(-1.5,2) +
  scale_fill_manual(values=c("royalblue4", "dodgerblue", "green4", "green1", "lavenderblush4", "lavenderblush1")) + 
  theme_bw() 
p1 + p2
aggregate((reads.B.lines.spades$PV04.read.cov+1e-4)/(reads.B.lines.spades$PV13.read.cov+1e-4)~b.status, FUN=mean, data = reads.B.lines.spades)
aggregate((reads.B.lines.spades$PV04.read.cov+1e-4)/(reads.B.lines.spades$PV13.read.cov+1e-4)~b.status, FUN=sd, data = reads.B.lines.spades)
