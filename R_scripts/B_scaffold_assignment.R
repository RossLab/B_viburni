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

# setwd("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/data/coverage_analysis")

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

# normalisation factor: median coverage difference between pairs of lines

cov.13v21_norm <- round((PV13.reads.mapped$PV13.mapped+1)/(PV21.reads.mapped$PV21.mapped+1),2)
cov.13v23_norm <- round((PV13.reads.mapped$PV13.mapped+1)/(PV23.reads.mapped$PV23.mapped+1),2)
cov.04v21_norm <- round((PV04.reads.mapped$PV04.mapped+1)/(PV21.reads.mapped$PV21.mapped+1),2)
cov.04v23_norm <- round((PV04.reads.mapped$PV04.mapped+1)/(PV23.reads.mapped$PV23.mapped+1),2)
cov.04v13_norm <- round((PV04.reads.mapped$PV04.mapped+1)/(PV13.reads.mapped$PV13.mapped+1),2)
cov.21v23_norm <- round((PV21.reads.mapped$PV21.mapped+1)/(PV23.reads.mapped$PV23.mapped+1),2)
median(cov.13v21_norm) # 2.08
median(cov.13v23_norm) # 1.91
median(cov.04v21_norm) # 1.87
median(cov.04v23_norm) # 1.77
median(cov.04v13_norm) # 0.95
median(cov.21v23_norm) # 0.94

# normalised coverage differences

reads.all.lines$cov.13v21 <- log2((reads.all.lines$PV13.mapped/median(cov.13v21_norm) + 1)/(reads.all.lines$PV21.mapped + 1))
reads.all.lines$cov.13v23 <- log2((reads.all.lines$PV13.mapped/median(cov.13v23_norm) + 1)/(reads.all.lines$PV23.mapped + 1))
reads.all.lines$cov.04v21 <- log2((reads.all.lines$PV04.mapped/median(cov.04v21_norm) + 1)/(reads.all.lines$PV21.mapped + 1))
reads.all.lines$cov.04v23 <- log2((reads.all.lines$PV04.mapped/median(cov.04v23_norm) + 1)/(reads.all.lines$PV23.mapped + 1))
reads.all.lines$cov.04v13 <- log2((reads.all.lines$PV04.mapped/median(cov.04v13_norm) + 1)/(reads.all.lines$PV13.mapped + 1))
reads.all.lines$cov.21v23 <- log2((reads.all.lines$PV21.mapped/median(cov.21v23_norm) + 1)/(reads.all.lines$PV23.mapped + 1))

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
sum(reads.all.lines[reads.all.lines$b.status == "B.loose",]$length)

# let's look at B candidates

reads.B.lines <- reads.all.lines[c(1,2,3,4,11,13)]
reads.B.lines$PV13.read.cov <- reads.B.lines$PV13.mapped*median(cov.04v13_norm)/reads.B.lines$length
reads.B.lines$PV04.read.cov <- reads.B.lines$PV04.mapped*1/reads.B.lines$length

p1 <- ggplot(reads.B.lines, aes(log10(PV13.read.cov+1e-4),log10(PV04.read.cov+1e-4))) + geom_point(aes(colour=b.status),size=1)  + scale_color_manual(values=c("azure3", "darkgreen", "deeppink3")) +
  labs(title="log10(norm read cov + 1e-4)", y="PV04", x = "PV13") + theme_bw()

reads.B.lines.cov <- reads.B.lines[c(1,6,7,8)]
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
nrow(B.strict)
B.all <- reads.all.lines[reads.all.lines$b.status != "A",]
nrow(B.all)

ggplot(B.strict, aes(length)) + geom_bar() + scale_x_binned(n.breaks = 20, limits = c(1,200000)) + labs(x="Length", y="Scaffold count") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(B.all, aes(length)) + geom_bar() + scale_x_binned(n.breaks = 20, limits = c(1,200000)) + labs(x="Length", y="Scaffold count") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# get alignments of SPAdes assemblies to the Pacbio reference

# import files: 1-to-1
pviburni.freeze <- read_delim("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/data/annotation/p.viburni.freeze.v0.softmasked.fa.fai",
                              "\t", escape_double = FALSE, col_names = FALSE,
                              trim_ws = TRUE)

colnames(pviburni.freeze)[1] <- "scaffold"
colnames(pviburni.freeze)[2] <- "len"

spades.nucmer.04 <- read_csv("spades/04.spades.v.freeze.v0.dnadiff.1coords.list", col_names = FALSE)
spades.nucmer.13 <- read_csv("spades/13.spades.v.freeze.v0.dnadiff.1coords.list", col_names = FALSE)
spades.nucmer.21 <- read_csv("spades/21.spades.v.freeze.v0.dnadiff.1coords.list", col_names = FALSE)
spades.nucmer.23 <- read_csv("spades/23.spades.v.freeze.v0.dnadiff.1coords.list", col_names = FALSE)

colnames(spades.nucmer.04)[1] <- "seq"
colnames(spades.nucmer.13)[1] <- "seq"
colnames(spades.nucmer.21)[1] <- "seq"
colnames(spades.nucmer.23)[1] <- "seq"
spades.nucmer.04$PV04.asn <- "Y"
spades.nucmer.13$PV13.asn <- "Y"
spades.nucmer.21$PV21.asn <- "Y"
spades.nucmer.23$PV23.asn <- "Y"

reads.spades.all.lines <- left_join(reads.all.lines,spades.nucmer.04,by="seq")
reads.spades.all.lines <- left_join(reads.spades.all.lines,spades.nucmer.13,by="seq")
reads.spades.all.lines <- left_join(reads.spades.all.lines,spades.nucmer.21,by="seq")
reads.spades.all.lines <- left_join(reads.spades.all.lines,spades.nucmer.23,by="seq")

# how many scaffolds are in the B+ lines and not in the B- lines?

reads.spades.all.lines$b.status.asn <- ifelse(!is.na(reads.spades.all.lines$PV13.asn) & !is.na(reads.spades.all.lines$PV04.asn) &
                                                is.na(reads.spades.all.lines$PV21.asn) & is.na(reads.spades.all.lines$PV23.asn),"B","A")
table(reads.spades.all.lines$b.status.asn)
sum(reads.spades.all.lines[reads.spades.all.lines$b.status.asn == "B",]$len)

# compare with coverage

reads.spades.all.lines$b.status.cov.asn <- ifelse((reads.spades.all.lines$b.status == "B.strict" & reads.spades.all.lines$b.status.asn == "B"),
                                                "B.strict.plus.asn", reads.spades.all.lines$b.status)
reads.spades.all.lines$b.status.cov.asn <- ifelse((reads.spades.all.lines$b.status == "B.loose" & reads.spades.all.lines$b.status.asn == "B"),
                                                  "B.loose.plus.asn", reads.spades.all.lines$b.status.cov.asn)
reads.spades.all.lines$b.status.cov.asn <- ifelse((reads.spades.all.lines$b.status == "A" & reads.spades.all.lines$b.status.asn == "B"),
                                                  "B.asn", reads.spades.all.lines$b.status.cov.asn)

table(reads.spades.all.lines$b.status.cov.asn)
ddply(reads.spades.all.lines,c("b.status.cov.asn"),summarise, N = length(seq), size = sum(length)/1000000) # get counts

reads.spades.all.lines$b.status.cov.asn <- factor(reads.spades.all.lines$b.status.cov.asn, levels = c("B.strict.plus.asn","B.strict","B.loose.plus.asn","B.loose","B.asn","A"))

p1 <- ggplot(reads.spades.all.lines, aes(log10(PV13.mapped*median(cov.04v13_norm)/length + 1e-4),
                                         log10(PV04.mapped*1/length + 1e-4))) + geom_point(aes(colour=b.status.cov.asn),size=1) +
  scale_color_manual(values=c("royalblue4", "dodgerblue", "green4", "green1", "lavenderblush4", "lavenderblush1")) +
  labs(title="log10(norm read cov + 1e-4)", y="PV04", x = "PV13") + theme_bw()

# import kmers

a.lo.kmers <- read_delim("a.candidates.loose.mapped.kmers",
                                "\t", escape_double = FALSE, col_names = FALSE,
                                trim_ws = TRUE)
b.lo.kmers <- read_delim("b.candidates.loose.mapped.kmers",
                          "\t", escape_double = FALSE, col_names = FALSE,
                          trim_ws = TRUE)
a.st.kmers <- read_delim("a.candidates.strict.mapped.kmers",
                          "\t", escape_double = FALSE, col_names = FALSE,
                          trim_ws = TRUE)
b.st.kmers <- read_delim("b.candidates.strict.mapped.kmers",
                          "\t", escape_double = FALSE, col_names = FALSE,
                          trim_ws = TRUE)

seq <- a.lo.kmers[ ,1]
length <- a.lo.kmers[ ,2]
A.lo <- a.lo.kmers[ ,3]
B.lo <- b.lo.kmers[ ,3]
A.st <- a.st.kmers[ ,3]
B.st <- b.st.kmers[ ,3]
kmers.all.lines0 <- data.frame(seq,length,A.lo,B.lo,A.st,B.st)
kmers.all.lines <-kmers.all.lines0[1:(nrow(kmers.all.lines0)-1),] # remove last line (not a scaffold)
colnames(kmers.all.lines) <- c("seq","length","A.lo","B.lo","A.st","B.st")

kmers.all.lines$AB.ratio.st <- log2((kmers.all.lines$B.st+1)/(kmers.all.lines$A.st+1))
kmers.all.lines$AB.ratio.lo <- log2((kmers.all.lines$B.lo+1)/(kmers.all.lines$A.lo+1))

scaffolds.final.assignment <- left_join(reads.spades.all.lines,kmers.all.lines)

table(scaffolds.final.assignment[scaffolds.final.assignment$AB.ratio.lo >= 0,]$b.status.cov.asn)
table(scaffolds.final.assignment[scaffolds.final.assignment$AB.ratio.st >= 0,]$b.status.cov.asn)

ggplot(scaffolds.final.assignment, aes(AB.ratio.st,cov.04v13)) + geom_point(aes(colour=b.status.cov.asn),size=1) +
  scale_color_manual(values=c("royalblue4", "dodgerblue", "green4", "green1", "lavenderblush4", "lavenderblush1")) +
  labs(title="log10(norm read cov + 1e-4)", y="PV04", x = "PV13") + theme_bw()

scaffolds.final.assignment$b.status.kmer <- ifelse((scaffolds.final.assignment$AB.ratio.lo > 0),"B","A")
table(scaffolds.final.assignment$b.status.kmer)

# final assignment

scaffolds.final.assignment$b.status.final <- 'A'

scaffolds.final.assignment$b.status.final <- ifelse(scaffolds.final.assignment$cov.04v13 > 0 &
                                                      (scaffolds.final.assignment$b.status == "B.strict") &
                                                      (scaffolds.final.assignment$b.status.asn == "B" | scaffolds.final.assignment$b.status.kmer == "B"),
                                                      "B1", scaffolds.final.assignment$b.status.final)

scaffolds.final.assignment$b.status.final <- ifelse(scaffolds.final.assignment$b.status.final == "A" & scaffolds.final.assignment$cov.04v13 > 0 &
                                                      (scaffolds.final.assignment$b.status != "A" |
                                                         scaffolds.final.assignment$b.status.asn == "B" |
                                                         scaffolds.final.assignment$b.status.kmer == "B"),
                                                    "B2", scaffolds.final.assignment$b.status.final)

scaffolds.final.assignment$b.status.final <- ifelse(scaffolds.final.assignment$b.status.final == "A" &
                                                      (scaffolds.final.assignment$b.status != "A" |
                                                         scaffolds.final.assignment$b.status.asn == "B" |
                                                         scaffolds.final.assignment$b.status.kmer == "B"),
                                                    "B3", scaffolds.final.assignment$b.status.final)
table(scaffolds.final.assignment$b.status.final)

# collect some stats

ddply(scaffolds.final.assignment,c("b.status.final"), summarise, N = length(seq), size = sum(length)/1000000) # get counts
aggregate((PV04.mapped*1/length + 1e-4)/(PV13.mapped*median(cov.04v13_norm)/length +1e-4)~b.status.final, FUN=mean, data = scaffolds.final.assignment)
aggregate((PV04.mapped*1/length + 1e-4)/(PV13.mapped*median(cov.04v13_norm)/length +1e-4)~b.status.final, FUN=sd, data = scaffolds.final.assignment)

# plot fig2

scaffolds.cov.sum <- ddply(scaffolds.final.assignment,c("b.status"),
                           summarise, N = length(seq), size = sum(length)) # get counts
scaffolds.cov.sum$b.status <- factor(scaffolds.cov.sum$b.status,
                                     levels = c("A","B.strict","B.loose"))
scaffolds.final.assignment$b.status <- factor(scaffolds.final.assignment$b.status,
                                     levels = c("A","B.strict","B.loose"))

log.labels <- expression(0,"10"^3, "10"^4, "10"^5, "10"^6, "10"^7,"10"^8)

fig2a <- ggplot(scaffolds.cov.sum, aes(b.status, log10(size),fill=b.status)) +
  geom_bar(stat="identity", position=position_dodge(),colour="black") +
  scale_fill_manual(values=c("gray85","gray50","gray60")) +
  scale_x_discrete(labels=c("A", "Putative B\n(strict)", "Additional B\n(loose)")) +
  geom_text(aes(label=N),position = position_dodge(0.9),vjust=-1) +
  geom_jitter(data=scaffolds.final.assignment, aes(b.status, log10(length), group=b.status),
              size=0.5, width = 0.4, alpha=1, show.legend=FALSE) +
  labs(title="Scaffold assignment\nby coverage differences",x="", y ="Mb (log10)") + coord_cartesian(ylim = c(3,8.7)) +
  scale_y_continuous(breaks = c(0, 3, 4, 5, 6, 7, 8), labels = log.labels) +
  theme_classic() + theme(legend.position="none") +
  theme(legend.position="none", plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)),
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13)))

shared.num <- nrow(scaffolds.final.assignment[(!is.na(scaffolds.final.assignment$PV04.asn) & !is.na(scaffolds.final.assignment$PV13.asn) &
                                                 !is.na(scaffolds.final.assignment$PV21.asn) & !is.na(scaffolds.final.assignment$PV23.asn)),])
b.num <- nrow(scaffolds.final.assignment[scaffolds.final.assignment$b.status.asn == "B",])
none.num <- nrow(scaffolds.final.assignment[(is.na(scaffolds.final.assignment$PV04.asn) & is.na(scaffolds.final.assignment$PV13.asn) &
                                   is.na(scaffolds.final.assignment$PV21.asn) & is.na(scaffolds.final.assignment$PV23.asn)),])
other.num <- nrow(scaffolds.final.assignment)-shared.num-b.num-none.num

mappings.plot <- data.frame(assembly=c("PV04","PV13","PV21","PV23","Reference","Reference","Reference","Reference"),
           N=c(nrow(spades.nucmer.04),nrow(spades.nucmer.13),nrow(spades.nucmer.21),nrow(spades.nucmer.13),
               shared.num, b.num, none.num, other.num),
           status=c(NA,NA,NA,NA,"Shared by all","In B lines only","In none","At least one"))

mappings.plot$assembly <- factor(mappings.plot$assembly,
                                 levels = c("PV04","PV13","PV21","PV23","Reference"))
mappings.plot$status <- factor(mappings.plot$status,
                                 levels = c(NA,"In none","At least one","In B lines only","Shared by all"))
fig2b <- ggplot(mappings.plot, aes(fill=status, y=N, x=assembly)) +
  geom_bar(stat="identity",colour="black") + coord_cartesian(ylim=c(0,2400)) +
  labs(title="Scaffold assignment\nby 1-to-1 mappings",x="", y ="Scaffolds") +
  scale_fill_manual(name="Mappings\nto reference scaffolds",values = c("gray75", "gray65", "gray50", "gray85"),
                    breaks=c("In none","At least one","In B lines only","Shared by all"),
                    labels=c("In none","At least one","In B lines only","Shared by all")) +
  scale_y_continuous(breaks = c(0, 500, 1000, 1500, 2000, 2392)) + coord_flip() +
  geom_text(aes(label = N), position = position_stack(vjust = 0.5)) +
  theme_classic() + theme(legend.position="bottom") +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)),
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13)))

fig2c <- ggplot(scaffolds.final.assignment, aes(color=b.status.kmer, y=log10(A.lo+B.lo), x=AB.ratio.lo)) +
  geom_point() + labs(title="Scaffold assignment\nby kmer alignment",
                      x="log2(candidate B/A kmers)", y ="Total number of mapped\n27-mers (log10)") +
  scale_colour_manual(values = c("gray85", "gray50")) +
  geom_vline(xintercept = 0, col = 2, lty = 2) +
  theme_classic() + theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)),
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)))

scaffolds.final.assignment$b.status.final <- factor(scaffolds.final.assignment$b.status.final,
                                 levels = c("B1","B2","B3","B4","A"))

ggplot(scaffolds.final.assignment, aes(log10(PV13.mapped*median(cov.04v13_norm)/length + 1e-4),
                                   log10(PV04.mapped*1/length + 1e-4))) +
  geom_point(aes(colour=b.status.final,size=log10(length)),alpha=0.75) +
  scale_size(range=c(0.5,2)) +
  scale_color_manual(values=c("royalblue4", "green", "royalblue2", "lavenderblush4","gray85")) +
  labs(title="B chromosome assignment", y="Cov in PV04 (log10)", x = "Cov in PV13 (log10)") +
  theme_classic() + theme(legend.position="right") +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)),
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13)))

scaffolds.final.assignment.b <- scaffolds.final.assignment[scaffolds.final.assignment$b.status.final != "A",]
scaffolds.final.assignment.a <- scaffolds.final.assignment[scaffolds.final.assignment$b.status.final == "A",]

fig2d <- ggplot(scaffolds.final.assignment.b) +
  geom_point(data=scaffolds.final.assignment.a,
             aes(log10(PV13.mapped*median(cov.04v13_norm)/length + 1e-4),
                 log10(PV04.mapped*1/length + 1e-4),
                 size=log10(length)),colour="gray90", show.legend=FALSE) +
  geom_point(aes(log10(PV13.mapped*median(cov.04v13_norm)/length + 1e-4),
                 log10(PV04.mapped*1/length + 1e-4),colour=b.status.final,size=log10(length))) +
  scale_size(range=c(0.5,2)) +
  scale_color_manual(values=c("royalblue4", "deepskyblue", "cadetblue", "lavenderblush4","gray85")) +
  labs(title="Combined B chromosome assignment", y="Cov in PV04 (log10)", x = "Cov in PV13 (log10)", color = "B set") +
  theme_classic() + theme(legend.position="right") +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)),
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) + guides(size = FALSE)
View(scaffolds.final.assignment)

#write.table(scaffolds.final.assignment, file = "output/scaffolds.preprint.assignment.csv",row.names = F,sep = ",")

library(patchwork)
jpeg("/Users/agarcia/Desktop/b.assignment.final.jpeg",
     width = 4200, height = 3800, units = 'px', res = 300)
fig2a + fig2b + fig2c + fig2d
dev.off()
