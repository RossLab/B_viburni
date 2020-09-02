rm(list=ls())
ls()
library(tidyverse)
library(gridExtra)
#devtools::install_github("thomasp85/patchwork")
library(patchwork)
library(lattice)
library(grid)
library(gridExtra)
library(reshape2)
library(plyr)

setwd("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/data/coverage_analysis")

# import number of reads mapped to scaffold

PV04.reads.mapped <- read_delim("PV_18-04.primary.reads.mapped.count", 
                                    "\t", escape_double = FALSE, col_names = FALSE, 
                                    trim_ws = TRUE)
PV13.reads.mapped <- read_delim("PV_18-13.primary.reads.mapped.count", 
                                "\t", escape_double = FALSE, col_names = FALSE, 
                                trim_ws = TRUE)
PV21.reads.mapped <- read_delim("PV_18-21.primary.reads.mapped.count", 
                                "\t", escape_double = FALSE, col_names = FALSE, 
                                trim_ws = TRUE)
PV23.reads.mapped <- read_delim("PV_18-23.primary.reads.mapped.count", 
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

sum(reads.all.lines[, 'PV04.mapped']) # 388044763
sum(reads.all.lines[, 'PV13.mapped']) # 363059148
sum(reads.all.lines[, 'PV21.mapped']) # 197844424
sum(reads.all.lines[, 'PV23.mapped']) # 208740794

# normalisation factor

norm.04 <- sum(reads.all.lines[, 'PV21.mapped']) / sum(reads.all.lines[, 'PV04.mapped']) # 0.51
norm.13 <- sum(reads.all.lines[, 'PV21.mapped']) / sum(reads.all.lines[, 'PV13.mapped']) # 0.55
norm.23 <- sum(reads.all.lines[, 'PV21.mapped']) / sum(reads.all.lines[, 'PV23.mapped']) # 0.95
norm.21 <- 1

# coverage differences (with normalised read counts)

reads.all.lines$cov.13v21 <- log2(((reads.all.lines$PV13.mapped)*norm.13 + 1)/((reads.all.lines$PV21.mapped)*norm.21 + 1))
reads.all.lines$cov.13v23 <- log2(((reads.all.lines$PV13.mapped)*norm.13 + 1)/((reads.all.lines$PV23.mapped)*norm.23 + 1))
reads.all.lines$cov.04v21 <- log2(((reads.all.lines$PV04.mapped)*norm.04 + 1)/((reads.all.lines$PV21.mapped)*norm.21 + 1))
reads.all.lines$cov.04v23 <- log2(((reads.all.lines$PV04.mapped)*norm.04 + 1)/((reads.all.lines$PV23.mapped)*norm.23 + 1))
reads.all.lines$cov.04v13 <- log2(((reads.all.lines$PV04.mapped)*norm.04 + 1)/((reads.all.lines$PV13.mapped)*norm.13 + 1))
reads.all.lines$cov.21v23 <- log2(((reads.all.lines$PV21.mapped)*norm.21 + 1)/((reads.all.lines$PV23.mapped)*norm.23 + 1))

# plot again

p1 <- ggplot(reads.all.lines, aes(cov.13v21)) + geom_histogram(bins=150) + theme_bw() + labs(title="PV13 v PV21",x="log2(norm read count + 1) ratio", y="Scaffold count")
p2 <- ggplot(reads.all.lines, aes(cov.13v23)) + geom_histogram(bins=150) + theme_bw() + labs(title="PV13 v PV23",x="log2(norm read count + 1) ratio", y="Scaffold count")
p3 <- ggplot(reads.all.lines, aes(cov.04v21)) + geom_histogram(bins=150) + theme_bw() + labs(title="PV04 v PV21",x="log2(norm read count + 1) ratio", y="Scaffold count")
p4 <- ggplot(reads.all.lines, aes(cov.04v23)) + geom_histogram(bins=150) + theme_bw() + labs(title="PV04 v PV23",x="log2(norm read count + 1) ratio", y="Scaffold count")
p5 <- ggplot(reads.all.lines, aes(cov.04v13)) + geom_histogram(bins=150) + theme_bw() + labs(title="PV04 v PV13",x="log2(norm read count + 1) ratio", y="Scaffold count")
p6 <- ggplot(reads.all.lines, aes(cov.21v23)) + geom_histogram(bins=150) + theme_bw() + labs(title="PV21 v PV23",x="log2(norm read count + 1) ratio", y="Scaffold count")

hist3 <- p1 + p2 + p3 + p4 + p5 + p6

# let's look at B candidates

reads.all.lines$b.status <- ifelse((reads.all.lines$cov.13v21 >= 0.58 & reads.all.lines$cov.13v23 >= 0.58 & reads.all.lines$cov.04v21 >= 0.58 & reads.all.lines$cov.04v23 >= 0.58), "B.loose", "A")
reads.all.lines$b.status <- ifelse((reads.all.lines$cov.13v21 >= 2 & reads.all.lines$cov.13v23 >= 2 & reads.all.lines$cov.04v21 >= 2 & reads.all.lines$cov.04v23 >= 2), "B.strict", reads.all.lines$b.status)

table(reads.all.lines$b.status)
sum(reads.all.lines[reads.all.lines$b.status != "A",]$length)

# let's look at B candidates

reads.B.lines <- reads.all.lines[c(1,2,3,4,13)]

norm2.04 <- sum(reads.all.lines[, 'PV13.mapped']) / sum(reads.all.lines[, 'PV04.mapped']) # 0.94
norm2.13 <- 1

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



