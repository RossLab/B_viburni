rm(list=ls())
ls()
library(tidyverse)
library(gridExtra)
#devtools::install_github("thomasp85/patchwork")
library(patchwork)
library(lattice)
library(grid)
library(gridExtra)

setwd("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/data/coverage_analysis")

# import number of reads mapped to scaffold

PV04.reads.mapped <- read_delim("PV_18-04.reads.mapped.count", 
                                    "\t", escape_double = FALSE, col_names = FALSE, 
                                    trim_ws = TRUE)
PV13.reads.mapped <- read_delim("PV_18-13.reads.mapped.count", 
                                "\t", escape_double = FALSE, col_names = FALSE, 
                                trim_ws = TRUE)
PV21.reads.mapped <- read_delim("PV_18-21.reads.mapped.count", 
                                "\t", escape_double = FALSE, col_names = FALSE, 
                                trim_ws = TRUE)
PV23.reads.mapped <- read_delim("PV_18-23.reads.mapped.count", 
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

sum(reads.all.lines[, 'PV04.mapped']) # 395178560
sum(reads.all.lines[, 'PV13.mapped']) # 367442808
sum(reads.all.lines[, 'PV21.mapped']) # 201661407
sum(reads.all.lines[, 'PV23.mapped']) # 212848740

# normalisation factor

norm.04 <- sum(reads.all.lines[, 'PV21.mapped']) / sum(reads.all.lines[, 'PV04.mapped']) # 0.51
norm.13 <- sum(reads.all.lines[, 'PV21.mapped']) / sum(reads.all.lines[, 'PV13.mapped']) # 0.59
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

reads.all.lines$b.status <- ifelse((reads.all.lines$cov.13v21 >= 0.58 & reads.all.lines$cov.13v23 >= 0.58 & reads.all.lines$cov.04v21 >= 0.58 & reads.all.lines$cov.04v23 >= 0.58), "b.loose", "no")
reads.all.lines$b.status <- ifelse((reads.all.lines$cov.13v21 >= 2 & reads.all.lines$cov.13v23 >= 2 & reads.all.lines$cov.04v21 >= 2 & reads.all.lines$cov.04v23 >= 2), "b.strict", reads.all.lines$b.status)

table(reads.all.lines$b.status)
sum(b.candidates$length)

nrow(b.candidates.strict)
sum(b.candidates.strict$length)

# let's look at B candidates

reads.B.lines <- reads.all.lines[c(1,2,3,4)]
reads.B.lines$PV13.read.cov <- reads.B.lines$PV13.mapped/reads.B.lines$length
reads.B.lines$PV04.read.cov <- reads.B.lines$PV04.mapped/reads.B.lines$length

ggplot(reads.B.lines, aes(log10(PV13.read.cov),log10(PV04.read.cov))) + geom_point()
                            
          
  geom_bar(stat="identity", position=position_dodge(),fill="plum4") +
  scale_x_discrete(limit = c("maternal.only", "maternal.bias", "biparental","paternal.bias","paternal.only"),
                   labels = c("M","MB","B","PB","P")) +
  geom_text(aes(label=N),position = position_dodge(0.9),vjust=-1) +
  labs(title="WC", y="%", x = "Category of bias") +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = (14)), 
        axis.title = element_text(family = "Helvetica", size = (13)),
        axis.text = element_text(family = "Helvetica", size = (12)),
        legend.text = element_text(family = "Helvetica", size = (12)),
        legend.title = element_text(family = "Helvetica", size = (13))) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
  theme()




