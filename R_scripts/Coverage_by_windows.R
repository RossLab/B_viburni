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

setwd("/Users/agarcia/Documents/genomics/B_viburni_ross_lab/data/coverage_analysis/windows")

##### Import and organise files

# chromosome assignment
scaffolds_final_assignment <- read_csv("~/Documents/genomics/B_viburni_ross_lab/output/scaffolds.final.assignment.csv")

# coverage by 1kb windows
PV_18_04_coverage_per_1kb <- read_delim("PV_18-04.coverage.per.1kb.window", "\t", escape_double = FALSE, trim_ws = TRUE)
PV_18_13_coverage_per_1kb <- read_delim("PV_18-13.coverage.per.1kb.window", "\t", escape_double = FALSE, trim_ws = TRUE)
PV_18_21_coverage_per_1kb <- read_delim("PV_18-21.coverage.per.1kb.window", "\t", escape_double = FALSE, trim_ws = TRUE)
PV_18_23_coverage_per_1kb <- read_delim("PV_18-23.coverage.per.1kb.window", "\t", escape_double = FALSE, trim_ws = TRUE)

# coverage by 1kb windows

seq <- PV_18_04_coverage_per_1kb[c(1)]
pos <- PV_18_04_coverage_per_1kb[c(3)]
PV04.cov <- PV_18_04_coverage_per_1kb[c(5)]
PV13.cov <- PV_18_13_coverage_per_1kb[c(5)]
PV21.cov <- PV_18_21_coverage_per_1kb[c(5)]
PV23.cov <- PV_18_23_coverage_per_1kb[c(5)]

windows.cov <- data.frame(seq, pos,PV04.cov,PV13.cov,PV21.cov,PV23.cov)
colnames(windows.cov) <- c("seq","window","PV04","PV13","PV21","PV23")

# estimate average coverage in B+/B- lines normalise by median coverage differences

windows.cov$B.avg <- rowMeans(windows.cov[3:4])
windows.cov$nonB.avg <- rowMeans(windows.cov[5:6])
cov.BvBminus_norm <- (windows.cov$B.avg+1)/(windows.cov$nonB.avg+1)
median(cov.BvBminus_norm) 
windows.cov$ratio.B <- log2((windows.cov$B.avg/median(cov.BvBminus_norm) + 1)/(windows.cov$nonB.avg + 1))

windows.cov$B.avg <- rowMeans(windows.cov[3:4])
windows.cov$nonB.avg <- rowMeans(windows.cov[5:6])
cov.04v13_norm <- (windows.cov$PV04+1)/(windows.cov$PV13+1)
median(cov.04v13_norm)
windows.cov$ratio.04v13 <- log2((windows.cov$PV04/median(cov.04v13_norm) + 1)/(windows.cov$PV13 + 1))
View(windows.cov)

windows.cov <- merge(windows.cov,scaffolds_final_assignment[c(1,27)],by=c("seq"))
windows.cov.A  <- windows.cov[windows.cov$"b.status.final" == "A",]
windows.cov.B1 <- windows.cov[windows.cov$"b.status.final" == "B1",]
windows.cov.B2 <- windows.cov[windows.cov$"b.status.final" == "B2",]
windows.cov.B3 <- windows.cov[windows.cov$"b.status.final" == "B3",]
windows.cov.B4 <- windows.cov[windows.cov$"b.status.final" == "B4",]

pdf("B1.scaffold.by.window.pdf", width=10, height=10, pointsize=12)
par(mfrow = c(1, 1)) # 2 rows, 2 columns
for (i in unique(windows.cov.B1$seq)) {
  a <- ggplot(windows.cov.B1[windows.cov.B1$seq == i,]) +
    geom_point(aes(window, ratio.B),colour="black",size=0.2) +
    geom_point(aes(window, ratio.04v13),colour="tomato",size=0.2) +
    ylim(c(-10,10)) +
    geom_hline(yintercept=1, linetype="dashed", color = "steelblue", size=0.3) +
    geom_hline(yintercept=0, linetype="dashed", color = "darkblue", size=0.3) +
    labs(title=i, y="log2 ratios", x = "Position") + theme_classic()
  print(a)
}
dev.off()


pdf("B2.scaffold.by.window.pdf", width=10, height=10, pointsize=12)
par(mfrow = c(1, 1)) # 2 rows, 2 columns
for (i in unique(windows.cov.B2$seq)) {
  a <- ggplot(windows.cov.B2[windows.cov.B2$seq == i,]) +
    geom_point(aes(window, ratio.B),colour="black",size=0.2) +
    geom_point(aes(window, ratio.04v13),colour="tomato",size=0.2) +
    ylim(c(-10,10)) +
    geom_hline(yintercept=1, linetype="dashed", color = "steelblue", size=0.3) +
    geom_hline(yintercept=0, linetype="dashed", color = "darkblue", size=0.3) +
    labs(title=i, y="log2 ratios", x = "Position") + theme_classic()
  print(a)
}
dev.off()



#pdf("rsem_gene_mdplot-prenorm.pdf", width=10, height=10, pointsize=12)
#par(mfrow = c(1, 1)) # 2 rows, 2 columns
for (i in 1:ncol(x)) {
  plotMD(cpm(x, log=TRUE), column=i, xlab = "Average log-expression", ylab = "Expression log-ratio (this sample vs others)", main = colnames(x)[i])
  abline(h=0, col="red", lty=2, lwd=2)
}
#dev.off()

ggplot(windows.cov[windows.cov$seq == "scaffold_1520",]) +
  geom_point(aes(window, ratio.B),colour="black",size=0.1) +
  geom_point(aes(window, ratio.04v13),colour="tomato",size=0.1) + ylim(c(-10,10)) +
  geom_hline(yintercept=1, linetype="dashed",  color = "darkblue", size=0.3) +
  geom_hline(yintercept=0, linetype="dashed",  color = "steelblue", size=0.3) +
  theme_classic()

View(windows.cov[windows.cov$seq == "scaffold_1520",])
