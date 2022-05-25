assignments <- read.csv('output/scaffolds.final.assignment.csv')
row.names(assignments) <- assignments$seq
head(assignments)

# removing all the old info about mapping that is obsolete now
assignments <- assignments[, c(1, 2, 13:27)]
colnames(assignments)[3] <- 'old.cov.b.status'

coverage_table <- read.table('data/4_cov_analysis/Comeplete_window_coverage_table.tsv', header = T)
scaffolds <- unique(coverage_table$scf)
scf2fraction_of_candidates <- function(scf){
	mean(coverage_table[coverage_table$scf == scf, 'B_candidate'])
}

scaffold_B_candidate_fraction <- sapply(scaffolds, scf2fraction_of_candidates)
# hist(scaffold_B_candidate_fraction)
B_scaffolds <- scaffolds[scaffold_B_candidate_fraction > 0.6]

assignments[scaffolds, 'fraction_of_B_windows'] <- round(scaffold_B_candidate_fraction, 3)
assignments[, 'new.cov.b.status'] <- 'A'
assignments[B_scaffolds, 'new.cov.b.status'] <- 'B'

assignments[assignments[, 'b.status.kmer'] != assignments[, 'old.cov.b.status'], ]

assignments[, c('A.lo_frac', 'B.lo_frac', 'A.st_frac', 'B.st_frac')] <- assignments[, c('A.lo', 'B.lo', 'A.st', 'B.st')] / assignments[, 'length']

naive_AB <- apply(assignments[, c('A.lo', 'B.lo')], 1, which.max)
score <- apply(assignments[, c('A.lo', 'B.lo')], 1, max) / rowSums(assignments[, c('A.lo', 'B.lo')])
B_score <- assignments[, 'B.lo'] / rowSums(assignments[, c('A.lo', 'B.lo')])

hist(score, breaks = 60)
hist(score[naive_AB == 2], breaks = 60, add = F)

hist(B_score[assignments[, 'old.cov.b.status'] != 'A'], breaks = 60)
hist(B_score[assignments[, 'new.cov.b.status'] != 'A'], breaks = 60, add = T, col = 'yellow')
hist(B_score[assignments[, 'b.status.kmer'] != 'A'], breaks = 60, add = T, col = 'purple')
hist(B_score[assignments[, 'b.status.asn'] != 'A'], breaks = 60, add = T, col = 'cyan')

plot(assignments[assignments[, 'old.cov.b.status'] != 'A', 'length'], B_score[assignments[, 'old.cov.b.status'] != 'A'], pch = 20, cex = 1.3)
points(assignments[assignments[, 'new.cov.b.status'] != 'A', 'length'], B_score[assignments[, 'new.cov.b.status'] != 'A'], pch = 20, cex = 1, col = 'yellow')
points(assignments[assignments[, 'b.status.kmer'] != 'A', 'length'], B_score[assignments[, 'b.status.kmer'] != 'A'], pch = 20, cex = 0.7, col = 'purple')
points(assignments[assignments[, 'b.status.asn'] != 'A', 'length'], B_score[assignments[, 'b.status.asn'] != 'A'], pch = 20, cex = 0.4, col = 'cyan')


plot(assignments[assignments[, 'old.cov.b.status'] != 'A', 'length'], assignments[assignments[, 'old.cov.b.status'] != 'A', 'fraction_of_B_windows'], pch = 20, cex = 1.3)
points(assignments[assignments[, 'new.cov.b.status'] != 'A', 'length'], assignments[assignments[, 'new.cov.b.status'] != 'A', 'fraction_of_B_windows'], pch = 20, cex = 1, col = 'yellow')
points(assignments[assignments[, 'b.status.kmer'] != 'A', 'length'], assignments[assignments[, 'b.status.kmer'] != 'A', 'fraction_of_B_windows'], pch = 20, cex = 0.7, col = 'purple')
points(assignments[assignments[, 'b.status.asn'] != 'A', 'length'], assignments[assignments[, 'b.status.asn'] != 'A', 'fraction_of_B_windows'], pch = 20, cex = 0.4, col = 'cyan')
