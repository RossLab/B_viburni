assignments <- read.csv('output/scaffolds.preprint.assignment.csv')
row.names(assignments) <- assignments$seq
head(assignments)

# removing all the old info about mapping that is obsolete now
assignments <- assignments[, c(1, 2, 14:18, 20:27)]
# head(assignments)
colnames(assignments)[15] <- 'b.status.old.final'

# generated together with plots in R_scripts/B_scaffolds_assignment_plots.R
coverage_table <- read.table('output/B_scaffold_assignment_comeplete_window_coverage_table.tsv', header = T)
scaffolds <- assignments$seq
scf2mean_property <- function(scf, prop = 'B_candidate'){
	mean(coverage_table[coverage_table$scf == scf, prop])
}

assignments[scaffolds, 'norm_cov_PV04'] <- sapply(scaffolds, scf2mean_property, 'norm_cov_04')
assignments[scaffolds, 'norm_cov_PV13'] <- sapply(scaffolds, scf2mean_property, 'norm_cov_13')
assignments[scaffolds, 'norm_cov_PV21'] <- sapply(scaffolds, scf2mean_property, 'norm_cov_21')
assignments[scaffolds, 'norm_cov_PV23'] <- sapply(scaffolds, scf2mean_property, 'norm_cov_23')

scaffold_B_candidate_fraction <- sapply(scaffolds, scf2mean_property)
# hist(scaffold_B_candidate_fraction)
potentially_B_scaffolds <- scaffolds[scaffold_B_candidate_fraction > 0.1]
B_scaffolds <- scaffolds[scaffold_B_candidate_fraction > 0.6]

assignments[scaffolds, 'fraction_of_B_windows'] <- round(scaffold_B_candidate_fraction, 4)
assignments[, 'cov.b.status'] <- 'A'
assignments[potentially_B_scaffolds, 'cov.b.status'] <- 'Bc'
assignments[B_scaffolds, 'cov.b.status'] <- 'B'

assignments$b.status.final <- ifelse(assignments$cov.b.status == "B" &
                                     (assignments$b.status.asn == "B" | assignments$b.status.kmer == "B"),
                                     "B", "A")
assignments$b.status.final <- ifelse(assignments$b.status.final == "A" &
                                     (assignments$b.status.asn == "B" |
                                      assignments$cov.b.status != "A" |
                                      assignments$b.status.kmer == "B"),
                                     "Bc", assignments$b.status.final)

assignments$assignment_changes <- paste(assignments$b.status.old.final, assignments$b.status.final, sep = ' -> ')
assignment_changes_table <- table(assignments$assignment_changes)

assignment_change_size <- sapply(names(assignment_changes_table), function(x){ round(sum(assignments[assignments$assignment_changes == x, 'length']) / 1e6, 2) } )
data.frame(change = names(assignment_changes_table), number_of_scaffolds = as.vector(assignment_changes_table), length = as.vector(assignment_change_size))

# assignments <- read.table('output/scaffolds.final.assignment.tsv', sep = '\t', header = T)
sapply(c('A', 'B1', 'B2', 'B3'), function(x){ round(sum(assignments[assignments$b.status.old.final == x, 'length']) / 1e6, 2) } )
sapply(c('A', 'B', 'Bc'), function(x){ round(sum(assignments[assignments$b.status.final == x, 'length']) / 1e6, 2) } )

write.table(assignments, file = 'output/scaffolds.final.assignment.tsv', quote = F, row.names = F, sep = '\t', col.names = T)
write.table(assignments[, c('seq', 'length', 'b.status.final')], 'output/scaffolds.final.assignment.table.csv', quote = F, row.names = F, sep = ';', col.names = T)
