strains <- c('PV_18-04', 'PV_18-13', 'PV_18-21', 'PV_18-23')

coverage_files <- paste0('data/4_cov_analysis/cov/', strains, '.freeze.v0.w10000_cov.bed')

load_coverage_table <- function(filename){
  tab <- read.table(filename, header = F, col.names = c('scf', 'from', 'to', 'cov'))
  return(tab)
}

coverage_tabs <- lapply(coverage_files, load_coverage_table)

monoploid_cov_estimates <- sapply(coverage_tabs, function(x){ median(x[, 4]) } ) / 2

coverage_table <- data.frame(scf = coverage_tabs[[1]][, 'scf'],
                             from = coverage_tabs[[1]][, 'from'],
                             to = coverage_tabs[[1]][, 'to'])

# the table is ready, now filling with individual data
coverage_table[, c('cov_04', 'cov_13', 'cov_21', 'cov_23')] <- 0 # defaut coverage is 0 (would be NA if not initiated by 0s)
coverage_table[, 'cov_04'] <- coverage_tabs[[1]][, 'cov']
coverage_table[, 'cov_13'] <- coverage_tabs[[2]][, 'cov']
coverage_table[, 'cov_21'] <- coverage_tabs[[3]][, 'cov']
coverage_table[, 'cov_23'] <- coverage_tabs[[4]][, 'cov']

coverage_table[, c('norm_cov_04', 'norm_cov_13', 'norm_cov_21', 'norm_cov_23')] <- t(t(coverage_table[, c('cov_04', 'cov_13', 'cov_21', 'cov_23')]) / monoploid_cov_estimates)

mean_autosomal_coverages <- rowMeans(coverage_table[, c('norm_cov_21', 'norm_cov_23')])

coverage_table$norm_B_cov_04 <- coverage_table$norm_cov_04 - mean_autosomal_coverages
coverage_table$norm_B_cov_13 <- coverage_table$norm_cov_13 - mean_autosomal_coverages

potentially_B <- coverage_table$norm_B_cov_04 > 1 & coverage_table$norm_B_cov_13 > 0.5

coverage_table$PV04_PV13_B_cov_ratio <- (coverage_table$norm_B_cov_13 / coverage_table$norm_B_cov_04)
coverage_table$B_candidate <- coverage_table$PV04_PV13_B_cov_ratio < 0.7 & potentially_B

scaffolds <- unique(coverage_table$scf)
scf2fraction_of_candidates <- function(scf){
	mean(coverage_table[coverage_table$scf == scf, 'B_candidate'])
}

scaffold_B_candidate_fraction <- sapply(scaffolds, scf2fraction_of_candidates)
# hist(scaffold_B_candidate_fraction)
B_scaffolds <- scaffolds[scaffold_B_candidate_fraction > 0.6]

coverage_table[, 'asn'] <- 'A'
coverage_table[coverage_table$scf %in% B_scaffolds, 'asn'] <- 'B'

write.table(coverage_table, 'output/B_scaffold_assignment_comeplete_window_coverage_table.tsv', sep = '\t', quote = F, row.names = F, col.names = T)
