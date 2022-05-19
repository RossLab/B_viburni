strains <- c('PV_18-04', 'PV_18-13', 'PV_18-21', 'PV_18-23')
coverage_files <- paste0('cov/', strains, '.freeze.v0.w1000_cov.bed')

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

hist(coverage_table$norm_cov_13, xlim = c(0,6), breaks = 2000)
plot(coverage_table$norm_cov_21, coverage_table$norm_cov_04)

plot(coverage_table$norm_cov_13, coverage_table$norm_cov_04, xlim = c(0, 300), ylim = c(0, 300))
lines(c(0, 1000), c(0, 1000))
# lines(c(0, 1000), c(0, 2000), lty = 2)
lines(c(0, 1000), c(0, 3000), lty = 3)

plot(coverage_table$norm_cov_13, coverage_table$norm_cov_04, xlim = c(0, 10), ylim = c(0, 10))
lines(c(0, 1000), c(0, 1000))
lines(c(0, 1000), c(0, 3000), lty = 3)

# coverage_table$consistent_no_B_ploidy <- Inf
# for (i in c(0, 2, 4, 6){
# 	coverage_table$consistent_no_B_ploidy[round(coverage_table$norm_cov_21) == i & round(coverage_table$norm_cov_23) == i] <- i
# }

# coverage_table[round(coverage_table$norm_cov_21) != round(coverage_table$norm_cov_23), ]

# informative_windows <- coverage_table[coverage_table$consistent_no_B_ploidy != Inf, ]

mean_autosomal_coverages <- rowMeans(coverage_table[, c('norm_cov_21', 'norm_cov_23')])
coverage_table$norm_B_cov_04 <- coverage_table$norm_cov_04 - mean_autosomal_coverages
coverage_table$norm_B_cov_13 <- coverage_table$norm_cov_13 - mean_autosomal_coverages

plot(coverage_table$norm_B_cov_04, coverage_table$norm_B_cov_13, xlim = c(-5, 200), ylim = c(-5, 100))

plot(coverage_table$norm_cov_04, coverage_table$norm_cov_13, xlim = c(-5, 200), ylim = c(-5, 100))

high_copy_number_B <- coverage_table[(coverage_table$norm_cov_13 / (coverage_table$norm_cov_04 + 1e-6)) < 0.65 & coverage_table$norm_cov_04 > 25 & coverage_table$norm_B_cov_13 > 5, ]

plot(high_copy_number_B$norm_cov_04 ~ high_copy_number_B$norm_cov_13)


coverage_table$PV04_PV13_B_cov_ratio <- (coverage_table$norm_B_cov_13 / coverage_table$norm_B_cov_04)
hist(coverage_table$PV04_PV13_B_cov_ratio[coverage_table$PV04_PV13_B_cov_ratio > -5 & coverage_table$PV04_PV13_B_cov_ratio < 5 & coverage_table$scf == 'scaffold_362'])

###############################

hist(coverage_table$norm_B_cov_04[coverage_table$norm_B_cov_04 > -2 & coverage_table$norm_B_cov_04 < 4], breaks = 50)
hist(coverage_table$norm_B_cov_13[coverage_table$norm_B_cov_13 > -2 & coverage_table$norm_B_cov_13 < 4], breaks = 50)

potentially_B <- coverage_table$norm_B_cov_04 > 1 & coverage_table$norm_B_cov_13 > 0.5

hist(coverage_table$PV04_PV13_B_cov_ratio[potentially_B], breaks = 500, xlim = c(0, 2))
# THIS MAKES A LOT OF SENSE!!!

coverage_table$B_candidate <- coverage_table$PV04_PV13_B_cov_ratio < 0.7 & potentially_B

scaffolds <- unique(coverage_table$scf)

scf <- 'scaffold_360'
scf2fraction_of_candidates <- function(scf){
	mean(coverage_table[coverage_table$scf == scf, 'B_candidate'])
}

scaffold_B_candidate_fraction <- sapply(scaffolds, scf2fraction_of_candidates)
hist(scaffold_B_candidate_fraction, ylim = c(0, 100))

B_scaffolds <- scaffolds[scaffold_B_candidate_fraction > 0.6]

old_asn <- read.csv('../5_B_char_heterozygosity_in_B04/scaffolds.final.assignment.table.csv')

old_asn$new_asn <- 'A'
old_asn[old_asn$seq %in% B_scaffolds, 'new_asn'] <- 'B'


old_asn[old_asn$new_asn == 'B', ]

old_asn[old_asn$b.status.final %in% c('B1', 'B2', 'B3') & old_asn$new_asn == 'A', ]

old_asn[old_asn$new_asn == 'B', 'length']

DE_genes <- read.csv('../../output/B_diff_expr/FB.vs.FnoB.de.annotated.csv')

DE_genes$new_asn <- 'A'
DE_genes[DE_genes$seq.y %in% B_scaffolds, 'new_asn'] <- 'B'
