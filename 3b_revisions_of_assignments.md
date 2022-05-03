#### Extracting per-window coverage

Using `PV_18-??.initial.sorted.primary.only.bam` files (`/data/ross/mealybugs/analyses/B_viburni_2020/4_cov_analysis/cov` on cluster). I extracted-per 10k window coverages using `samtools depth` and I perl script. Like this (for all 4 strains)

```bash
qsub -o logs -e logs -cwd -b yes -N depth -V -pe smp64 1 'samtools depth PV_18-23.initial.sorted.primary.only.bam | perl scripts/depth2windows.pl 10000 > PV_18-23.initial.sorted.primary.only_window_coverages.tsv'
```

#### Creating coverage-assignments

I will load the tables in R and merge the four strains in a single table (`scf`, `window`, `cov_04`, ...)

```R
strains <- c('PV_18-04', 'PV_18-13', 'PV_18-21', 'PV_18-23')

coverage_files <- paste0('cov/', strains, '.initial.sorted.primary.only_window_coverages.tsv')

load_coverage_table <- function(filename){
  tab <- read.table(filename, header = F, col.names = c('scf', 'range', 'cov'))
  return(tab)
}

coverage_tabs <- lapply(coverage_files, load_coverage_table)
monoploid_cov_estimates <- sapply(coverage_tabs, function(x){ median(x[, 3]) } ) / 2
```

Now I have list with the coverage of all windows, all I need to do now is to merge them.

```{R}
row_names <- lapply(coverage_tabs, function(x){ paste(x[, 'scf'], x[, 'range']) } ) # creating a list of row names of the four table

all_windows <- unique(unlist(row_names)) # merging them into one unique set (the window is not present if no read maps there, which is the reason why i must merge them first)

scfs <- sapply(strsplit(all_windows, split = ' '), function(x){ return(x[1]) } ) # then we cut the scaffold name out of the IDs we just created
windows <- sapply(strsplit(all_windows, split = ' '), function(x){ return(x[2]) } ) # the same with windows

coverage_table <- data.frame(scf = scfs, window = windows) # making the table
rownames(coverage_table) <- all_windows # naming the rows so we can easily add there the infomation
coverage_table <- coverage_table[order(sapply(strsplit(scfs, '_'), function(x){as.numeric(x[2])} )), ] # just reordering the table so the rows are sorted by the scaffold ID

# the table is ready, now filling with individual data
coverage_table[, c('cov_04', 'cov_13', 'cov_21', 'cov_23')] <- 0 # defaut coverage is 0 (would be NA if not initiated by 0s)

coverage_table[row_names[[1]], 'cov_04'] <- coverage_tabs[[1]][, 'cov']
coverage_table[row_names[[2]], 'cov_13'] <- coverage_tabs[[2]][, 'cov']
coverage_table[row_names[[3]], 'cov_21'] <- coverage_tabs[[3]][, 'cov']
coverage_table[row_names[[4]], 'cov_23'] <- coverage_tabs[[4]][, 'cov']
rownames(coverage_table) <- 1:nrow(coverage_table) # the row names are not needed anymore

coverage_table[, c('norm_cov_04', 'norm_cov_13', 'norm_cov_21', 'norm_cov_23')] <- t(t(coverage_table[, c('cov_04', 'cov_13', 'cov_21', 'cov_23')]) / monoploid_cov_estimates)

hist(coverage_table$norm_cov_13, xlim = c(0,6), breaks = 2000)

plot(coverage_table$norm_cov_21, coverage_table$norm_cov_04)

plot(coverage_table$norm_cov_13, coverage_table$norm_cov_04)
lines(c(0, 1000), c(0, 1000))
lines(c(0, 1000), c(0, 2000), lty = 2)
lines(c(0, 1000), c(0, 3000), lty = 3)

plot(coverage_table$norm_cov_13, coverage_table$norm_cov_04, xlim = c(0, 10), ylim = c(0, 10))
lines(c(0, 1000), c(0, 1000))
lines(c(0, 1000), c(0, 2000), lty = 2)
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

head(coverage_table[coverage_table$scf == 'scaffold_360', ])

coverage_table$PV04_PV13_B_cov_ratio <- (coverage_table$norm_B_cov_13 / coverage_table$norm_B_cov_04)

coverage_table$scf == 'scaffold_360'

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

B_scaffolds <- scaffolds[scaffold_B_candidate_fraction > 0.6]

old_asn <- read.csv('../5_B_char_heterozygosity_in_B04/scaffolds.final.assignment.table.csv')

old_asn$new_asn <- 'A'
old_asn[old_asn$seq %in% B_scaffolds, 'new_asn'] <- 'B'


old_asn[old_asn$new_asn == 'B', ]

old_asn[old_asn$b.status.final %in% c('B1', 'B2', 'B3') & old_asn$new_asn == 'A', ]

old_asn[old_asn$new_asn == 'B', 'length']
```
