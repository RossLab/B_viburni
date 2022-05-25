#### Extracting per-window coverage

Using `PV_18-??.initial.sorted.primary.only.bam` files (`/data/ross/mealybugs/analyses/B_viburni_2020/4_cov_analysis/cov` on cluster). I extracted-per 10k window coverages using `samtools depth` and I perl script. Like this (for all 4 strains)

```bash
for strain in 04 13 21 23; do
 qsub -o logs -e logs -cwd -b yes -N depth -V -pe smp64 1 "samtools depth -aa PV_18-$strain.freeze.v0.sorted.bam | scripts/depth2bed_coverage.py -b PV_18-$strain.freeze.v0.sorted.bam -w 10000 > PV_18-$strain.freeze.v0.w10000_cov.bed"
done

for strain in 04 13 21 23; do
 qsub -o logs -e logs -cwd -b yes -N depth -V -pe smp64 1 "samtools depth -aa PV_18-$strain.freeze.v0.sorted.bam | scripts/depth2bed_coverage.py -b PV_18-$strain.freeze.v0.sorted.bam -w 1000 > /scratch/$USER/PV_18-$strain.freeze.v0.w1000_cov.bed && rsync -av --remove-source-files /scratch/$USER/PV_18-$strain.freeze.v0.w1000_cov.bed ."
done
```

I will proceed with `w10000` files. using 1000bp long windows (`w1000`) is just a bit messier but resulting in moreless the same assignments (Exploration here: `R_scripts/w1000_coverage_files.R`)

#### Creating coverage-assignments

I will load the tables in R and merge the four strains in a single table (`scf`, `window`, `cov_04`, ...)

```R
strains <- c('PV_18-04', 'PV_18-13', 'PV_18-21', 'PV_18-23')

coverage_files <- paste0('cov/', strains, '.freeze.v0.w10000_cov.bed')

load_coverage_table <- function(filename){
  tab <- read.table(filename, header = F, col.names = c('scf', 'from', 'to', 'cov'))
  return(tab)
}

coverage_tabs <- lapply(coverage_files, load_coverage_table)
monoploid_cov_estimates <- sapply(coverage_tabs, function(x){ median(x[, 4]) } ) / 2
```

Now I have list with the coverage of all windows, all I need to do now is to merge them. They are constructed with the same reference, so there is a guarantee the order and number of windows is the same in all four files.

```{R}
coverage_table <- data.frame(scf = coverage_tabs[[1]][, 'scf'],
                             from = coverage_tabs[[1]][, 'from'],
                             to = coverage_tabs[[1]][, 'to'])

# the table is ready, now filling with individual data
coverage_table[, c('cov_04', 'cov_13', 'cov_21', 'cov_23')] <- 0 # defaut coverage is 0 (would be NA if not initiated by 0s)
coverage_table[, 'cov_04'] <- coverage_tabs[[1]][, 'cov']
coverage_table[, 'cov_13'] <- coverage_tabs[[2]][, 'cov']
coverage_table[, 'cov_21'] <- coverage_tabs[[3]][, 'cov']
coverage_table[, 'cov_23'] <- coverage_tabs[[4]][, 'cov']

# normalized coverage (corresponding moreless to copy-number)
coverage_table[, c('norm_cov_04', 'norm_cov_13', 'norm_cov_21', 'norm_cov_23')] <- t(t(coverage_table[, c('cov_04', 'cov_13', 'cov_21', 'cov_23')]) / monoploid_cov_estimates)

# hist(coverage_table$norm_cov_13, xlim = c(0,6), breaks = 2000)
#
# plot(coverage_table$norm_cov_21, coverage_table$norm_cov_04)

pdf('PV04_PV13_nomalised_coverages.pdf')
  plot(coverage_table$norm_cov_13, coverage_table$norm_cov_04, xlim = c(0, 300), ylim = c(0, 300), pch = 20, cex = 0.6, xlab = 'PV_18-13 nomalized coverage', ylab = 'PV_18-04 nomalized coverage', cex.axis = 1.3, cex.lab = 1.3)
  lines(c(0, 1000), c(0, 1000))
  # lines(c(0, 1000), c(0, 2000), lty = 2)
  lines(c(0, 1000), c(0, 3000), lty = 3)
  legend('topright', bty = 'n', c('1:1 (repetitive autosomal)', '1:3 (putatively B-linked)'), title = 'coverage ratio', lty = c(1, 3), cex = 1.3)
dev.off()

pdf('PV04_PV13_nomalised_coverages_zoomin.pdf')
  subset_to_plot <- sample(1:nrow(coverage_table), 5000)
  plot(coverage_table$norm_cov_13[subset_to_plot], coverage_table$norm_cov_04[subset_to_plot], xlim = c(0, 6), ylim = c(0, 6), cex = 0.6, xlab = 'PV_18-13 nomalized coverage', ylab = 'PV_18-04 nomalized coverage', cex.axis = 1.3, cex.lab = 1.3)
  lines(c(0, 1000), c(0, 1000))
  lines(c(0, 1000), c(0, 3000), lty = 3)
  legend('topright', bty = 'n', c('1:1 (repetitive autosomal)', '1:3 (putatively B-linked)'), title = 'coverage ratio', lty = c(1, 3), cex = 1.3)
dev.off()

mean_autosomal_coverages <- rowMeans(coverage_table[, c('norm_cov_21', 'norm_cov_23')])

coverage_table$norm_B_cov_04 <- coverage_table$norm_cov_04 - mean_autosomal_coverages
coverage_table$norm_B_cov_13 <- coverage_table$norm_cov_13 - mean_autosomal_coverages

plot(coverage_table$norm_B_cov_13, coverage_table$norm_B_cov_04, xlim = c(-5, 200), ylim = c(-5, 200))

# high_copy_number_B <- coverage_table[(coverage_table$norm_cov_13 / (coverage_table$norm_cov_04 + 1e-6)) < 0.65 & coverage_table$norm_cov_04 > 25 & coverage_table$norm_B_cov_13 > 5, ]
#
# plot(high_copy_number_B$norm_cov_04 ~ high_copy_number_B$norm_cov_13)

coverage_table$PV04_PV13_B_cov_ratio <- (coverage_table$norm_B_cov_13 / coverage_table$norm_B_cov_04)

hist(coverage_table$PV04_PV13_B_cov_ratio[coverage_table$PV04_PV13_B_cov_ratio > -5 & coverage_table$PV04_PV13_B_cov_ratio < 5 & coverage_table$scf == 'scaffold_362'])

###############################

hist(coverage_table$norm_B_cov_04[coverage_table$norm_B_cov_04 > -2 & coverage_table$norm_B_cov_04 < 4], breaks = 50)

hist(coverage_table$norm_B_cov_13[coverage_table$norm_B_cov_13 > -2 & coverage_table$norm_B_cov_13 < 4], breaks = 50)

points_to_plot <- coverage_table$norm_B_cov_04 > 1 | coverage_table$norm_B_cov_13 > 0.5
potentially_B <- coverage_table$norm_B_cov_04 > 1 & coverage_table$norm_B_cov_13 > 0.5

# hist(coverage_table$PV04_PV13_B_cov_ratio[potentially_B], breaks = 500, xlim = c(0, 2))
# THIS MAKES A LOT OF SENSE!!!

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

pdf('window_coverages_Bline_specific_vs_Blinked.pdf')
  plot(coverage_table$norm_B_cov_13[points_to_plot], coverage_table$norm_B_cov_04[points_to_plot], xlim = c(-5, 200), ylim = c(-5, 200), xlab = 'PV_18-13 nomalized coverage', ylab = 'PV_18-04 nomalized coverage', cex.axis = 1.3, cex.lab = 1.3)
  points(coverage_table$norm_B_cov_13[potentially_B], coverage_table$norm_B_cov_04[potentially_B], col = 'red', pch = 20)
  points(coverage_table[coverage_table$asn == 'B', 'norm_B_cov_13'], coverage_table[coverage_table$asn == 'B', 'norm_B_cov_04'], col = 'yellow', pch = 20, cex = 0.7)
  legend('topright', bty = 'n', pch = c(1, 20, 20), col = c(1, 'red', 'yellow'), c('considered windows', 'B-lines associated windows', 'Windows of B-linked scaffolds'), cex = 1.3)
dev.off()

pdf('window_coverages_Bline_specific_vs_Blinked_zoomed.pdf')
  plot(coverage_table$norm_B_cov_13[points_to_plot], coverage_table$norm_B_cov_04[points_to_plot], xlim = c(-1, 8), ylim = c(-1, 8), xlab = 'PV_18-13 nomalized coverage', ylab = 'PV_18-04 nomalized coverage', cex.axis = 1.3, cex.lab = 1.3)
  points(coverage_table$norm_B_cov_13[potentially_B], coverage_table$norm_B_cov_04[potentially_B], col = 'red', pch = 20)
  points(coverage_table[coverage_table$asn == 'B', 'norm_B_cov_13'], coverage_table[coverage_table$asn == 'B', 'norm_B_cov_04'], col = 'yellow', pch = 20, cex = 0.7)
  legend('topright', bty = 'n', pch = c(1, 20, 20), col = c(1, 'red', 'yellow'), c('considered windows', 'B-lines associated windows', 'Windows of B-linked scaffolds'), cex = 1.3)
dev.off()

write.table(coverage_table, 'Comeplete_window_coverage_table.tsv', sep = '\t', quote = F, row.names = F, col.names = T)
#
# old_asn <- read.csv('../5_B_char_heterozygosity_in_B04/scaffolds.final.assignment.table.csv')
#
# old_asn$new_asn <- 'A'
# old_asn[old_asn$seq %in% B_scaffolds, 'new_asn'] <- 'B'
#
#
# old_asn[old_asn$new_asn == 'B', ]
#
# old_asn[old_asn$b.status.final %in% c('B1', 'B2', 'B3') & old_asn$new_asn == 'A', ]
#
# old_asn[old_asn$new_asn == 'B', 'length']
#
# DE_genes <- read.csv('../../output/B_diff_expr/FB.vs.FnoB.de.annotated.csv')
#
# DE_genes$new_asn <- 'A'
# DE_genes[DE_genes$seq.y %in% B_scaffolds, 'new_asn'] <- 'B'

```

### K-mers

I followed the coverage analysis workflow and it seems there is the same problem as with mapping coverages. The B kmers are simply B-line enriched, but not necesarily B specific.

I will redo this analysis too with subtracting:
 - B-specific (absent in B- line, Enrichment in PV04 line compared to PV13)
 - B-nonspecific (B-enriched, Enrichment in PV04 line compared to PV13)
 - B-line-associated (B-enriched but with not cov. ratio skew)
 - Reliably Autosomal (similar cov levels in all samples)

But first I will take a sample of kmers and explore if it actually make sense.

```R
kmer_tab <- read.table('kmer/decon_kmers/kmer_dump_sample.dump')
colnames(kmer_tab) <- c('kmer', 'PV04', 'PV13', 'PV21', 'PV23')

# the kmer coverages are 1n estimates from GenomeScope
# however that was for k = 21, for k = 27, we need to correct by factor (R - 27 + 1) / (R - 21 + 1), where R is readlength and = 150
kmer_coverages <- c('PV04' = 42.5, "PV13" = 40.1, "PV21" = 22.7, "PV23" = 24.7) * (150 - 27 + 1) / (150 - 21 + 1)

# replacing all the NAs with 0s
kmer_tab[is.na(kmer_tab)] <- 0
# normalising kmer coverages
kmer_tab[, 2:5] <- t(t(kmer_tab[, 2:5]) / kmer_coverages)

kmer_tab[, c('PV04_B', 'PV13_B')] <- kmer_tab[, c('PV04', 'PV13')] - rowMeans(kmer_tab[, c('PV21', 'PV23')])

kmer_tab[, 2:5] <- round(kmer_tab[, 2:5])


# removing kmers that are likely errors (those were all of them have normalised rounded coverage == 0)
totally_absent <- rowSums(kmer_tab[, 2:5]) == 0
kmer_tab <- kmer_tab[!totally_absent, ]

# now, 20% of kmers are B specific
B_line_exclusive <- rowSums(kmer_tab[, c('PV21', 'PV23')]) == 0
B_minus_exclusive <- rowSums(kmer_tab[, c('PV04', 'PV13')]) == 0
B_enriched <- kmer_tab[, 'PV04_B'] > 0.5 & kmer_tab[, 'PV13_B'] > 0.5

hist(log2(kmer_tab[B_enriched, 'PV13_B'] / kmer_tab[B_enriched, 'PV04_B']), breaks = 100)

plot(kmer_tab[B_enriched, 'PV13_B'], kmer_tab[B_enriched, 'PV04_B'], xlim = c(0, 10), ylim = c(0, 10))
```

My conclusion here is that so much data transformation (subtraction, normalisation, coverage ratio) will cause too much coverage variation to measure clear peak of kmers that are ~3x more frequent in PV04 compared to PV13. So perhaps the way out will be use Andres' mapping and explore how the "enriched kmers" distribute on the assembled sequences.

### Re-evaluating fractions of assigned kmers

So instead I will at least check for considency of individual assignments with the suboptimal kmer approach we already have implemented (`R_scripts/explore_assignments.R`).

All those plots seem quite alright. So, I will keep those.
