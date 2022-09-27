library('plotrix')

strains <- c('PV_18-04', 'PV_18-13', 'PV_18-21', 'PV_18-23')
coverage_table <- read.table('output/B_scaffold_assignment_comeplete_window_coverage_table.tsv', header = T, sep = '\t')
asn_table <- read.table('output/scaffolds.final.assignment.tsv', header = T, sep = '\t')

# general graph properties
general_cex = 1.2
B_candidate_col = 'darkgoldenrod1'
B_conf_candidate_col = "cadetblue"
B_col = "royalblue4"
B_c_col = "deepskyblue"

plot_pdf = T


###########
# Panel A #
###########
# for publication
Bc_windows <- coverage_table$scf %in% asn_table[asn_table$b.status.final == 'Bc', 'seq']
B_windows <- coverage_table$scf %in% asn_table[asn_table$b.status.final == 'B', 'seq']

if ( plot_pdf ){
  pdf('misc/asn_PV04_PV13_nomalised_coverage_windows_pub.pdf')
} else {
  tiff('misc/asn_PV04_PV13_nomalised_coverages_windows_pub.tiff', width = 8, height = 8, units = 'in', res = 150)
}

    plot(coverage_table$norm_cov_13, coverage_table$norm_cov_04, xlim = c(0, 300), ylim = c(0, 300), pch = 20, cex = 0.6, xlab = 'PV_18-13 nomalized coverage', ylab = 'PV_18-04 nomalized coverage', cex.axis = general_cex, cex.lab = general_cex, bty = 'n', col = 'grey80')

    points(coverage_table$norm_cov_13[Bc_windows], coverage_table$norm_cov_04[Bc_windows], col = B_c_col, pch = 20, cex = 0.6)
    points(coverage_table$norm_cov_13[B_windows], coverage_table$norm_cov_04[B_windows], col = B_col, pch = 20, cex = 0.6)

    lines(c(0, 1000), c(0, 1000))
    # lines(c(0, 1000), c(0, 2000), lty = 2)
    lines(c(0, 1000), c(0, 3000), lty = 3)
    # legend('topright', bty = 'n', c('1:1 (repetitive autosomal)', '1:3 (putatively B-linked)'), title = 'coverage ratio', lty = c(1, 3), cex = 1.3)
    legend(105, 315, bty = 'n', c('1:1 (repetitive autosomal)', '1:3 (putatively B-linked)'), title = 'coverage ratio', lty = c(1, 3), cex = general_cex)
dev.off()

Bc_scaffolds <- asn_table$b.status.final == 'Bc'
B_scaffolds <- asn_table$b.status.final == 'B'

if ( plot_pdf ){
  pdf('manuscript/figures_revision/asn_PV04_PV13_nomalised_coverages_pub.pdf')
} else {
  tiff('manuscript/figures_revision/asn_PV04_PV13_nomalised_coverages_pub.tiff', width = 8, height = 8, units = 'in', res = 150)
}
    plot(asn_table$norm_cov_PV13, asn_table$norm_cov_PV04, xlim = c(0, 100), ylim = c(0, 100), pch = 20, cex = 2, xlab = 'PV_18-13 nomalized coverage', ylab = 'PV_18-04 nomalized coverage', cex.axis = general_cex, cex.lab = general_cex, bty = 'n', col = 'grey80')


    points(asn_table$norm_cov_PV13[Bc_scaffolds], asn_table$norm_cov_PV04[Bc_scaffolds], col = B_c_col, pch = 20, cex = 2)
    points(asn_table$norm_cov_PV13[B_scaffolds], asn_table$norm_cov_PV04[B_scaffolds], col = B_col, pch = 20, cex = 2)

    lines(c(0, 1000), c(0, 1000))
    # lines(c(0, 1000), c(0, 2000), lty = 2)
    lines(c(0, 1000), c(0, 3000), lty = 3)
    # legend('topright', bty = 'n', c('1:1 (repetitive autosomal)', '1:3 (putatively B-linked)'), title = 'coverage ratio', lty = c(1, 3), cex = 1.3)
    legend(105, 315, bty = 'n', c('1:1 (repetitive autosomal)', '1:3 (putatively B-linked)'), title = 'coverage ratio', lty = c(1, 3), cex = general_cex)
dev.off()

# low resultion plot
png('manuscript/figures_revision/asn_PV04_PV13_nomalised_coverages.png')
    plot(coverage_table$norm_cov_13, coverage_table$norm_cov_04, xlim = c(0, 300), ylim = c(0, 300), pch = 20, cex = 0.6, xlab = 'PV_18-13 nomalized coverage', ylab = 'PV_18-04 nomalized coverage', cex.axis = 1.3, cex.lab = 1.3)
    lines(c(0, 1000), c(0, 1000))
    # lines(c(0, 1000), c(0, 2000), lty = 2)
    lines(c(0, 1000), c(0, 3000), lty = 3)
    legend('topright', bty = 'n', c('1:1 (repetitive autosomal)', '1:3 (putatively B-linked)'), title = 'coverage ratio', lty = c(1, 3), cex = 1.3)
dev.off()


png('manuscript/figures_revision/asn_PV04_PV13_nomalised_coverages_zoomin.png')
  subset_to_plot <- sample(1:nrow(coverage_table), 5000)
  plot(coverage_table$norm_cov_13[subset_to_plot], coverage_table$norm_cov_04[subset_to_plot], xlim = c(0, 10), ylim = c(0, 10), cex = 0.6, xlab = 'PV_18-13 nomalized coverage', ylab = 'PV_18-04 nomalized coverage', cex.axis = 1.3, cex.lab = 1.3)
  lines(c(0, 1000), c(0, 1000))
  lines(c(0, 1000), c(0, 3000), lty = 3)
  legend('topright', bty = 'n', c('1:1 (repetitive autosomal)', '1:3 (putatively B-linked)'), title = 'coverage ratio', lty = c(1, 3), cex = 1.3)
dev.off()

# plot(coverage_table$norm_B_cov_13, coverage_table$norm_B_cov_04, xlim = c(-5, 200), ylim = c(-5, 200))

# high_copy_number_B <- coverage_table[(coverage_table$norm_cov_13 / (coverage_table$norm_cov_04 + 1e-6)) < 0.65 & coverage_table$norm_cov_04 > 25 & coverage_table$norm_B_cov_13 > 5, ]
#
# plot(high_copy_number_B$norm_cov_04 ~ high_copy_number_B$norm_cov_13)



# hist(coverage_table$PV04_PV13_B_cov_ratio[coverage_table$PV04_PV13_B_cov_ratio > -5 & coverage_table$PV04_PV13_B_cov_ratio < 5 & coverage_table$scf == 'scaffold_362'])
#
# hist(coverage_table$norm_B_cov_04[coverage_table$norm_B_cov_04 > -2 & coverage_table$norm_B_cov_04 < 4], breaks = 50)
#
# hist(coverage_table$norm_B_cov_13[coverage_table$norm_B_cov_13 > -2 & coverage_table$norm_B_cov_13 < 4], breaks = 50)
points_to_plot <- coverage_table$norm_B_cov_04 > 1 | coverage_table$norm_B_cov_13 > 0.5
potentially_B <- coverage_table$norm_B_cov_04 > 1 & coverage_table$norm_B_cov_13 > 0.5

pdf('manuscript/figures_revision/asn_window_coverages_Bline_specific_vs_Blinked.pdf')
  plot(coverage_table$norm_B_cov_13[points_to_plot], coverage_table$norm_B_cov_04[points_to_plot], xlim = c(-5, 200), ylim = c(-5, 200), xlab = 'PV_18-13 nomalized coverage', ylab = 'PV_18-04 nomalized coverage', cex.axis = 1.3, cex.lab = 1.3)
  points(coverage_table$norm_B_cov_13[potentially_B], coverage_table$norm_B_cov_04[potentially_B], col = 'red', pch = 20)
  points(coverage_table[coverage_table$asn == 'B', 'norm_B_cov_13'], coverage_table[coverage_table$asn == 'B', 'norm_B_cov_04'], col = 'yellow', pch = 20, cex = 0.7)
  legend('topright', bty = 'n', pch = c(1, 20, 20), col = c(1, 'red', 'yellow'), c('considered windows', 'B-lines associated windows', 'Windows of B-linked scaffolds'), cex = 1.3)
dev.off()

pdf('manuscript/figures_revision/asn_window_coverages_Bline_specific_vs_Blinked_zoomed.pdf')
  plot(coverage_table$norm_B_cov_13[points_to_plot], coverage_table$norm_B_cov_04[points_to_plot], xlim = c(-1, 8), ylim = c(-1, 8), xlab = 'PV_18-13 nomalized coverage', ylab = 'PV_18-04 nomalized coverage', cex.axis = 1.3, cex.lab = 1.3)
  points(coverage_table$norm_B_cov_13[potentially_B], coverage_table$norm_B_cov_04[potentially_B], col = 'red', pch = 20)
  points(coverage_table[coverage_table$asn == 'B', 'norm_B_cov_13'], coverage_table[coverage_table$asn == 'B', 'norm_B_cov_04'], col = 'yellow', pch = 20, cex = 0.7)
  legend('topright', bty = 'n', pch = c(1, 20, 20), col = c(1, 'red', 'yellow'), c('considered windows', 'B-lines associated windows', 'Windows of B-linked scaffolds'), cex = 1.3)
dev.off()

###########
# Panel B #
###########

strains_short <- paste0('PV', rev(c('04', '13', '21', '23')))
white_bar_sizes <- c(sapply(paste0(strains_short, '.asn'), function(x){table(asn_table[, x])}), nrow(asn_table))
names(white_bar_sizes) <- c(strains_short, 'Reference')

B_scaffolds <- sum(asn_table$b.status.asn == 'B')
shared_by <- 4 - rowSums(is.na(asn_table[, paste0(strains_short, '.asn')]))
shared_by_all <- sum(shared_by == 4)
unique_to_reference <- sum(shared_by == 0)
in_at_least_one <- sum(shared_by > 0)

# moreless the original one
# tiff('misc/asn_asm_B_status_scaffolds.tiff', width = 8, height = 8, units = 'in', res = 150)
#     par(mar=c(5.1,6.1,4.1,2.1))
#     bar_plot_info <- barplot(white_bar_sizes, horiz = T, col = 'white', xlab = 'Scaffolds', cex.axis = general_cex, cex.lab = general_cex, axisnames = F)
#     axis(2, at = bar_plot_info, labels = names(white_bar_sizes), las = 2, cex.axis = general_cex, lwd = 0)
#     barplot(c(rep(NA, 4), in_at_least_one), horiz = T, col = 'gray40', add = T, axes = F)
#     barplot(c(rep(NA, 4), shared_by_all + B_scaffolds), horiz = T, col = B_candidate_col, add = T, axes = F)
#     barplot(c(rep(NA, 4), shared_by_all), horiz = T, col = 'gray80', add = T, axes = F)
# dev.off()

if ( plot_pdf ){
  pdf('manuscript/figures_revision/asn_asm_B_status_scaffolds.pdf')
} else {
  tiff('manuscript/figures_revision/asn_asm_B_status_scaffolds.tiff', width = 8, height = 8, units = 'in', res = 150)
}

    par(mar=c(5.1,6.1,4.1,2.1))
    bar_plot_info <- barplot(white_bar_sizes, horiz = T, col = c(rep('white', 4), 'gray40'), xlab = 'Scaffolds', cex.axis = general_cex, cex.lab = general_cex, axisnames = F)
    axis(2, at = bar_plot_info, labels = names(white_bar_sizes), las = 2, cex.axis = general_cex, lwd = 0)
    barplot(c(rep(NA, 4), in_at_least_one), horiz = T, col = 'white', add = T, axes = F)
    barplot(c(rep(NA, 2), rep(shared_by_all + B_scaffolds, 3)), horiz = T, col = B_candidate_col, add = T, axes = F)
    barplot(c(rep(shared_by_all, 5)), horiz = T, col = 'gray80', add = T, axes = F)
dev.off()

###########
# Panel C #
###########

total_kmers <- log10(asn_table[, 'A.lo'] + asn_table[, 'B.lo'])
kmer_ratio  <- asn_table[, 'AB.ratio.lo']
B_candidates <- (kmer_ratio > 0)

if ( plot_pdf ){
  pdf('manuscript/figures_revision/asn_kmer_B_status_scaffolds.pdf')
} else {
  tiff('manuscript/figures_revision/asn_kmer_B_status_scaffolds.tiff', width = 8, height = 8, units = 'in', res = 150)
}
  plot(total_kmers ~ kmer_ratio, bty = 'n', xlab = 'log2 (candidate B/A k-mers)', ylab = 'Total number of mapped 27-mers (log10)',
       cex.axis = general_cex, cex.lab = general_cex, pch = 20, cex = 2, col = 'grey80')
  points(total_kmers[B_candidates] ~ kmer_ratio[B_candidates], pch = 20, cex = 2, col = B_candidate_col)
  lines(c(0, 0), c(-10, 10), lty = 2, lwd = 1.5, col = 'red')
dev.off()

###########
# Panel D #
###########
scaffolds <- unique(coverage_table$scf)
scf2fraction_of_candidates <- function(scf){
	mean(coverage_table[coverage_table$scf == scf, 'B_candidate'])
}
scaffold_B_candidate_fraction <- sapply(scaffolds, scf2fraction_of_candidates)

fraction_of_w_windows_histogram <- hist(scaffold_B_candidate_fraction, plot = F)
xlab <- 'Fraction of B-supported windows'
ylab <- 'Number of scaffolds'

if ( plot_pdf ){
  pdf('manuscript/figures_revision/asn_fraction_of_B_windows_in_a_scaffold.pdf')
} else {
  tiff('manuscript/figures_revision/asn_fraction_of_B_windows_in_a_scaffold.tiff', width = 8, height = 8, units = 'in', res = 150)
}
    # native gap.barplot
    # gap.barplot(fraction_of_w_windows_histogram$counts, gap = c(220, 2100), ylim = c(0, 400), xtics = seq(1,10, by = 1), xaxlab=NA, ytics = c(0, 100, 200, 300, 500, 2200), xlab = xlab, ylab = ylab, col = c('grey', rep('yellow', 5), rep('green', 4)))

    # home-implemented barplot
    g_from <- 220
    g_to <- 2170
    values_to_shift <- which(fraction_of_w_windows_histogram$counts > g_to)
    fraction_of_w_windows_histogram$counts[values_to_shift] <- fraction_of_w_windows_histogram$counts[values_to_shift] - (g_to - g_from)
    ymax <- 2300

    gap.plot(100, c(-500), gap=c(g_from, g_to), xlim = c(0, 1), ylim = c(0, ymax),
    xlab = xlab, ylab = ylab, cex.lab = general_cex)
    rect(par("usr")[1], par("usr")[3], par("usr")[2], g_from, col=rgb(0.92,0.92,0.92))
    rect(par("usr")[1], g_from*(1+0.02), par("usr")[2], par("usr")[4], col=rgb(0.92,0.92,0.92))
    plot(fraction_of_w_windows_histogram, add = T, col = c('grey', rep(B_candidate_col, 5), rep(B_conf_candidate_col, 4)))
    axis.break(2, g_from, breakcol="snow", style="gap")
    axis.break(2, g_from*(1+0.02), breakcol="black", style="slash")
    axis.break(4, g_from*(1+0.02), breakcol="black", style="slash")
    axis(2, seq(0, g_from, by = 100), col = NA, col.ticks = 1, cex.axis = general_cex)
    axis(2, seq(g_from + 30, g_from + ymax - g_to + 30, by = 100), col = NA, col.ticks = 1, labels = seq(g_to + 30, ymax, by = 100), cex.axis = general_cex)
    axis(1, seq(0, 1, by = 0.1), col = NA, col.ticks = 1, cex.axis = general_cex)

    # Log plot
    # fraction_of_w_windows_histogram$counts <- sapply(fraction_of_w_windows_histogram$counts, function(x){ifelse(x == 0, 0, log10(x))})
    # plot(fraction_of_w_windows_histogram)
dev.off()
