strains <- c('PV_18-04', 'PV_18-13', 'PV_18-21', 'PV_18-23')
coverage_table <- read.table('output/B_scaffold_assignment_comeplete_window_coverage_table.tsv', header = T, sep = '\t')

make_pub <- T

if (make_pub){
  tiff('misc/asn_PV04_PV13_nomalised_coverages_pub.tiff', width = 8, height = 8, units = 'in', res = 150)
  plot(coverage_table$norm_cov_13, coverage_table$norm_cov_04, xlim = c(0, 300), ylim = c(0, 300), pch = 20, cex = 0.6, xlab = 'PV_18-13 nomalized coverage', ylab = 'PV_18-04 nomalized coverage', cex.axis = 1.2, cex.lab = 1.2)
  lines(c(0, 1000), c(0, 1000))
  # lines(c(0, 1000), c(0, 2000), lty = 2)
  lines(c(0, 1000), c(0, 3000), lty = 3)
  # legend('topright', bty = 'n', c('1:1 (repetitive autosomal)', '1:3 (putatively B-linked)'), title = 'coverage ratio', lty = c(1, 3), cex = 1.3)
  legend(105, 315, bty = 'n', c('1:1 (repetitive autosomal)', '1:3 (putatively B-linked)'), title = 'coverage ratio', lty = c(1, 3), cex = 1.2)
} else {
  png('misc/asn_PV04_PV13_nomalised_coverages.png')
    plot(coverage_table$norm_cov_13, coverage_table$norm_cov_04, xlim = c(0, 300), ylim = c(0, 300), pch = 20, cex = 0.6, xlab = 'PV_18-13 nomalized coverage', ylab = 'PV_18-04 nomalized coverage', cex.axis = 1.3, cex.lab = 1.3)
    lines(c(0, 1000), c(0, 1000))
    # lines(c(0, 1000), c(0, 2000), lty = 2)
    lines(c(0, 1000), c(0, 3000), lty = 3)
    legend('topright', bty = 'n', c('1:1 (repetitive autosomal)', '1:3 (putatively B-linked)'), title = 'coverage ratio', lty = c(1, 3), cex = 1.3)
}
dev.off()


png('misc/asn_PV04_PV13_nomalised_coverages_zoomin.png')
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

pdf('misc/asn_window_coverages_Bline_specific_vs_Blinked.pdf')
  plot(coverage_table$norm_B_cov_13[points_to_plot], coverage_table$norm_B_cov_04[points_to_plot], xlim = c(-5, 200), ylim = c(-5, 200), xlab = 'PV_18-13 nomalized coverage', ylab = 'PV_18-04 nomalized coverage', cex.axis = 1.3, cex.lab = 1.3)
  points(coverage_table$norm_B_cov_13[potentially_B], coverage_table$norm_B_cov_04[potentially_B], col = 'red', pch = 20)
  points(coverage_table[coverage_table$asn == 'B', 'norm_B_cov_13'], coverage_table[coverage_table$asn == 'B', 'norm_B_cov_04'], col = 'yellow', pch = 20, cex = 0.7)
  legend('topright', bty = 'n', pch = c(1, 20, 20), col = c(1, 'red', 'yellow'), c('considered windows', 'B-lines associated windows', 'Windows of B-linked scaffolds'), cex = 1.3)
dev.off()

pdf('misc/asn_window_coverages_Bline_specific_vs_Blinked_zoomed.pdf')
  plot(coverage_table$norm_B_cov_13[points_to_plot], coverage_table$norm_B_cov_04[points_to_plot], xlim = c(-1, 8), ylim = c(-1, 8), xlab = 'PV_18-13 nomalized coverage', ylab = 'PV_18-04 nomalized coverage', cex.axis = 1.3, cex.lab = 1.3)
  points(coverage_table$norm_B_cov_13[potentially_B], coverage_table$norm_B_cov_04[potentially_B], col = 'red', pch = 20)
  points(coverage_table[coverage_table$asn == 'B', 'norm_B_cov_13'], coverage_table[coverage_table$asn == 'B', 'norm_B_cov_04'], col = 'yellow', pch = 20, cex = 0.7)
  legend('topright', bty = 'n', pch = c(1, 20, 20), col = c(1, 'red', 'yellow'), c('considered windows', 'B-lines associated windows', 'Windows of B-linked scaffolds'), cex = 1.3)
dev.off()
