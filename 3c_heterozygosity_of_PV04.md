TODO: update this analysis or remove completelly

### The B heterozygosity of PV_04

If the two Bs are homologous, and if the higher coverage is due to reads of both Bs map to the same scaffolds, we should be able to find at least some loci that are heterozygous (hopefully). That would help us confirm the two Bs are actually homologous chromosomes and to certain extend quantify their divergence (although we need to understand how crute the estimate will be!). I will use [these assignments](output/scaffolds.final.assignment.table.csv).

```bash
# /data/ross/mealybugs/analyses/B_viburni_2020/5_B_char/heterozygosity_in_B04

# ln -s /data/ross/mealybugs/analyses/B_viburni_2020/2_short_read_DNA_seq/1_mapping/PV_18-13.Illumina.350.sorted.bam* .
# ln -s /data/ross/mealybugs/analyses/B_viburni_2020/2_short_read_DNA_seq/1_mapping/PV_18-13.Illumina.550.sorted.bam* .

# Creating lists of scaffolds from scaffolds.final.assignment.table.csv
grep "B1" scaffolds.final.assignment.table.csv | cut -f 4 -d '"' > B1_scaffolds.list
grep "B2" scaffolds.final.assignment.table.csv | cut -f 4 -d '"' > B2_scaffolds.list
grep "B3" scaffolds.final.assignment.table.csv | cut -f 4 -d '"' > B3_scaffolds.list

qsub -o logs -e logs -cwd -N bam2sync -V -pe smp64 16 -b yes 'samtools mpileup -a --no-BAQ --fasta-ref p.viburni.freeze.v0.fa --output /scratch/$USER/PV_18-13.Illumina.550.mpileup PV_18-13.Illumina.550.sorted.bam && java -jar ~/src/popoolation2/mpileup2sync.jar --input /scratch/$USER/PV_18-13.Illumina.550.mpileup --threads 16 --output /scratch/$USER/PV_18-13.Illumina.550.mpileup.sync && rsync -av --remove-source-files /scratch/$USER/PV_18-13.Illumina.550.mpileup.sync .'

qsub -o logs -e logs -cwd -N bam2sync -V -pe smp64 16 -b yes 'samtools mpileup -a --no-BAQ --fasta-ref p.viburni.freeze.v0.fa --output /scratch/$USER/PV_18-13.Illumina.350.mpileup PV_18-13.Illumina.350.sorted.bam && java -jar ~/src/popoolation2/mpileup2sync.jar --input /scratch/$USER/PV_18-13.Illumina.350.mpileup --threads 16 --output /scratch/$USER/PV_18-13.Illumina.350.mpileup.sync && rsync -av --remove-source-files /scratch/$USER/PV_18-13.Illumina.350.mpileup.sync .'

qsub -o logs -e logs -cwd -N bam2sync -V -pe smp64 16 -b yes 'samtools mpileup -a --no-BAQ --fasta-ref p.viburni.freeze.v0.fa --output /scratch/$USER/18-04.freeze.v0.mpileup PV_18-04.freeze.v0.sorted.bam && java -jar ~/src/popoolation2/mpileup2sync.jar --input /scratch/$USER/18-04.freeze.v0.mpileup --threads 16 --output /scratch/$USER/PV_18-04.freeze.v0.mpileup.sync && rsync -av --remove-source-files /scratch/$USER/PV_18-04.freeze.v0.mpileup.sync .'
```

The `sync` files contain for each genomic position coverage supports for each of 4 nucleotides. Most of them will have support of only one, but we are not interested in those. Let's select only those positions with two states where the less covered one has at least 3x coverage.

```
cat PV_18-13.Illumina.550.mpileup.sync | ./sync2variable_sites.py > PV_18-13.Illumina.550_bistates.tsv
# discarded in total 978040 sites

cat PV_18-13.Illumina.350.mpileup.sync | ./sync2variable_sites.py > PV_18-13.Illumina.350_bistates.tsv
# discarded in total 821457 sites

cat PV_18-04.freeze.v0.mpileup.sync | ./sync2variable_sites.py > PV_18-04.freeze.v0_bistates.tsv
# discarded in total 1765812 sites
```

Notice that `PV_18-04` had a lot more tri/tetra states. They still might be relevant for Bs are actually really repetitive. So, let's not forget we might have deiscarded a lot of signal here.

```
cat PV_18-04.freeze.v0_bistates.tsv | ./annotate_bistates.py > PV_18-04.freeze.v0_B_bistates.tsv
cat PV_18-13.Illumina.350_bistates.tsv | ./annotate_bistates.py > PV_18-13.Illumina.350_B_bistates.tsv
cat PV_18-13.Illumina.550_bistates.tsv | ./annotate_bistates.py > PV_18-13.Illumina.550_B_bistates.tsv
```

Note, there are so many multi-states for Bs, perhaps I should do the "annotate_bistates" step on the raw sync files

### Plotting the distributions in R

For this section you need to install in R [smudgeplot](https://github.com/KamilSJaron/smudgeplot) library. That is not so difficult, you can either install the whole package via conda, but that is unnecesarily headvy, we use only the plotting part of the program, which is very easy to install on its own. If there is nothing unusual on your computational settings, this should work: cloning the smudgeplot repo, `cd smudgeplot`, open `R` and run `install.packages(".", repos = NULL, type="source")`.

Once that is sone we can explore our data.

```{R}
library(smudgeplot)

PV04 <- read.table('PV_18-04.freeze.v0_B_bistates.tsv', col.names = c('scf', 'pos', 'covA', 'covB', 'chr'))
PV13_350 <- read.table('PV_18-13.Illumina.350_B_bistates.tsv', col.names = c('scf', 'pos', 'covA', 'covB', 'chr'))
PV13_550 <- read.table('PV_18-13.Illumina.550_B_bistates.tsv', col.names = c('scf', 'pos', 'covA', 'covB', 'chr'))

PV04_Autosome <- read.table('PV_18-04.freeze.v0_Autosome_sample.tsv', col.names = c('scf', 'pos', 'covA', 'covB'))
PV13_350_Autosome <- read.table('PV_18-13.Illumina.350_bistates_Autosome_sample.tsv', col.names = c('scf', 'pos', 'covA', 'covB'))
PV13_550_Autosome <- read.table('PV_18-13.Illumina.550_bistates_Autosome_sample.tsv', col.names = c('scf', 'pos', 'covA', 'covB'))

# this is for smudgeplot-like plot, extracted from smudgeplot codebase (https://github.com/KamilSJaron/smudgeplot)
smudgelike_plot <- function(minor_variant_rel_cov, total_pair_cov, ymax = 250, nbins = 30, draft_n = 50){
	smudge_container <- get_smudge_container(minor_variant_rel_cov, total_pair_cov,
																					 .nbins = nbins, .ylim = c(0, ymax))
	colour_ramp <- get_default_col_ramp() # get the default colour ramp (Spectral, 11)
	plot_smudgeplot(smudge_container, draft_n, colour_ramp)
	return(smudge_container)
}

# this is to extract the global maximum using kernel smoothing
kernel_smoothing_get_peak <- function(coverages, which_peak = 1, adjust = 1){
	ks <- density(coverages, bw = "SJ", adjust = adjust)
  second_deriv <- diff(sign(diff(ks$y)))

  peak_covs <- ks$x[which(second_deriv == -2) + 1]
  peak_heights <- ks$y[which(second_deriv == -2) + 1]
  peak_covs[order(peak_heights, decreasing=T)][which_peak]
}

PV04_minor_variant_rel_cov <- PV04$covB / (PV04$covA + PV04$covB)
PV04_total_pair_cov <- PV04$covA + PV04$covB
PV04_minor_variant_rel_cov_A <- PV04_Autosome$covB / (PV04_Autosome$covA + PV04_Autosome$covB)
PV04_total_pair_cov_A <- PV04_Autosome$covA + PV04_Autosome$covB

pdf('allelic_smudgeplot-like-plots.pdf', width = 16, height = 10)

par(mfrow = c(3, 4))

############
### PV04 ###
############

# plot the individual smudgeplots, the containers will contain the plotted matrix (we don't do anything with is atm, I explored it a bit, see commented code)
PV04_A_container <- smudgelike_plot(PV04_minor_variant_rel_cov_A, PV04_total_pair_cov_A, ymax = 250, nbin = 30, draft_n = 48)
PV04_B1_container <- smudgelike_plot(PV04_minor_variant_rel_cov[PV04$chr == 'B1'], PV04_total_pair_cov[PV04$chr == 'B1'], ymax = 250, nbin = 30, draft_n = 48)
PV04_B2_container <- smudgelike_plot(PV04_minor_variant_rel_cov[PV04$chr == 'B2'], PV04_total_pair_cov[PV04$chr == 'B2'], ymax = 250, nbin = 30, draft_n = 48)
PV04_B3_container <- smudgelike_plot(PV04_minor_variant_rel_cov[PV04$chr == 'B3'], PV04_total_pair_cov[PV04$chr == 'B3'], ymax = 250, nbin = 30, draft_n = 48)

# I remove the really messy parts of the graph and use kernel smoothing to estimate
PV04_cov_sums_A <- PV04_total_pair_cov_A[PV04_minor_variant_rel_cov_A > 0.25 & PV04_total_pair_cov_A < 250]
PV04_sane_filter <- PV04_minor_variant_rel_cov > 0.25 & PV04_total_pair_cov < 250
PV04_cov_sums_B1 <- PV04_total_pair_cov[PV04_sane_filter & PV04$ch == 'B2']
PV04_cov_sums_B2 <- PV04_total_pair_cov[PV04_sane_filter & PV04$ch == 'B3']

sapply(list(PV04_cov_sums_A, PV04_cov_sums_B1, PV04_cov_sums_B2), kernel_smoothing_get_peak) # returns a vector of 2n coverage estimates from A, B2 and B3 subsets
PV04_1n_cov_est <- 48.87182 # I adjusted "draft_n" above with this value and replotted the image

###############
##### PV13 ####
# library 350 #
###############

PV13_350_minor_variant_rel_cov <- PV13_350$covB / (PV13_350$covA + PV13_350$covB)
PV13_350_total_pair_cov <- PV13_350$covA + PV13_350$covB
PV13_350_minor_variant_rel_cov_A <- PV13_350_Autosome$covB / (PV13_350_Autosome$covA + PV13_350_Autosome$covB)
PV13_350_total_pair_cov_A <- PV13_350_Autosome$covA + PV13_350_Autosome$covB

PV13_350_A_container <- smudgelike_plot(PV13_350_minor_variant_rel_cov_A, PV13_350_total_pair_cov_A, draft_n = 26, ymax = 125)
PV13_350_B1_container <- smudgelike_plot(PV13_350_minor_variant_rel_cov[PV13_350$chr == 'B1'], PV13_350_total_pair_cov[PV13_350$chr == 'B1'], draft_n = 26, ymax = 125)
PV13_350_B2_container <- smudgelike_plot(PV13_350_minor_variant_rel_cov[PV13_350$chr == 'B2'], PV13_350_total_pair_cov[PV13_350$chr == 'B2'], draft_n = 26, ymax = 125)
PV13_350_B3_container <- smudgelike_plot(PV13_350_minor_variant_rel_cov[PV13_350$chr == 'B3'], PV13_350_total_pair_cov[PV13_350$chr == 'B3'], draft_n = 26, ymax = 125)

PV13_350_cov_sums_A <- PV13_350_total_pair_cov_A[PV13_350_minor_variant_rel_cov_A > 0.25 & PV13_350_total_pair_cov_A < 250]
PV13_350_sane_filter <- PV13_350_minor_variant_rel_cov > 0.25 & PV13_350_total_pair_cov < 250
PV13_350_cov_sums_B2 <- PV13_350_total_pair_cov[PV13_350_sane_filter & PV13_350$ch == 'B2']
PV13_350_cov_sums_B3 <- PV13_350_total_pair_cov[PV13_350_sane_filter & PV13_350$ch == 'B3']
sapply(list(PV13_350_cov_sums_A, PV13_350_cov_sums_B2, PV13_350_cov_sums_B3), kernel_smoothing_get_peak) / 2
# [1] 25.92866 25.57749 26.83936
# against the expecation in this line it seems that the B2 and B3 seqeunces have the same coverage as the autosome

###############
##### PV13 ####
# library 550 #
###############

# transformation of coverages of the two states
PV13_550_minor_variant_rel_cov <- PV13_550$covB / (PV13_550$covA + PV13_550$covB)
PV13_550_total_pair_cov <- PV13_550$covA + PV13_550$covB
PV13_550_minor_variant_rel_cov_A <- PV13_550_Autosome$covB / (PV13_550_Autosome$covA + PV13_550_Autosome$covB)
PV13_550_total_pair_cov_A <- PV13_550_Autosome$covA + PV13_550_Autosome$covB

# plotting
PV13_550_A_container <- smudgelike_plot(PV13_550_minor_variant_rel_cov_A, PV13_550_total_pair_cov_A, draft_n = 19, ymax = 125)
PV13_550_B1_container <- smudgelike_plot(PV13_550_minor_variant_rel_cov[PV13_550$chr == 'B1'], PV13_550_total_pair_cov[PV13_550$chr == 'B1'], draft_n = 19, ymax = 125, nbin = 25)
PV13_550_B2_container <- smudgelike_plot(PV13_550_minor_variant_rel_cov[PV13_550$chr == 'B2'], PV13_550_total_pair_cov[PV13_550$chr == 'B2'], draft_n = 19, ymax = 125, nbin = 25)
PV13_550_B3_container <- smudgelike_plot(PV13_550_minor_variant_rel_cov[PV13_550$chr == 'B3'], PV13_550_total_pair_cov[PV13_550$chr == 'B3'], draft_n = 19, ymax = 125, nbin = 25)

# coverage estimate
PV13_550_cov_sums_A <- PV13_550_total_pair_cov_A[PV13_550_minor_variant_rel_cov_A > 0.25 & PV13_550_total_pair_cov_A < 250]
PV13_550_sane_filter <- PV13_550_minor_variant_rel_cov > 0.25 & PV13_550_total_pair_cov < 250
PV13_550_cov_sums_B2 <- PV13_550_total_pair_cov[PV13_550_sane_filter & PV13_550$ch == 'B2']
PV13_550_cov_sums_B3 <- PV13_550_total_pair_cov[PV13_550_sane_filter & PV13_550$ch == 'B3']
sapply(list(PV13_550_cov_sums_A, PV13_550_cov_sums_B2, PV13_550_cov_sums_B3), kernel_smoothing_get_peak) / 2

dev.off()
```

There are two more sanity checks to be made. One to check the consistency of the individual B-linked scaffolds in PV13 line. And second, to check consistency between 350 and 550 libraries. Are the same sites "heterozygous"?

```{R}
PV04$cov_sum <- PV04_total_pair_cov
PV04$cov_ratio <- PV04_minor_variant_rel_cov
PV04$sane <- PV04_sane_filter

scfs_tab <- data.frame(scf = unique(c(PV04$scf, PV13_350$scf, PV13_550$scf)))
rownames(scfs_tab) <- scfs_tab$scf

sites_per_scf <- table(PV04$scf)
fraction_of_sane <- sapply(names(sites_per_scf), function(scf){ mean(PV04[PV04$scf == scf, 'sane']) })

scfs_tab[names(sites_per_scf), 'PV04_bistates'] <- as.numeric(sites_per_scf)
scfs_tab[names(fraction_of_sane), 'PV04_sane_fraction'] <- as.numeric(fraction_of_sane)

PV13_350$sane <- PV13_350_sane_filter
sites_per_scf <- table(PV13_350$scf)
fraction_of_sane <- sapply(names(sites_per_scf), function(scf){ mean(PV13_350[PV13_350$scf == scf, 'sane']) })

scfs_tab[names(sites_per_scf), 'PV13_350_bistates'] <- as.numeric(sites_per_scf)
scfs_tab[names(fraction_of_sane), 'PV13_350_sane_fraction'] <- as.numeric(fraction_of_sane)

PV13_550$sane <- PV13_550_sane_filter
sites_per_scf <- table(PV13_550$scf)
fraction_of_sane <- sapply(names(sites_per_scf), function(scf){ mean(PV13_550[PV13_550$scf == scf, 'sane']) })

scfs_tab[names(sites_per_scf), 'PV13_550_bistates'] <- as.numeric(sites_per_scf)
scfs_tab[names(fraction_of_sane), 'PV13_550_sane_fraction'] <- as.numeric(fraction_of_sane)

scfs_tab[is.na(scfs_tab)] <- 0
# if there are no bistates, it should be a 0, not NA
```

And now I will just add there the scaffold lengths using `p.viburni.freeze.v0_scf_sizes.tsv` table of scaffold sizes.

```{R}
B1_col <- 'dodgerblue4'
B2_col <- 'deepskyblue1'
B3_col <- 'darkcyan'

scf_lengths <- read.table('p.viburni.freeze.v0_scf_sizes.tsv', col.names = c('scf', 'len'))
rownames(scf_lengths) <- scf_lengths$scf

scfs_tab$length <- scf_lengths[scfs_tab$scf, 'len']

PV13_350_sane_per_base <- scfs_tab$PV13_350_bistates * scfs_tab$PV13_350_sane_fraction / scfs_tab$length
PV13_550_sane_per_base <- scfs_tab$PV13_550_bistates * scfs_tab$PV13_550_sane_fraction / scfs_tab$length
PV04_sane_per_base <- scfs_tab$PV04_bistates * scfs_tab$PV04_sane_fraction / scfs_tab$length

plot(PV13_350_sane_per_base ~ scfs_tab$length)
plot(PV13_550_sane_per_base ~ scfs_tab$length)

plot(PV13_350_sane_per_base ~ PV13_550_sane_per_base)
# THIS shows that 350 and 550 libraries are actually really consistent in the "amount of detected heterozygosity", I suppose they would be probably every consistent regarding the positions too!

chormosome_assignments <- rbind(PV13_550[, c('scf', 'chr')], PV13_350[, c('scf', 'chr')], PV04[, c('scf', 'chr')])
chormosome_assignments <- chormosome_assignments[!duplicated(chormosome_assignments), ]
rownames(chormosome_assignments) <- chormosome_assignments$scf

scfs_tab$asn <- chormosome_assignments[scfs_tab$scf, 'chr']

plot(PV04_sane_per_base[scfs_tab$length > 20e3], PV13_350_sane_per_base[scfs_tab$length > 20e3], xlim = c(0, 0.02), ylim = c(0, 0.02), xlab = 'PV04 SNP-candidate density', ylab = 'PV13 SNP-candidate density')
points(PV04_sane_per_base[scfs_tab$asn == 'B1'], PV13_350_sane_per_base[scfs_tab$asn == 'B1'], pch = 20, col = B1_col)
points(PV04_sane_per_base[scfs_tab$asn == 'B2'], PV13_350_sane_per_base[scfs_tab$asn == 'B2'], pch = 20, col = B2_col)
points(PV04_sane_per_base[scfs_tab$asn == 'B3'], PV13_350_sane_per_base[scfs_tab$asn == 'B3'], pch = 20, col = B3_col)
lines(c(0, 0.019), c(0, 0.019))

legend('topright', c('B1', 'B2', 'B3'), col = c('orange', 'cyan', 'magenta'), pch = 20, bty = 'n')
```
