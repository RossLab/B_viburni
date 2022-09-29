### B heterozygosity of scaffold_360 and scaffold_957

Two scaffolds (`scaffold_360` and `scaffold_957`) are present on both autosomes and Bs. Namely, PV13 has 4 copies (2 autosomal, 2 B-copies) and PV04 has 8 copies (2 autosomal, 6 B-copies) of `scaffold_360` and `scaffold_957` should have 9 (8A, 1B) and 11 (8A, 3B) in PV13 and PV04 respectivelly.

Furthermore, we might also detect some loci that are heterozygous within B and A. I will use [these assignments](output/scaffolds.final.assignment.table.csv).

Following block of code will create pileup files in human redable format (`sync`). Which will allow us then to subtrack these two intereting scaffolds.

```bash
# /data/ross/mealybugs/analyses/B_viburni_2020/5_B_char/heterozygosity_in_B04

# ln -s /data/ross/mealybugs/analyses/B_viburni_2020/2_short_read_DNA_seq/1_mapping/PV_18-13.Illumina.350.sorted.bam* .
# ln -s /data/ross/mealybugs/analyses/B_viburni_2020/2_short_read_DNA_seq/1_mapping/PV_18-13.Illumina.550.sorted.bam* .

qsub -o logs -e logs -cwd -N bam2sync -V -pe smp64 16 -b yes 'samtools mpileup -a --no-BAQ --fasta-ref p.viburni.freeze.v0.fa --output /scratch/$USER/PV_18-13.Illumina.550.mpileup PV_18-13.Illumina.550.sorted.bam && java -jar ~/src/popoolation2/mpileup2sync.jar --input /scratch/$USER/PV_18-13.Illumina.550.mpileup --threads 16 --output /scratch/$USER/PV_18-13.Illumina.550.mpileup.sync && rsync -av --remove-source-files /scratch/$USER/PV_18-13.Illumina.550.mpileup.sync .'

qsub -o logs -e logs -cwd -N bam2sync -V -pe smp64 16 -b yes 'samtools mpileup -a --no-BAQ --fasta-ref p.viburni.freeze.v0.fa --output /scratch/$USER/PV_18-13.Illumina.350.mpileup PV_18-13.Illumina.350.sorted.bam && java -jar ~/src/popoolation2/mpileup2sync.jar --input /scratch/$USER/PV_18-13.Illumina.350.mpileup --threads 16 --output /scratch/$USER/PV_18-13.Illumina.350.mpileup.sync && rsync -av --remove-source-files /scratch/$USER/PV_18-13.Illumina.350.mpileup.sync .'

qsub -o logs -e logs -cwd -N bam2sync -V -pe smp64 16 -b yes 'samtools mpileup -a --no-BAQ --fasta-ref p.viburni.freeze.v0.fa --output /scratch/$USER/18-04.freeze.v0.mpileup PV_18-04.freeze.v0.sorted.bam && java -jar ~/src/popoolation2/mpileup2sync.jar --input /scratch/$USER/18-04.freeze.v0.mpileup --threads 16 --output /scratch/$USER/PV_18-04.freeze.v0.mpileup.sync && rsync -av --remove-source-files /scratch/$USER/PV_18-04.freeze.v0.mpileup.sync .'
```

The `sync` files contain for each genomic position coverage supports for each of 4 nucleotides.

```bash
grep -w "^scaffold_360" PV_18-04.freeze.v0.mpileup.sync > PV_18-04.scaffold_360.mpileup.sync
grep -w "^scaffold_360" PV_18-13.Illumina.350.mpileup.sync > PV_18-13_350.scaffold_360.mpileup.sync
grep -w "^scaffold_360" PV_18-13.Illumina.550.mpileup.sync > PV_18-13_550.scaffold_360.mpileup.sync
```

or possibly using only bistates (`PV_18-13.Illumina.550_bistates.tsv`)

```bash
cat PV_18-04.scaffold_360.mpileup.sync | python3 sync2multiallelic.py > PV_18-04.scaffold_360.mpileup.multiallelic.tsv
cat PV_18-13_350.scaffold_360.mpileup.sync | python3 sync2multiallelic.py > PV_18-13_350.scaffold_360.mpileup.multiallelic.tsv
cat PV_18-13_550.scaffold_360.mpileup.sync | python3 sync2multiallelic.py > PV_18-13_550.scaffold_360.mpileup.multiallelic.tsv
```

Now that I have three files with all multistates (i.e. more than one nucleotide mapping at the position), we can plot it.

### Plotting the distributions in R

For this section you need to install in R [smudgeplot](https://github.com/KamilSJaron/smudgeplot) library. That is not so difficult, you can either install the whole package via conda, but that is unnecesarily headvy, we use only the plotting part of the program, which is very easy to install on its own. If there is nothing unusual on your computational settings, this should work: cloning the smudgeplot repo, `cd smudgeplot`, open `R` and run `install.packages(".", repos = NULL, type="source")`.

Once that is sone we can explore our data.

```{R}
library(smudgeplot)

PV04 <- read.table('PV_18-04.scaffold_360.mpileup.multiallelic.tsv', col.names = c('scf', 'pos', 'covA', 'covB'))
PV13_350 <- read.table('PV_18-13_350.scaffold_360.mpileup.multiallelic.tsv', col.names = c('scf', 'pos', 'covA', 'covB'))
PV13_550 <- read.table('PV_18-13_550.scaffold_360.mpileup.multiallelic.tsv', col.names = c('scf', 'pos', 'covA', 'covB'))

PV04 <- PV04[PV04$covB > 20, ]

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

# pdf('allelic_smudgeplot-like-plots.pdf', width = 16, height = 10)

par(mfrow = c(3, 4))

############
### PV04 ###
############

# plot the individual smudgeplots, the containers will contain the plotted matrix (we don't do anything with is atm, I explored it a bit, see commented code)
PV04_B_container <- smudgelike_plot(PV04_minor_variant_rel_cov, PV04_total_pair_cov, ymax = 1000, nbin = 50, draft_n = 53)

###############
##### PV13 ####
# library 350 #
###############

PV13_350_minor_variant_rel_cov <- PV13_350$covB / (PV13_350$covA + PV13_350$covB)
PV13_350_total_pair_cov <- PV13_350$covA + PV13_350$covB

PV13_350_B_container <- smudgelike_plot(PV13_350_minor_variant_rel_cov, PV13_350_total_pair_cov, draft_n = 26, ymax = 300)

###############
##### PV13 ####
# library 550 #
###############

# transformation of coverages of the two states
PV13_550_minor_variant_rel_cov <- PV13_550$covB / (PV13_550$covA + PV13_550$covB)
PV13_550_total_pair_cov <- PV13_550$covA + PV13_550$covB

# plotting
PV13_550_A_container <- smudgelike_plot(PV13_550_minor_variant_rel_cov, PV13_550_total_pair_cov, draft_n = 19, ymax = 250)
```

There seems to be no SNPs on the `scaffold_360` that would follow any meaningful coverage patters (most of this seems to be just a mapping noise). This would explain a, why the coverage is a bit inflated compared the expecation across the samples b, the translocation is very recent and not enough singal has accumulated since.
