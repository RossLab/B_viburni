### Outlier analysis

Vast majority of B genes are either with very low expression in B- lines, or present on scaffolds that happen to be duplicated between the core genome and the B chromosome. All the works out with one exception: `g2644`, located between 25304 and 32976 base of `scaffold_786` (total 48770 nt).

The weird outlier is on a scaffold that is certainly B-linked but a bit iffy.
The normalised gneomic coverages (i.e. estimated number of genomic copies) of its scaffold are 10.01, 4.36,0.67, 0.72 in 04 (BBB), 13 (B), 21 (-), 23 (-). The mean coverage in B- lines is too high for "nothing", but too low to be reliably existent. I suspect there will be a segment of that scaffold that will be again on the core genome, but not all, that's why the coverage is somewhere in between. And if that's the case, the gene would be in the shared location. We will verify by looking at the coverage distribution of that scaffold.

### Extracting coverage of the funky scaffold from all 4 lines

It's `scaffold_786`, se let's get per-base coverages

```bash
# cd /data/ross/mealybugs/analyses/B_viburni_2020/5_B_char/heterozygosity_in_B04
for strain in PV_18-04 PV_18-13 PV_18-21 PV_18-23; do
	qsub -o logs -e logs -cwd -N dep_s786 -V -pe smp64 3 -b yes "samtools view -h $strain.freeze.v0.sorted.bam 'scaffold_786' | samtools depth -aa - | grep '^scaffold_786' > /scratch/kjaron/"$strain"_scaffold_786.depth && rsync -av --remove-source-files /scratch/$USER/"$strain"_scaffold_786.depth funky_scaffolds/"
done
```

### Plotting genomic coverage of the funky scaffold

```{R}
strains <- c("04", "13", "21", "23")
depth_files <- paste0("data/4_cov_analysis/cov/funky_scaffolds/PV_18-", strains ,"_scaffold_786.depth")

depth_tabs <- lapply(depth_files, read.table, col.names = c('scf', 'pos', 'cov'))

monoploid_coverages <- c(53.28300, 51.58885, 28.51315, 30.68950)
normalised_depths <- lapply(1:4, function(i){depth_tabs[[i]][, 3] / monoploid_coverages[i]})
names(normalised_depths) <- strains
xrange <- 1:48770

png('output/funky_gene_coverage_plot.png')

plot(normalised_depths[["04"]], type = 'l', lwd = 3, col = 'dodgerblue3', ylim = c(0, 6), xlim = c(20000, 35000), xlab = 'position on scaffold_786', ylab = 'normalised genomic coverage')
lines(normalised_depths[["13"]], lwd = 3, col = 'limegreen')
lines(normalised_depths[["21"]], lwd = 3, col = 'deeppink')
lines(normalised_depths[["23"]], type = 'l', lwd = 3, col = 'goldenrod2')

g2644_exons <- read.table('data/4_cov_analysis/cov/funky_scaffolds/g2644_exon_coordinates.tsv', header = F, col.names = c('from', 'to'))

lines(c(25304, 32976), c(-0.1, -0.1), lwd = 5)
for (i in 1:nrow(g2644_exons)){
	lines(c(g2644_exons[i,1], g2644_exons[i,2]), c(-0.1, -0.1), lwd = 3, col = 'red')
}

legend('topleft', bty = 'n', col = c('dodgerblue3', 'limegreen', 'deeppink', 'goldenrod2'), lty = 1, lwd = 3, c("04 - BBB", "13 - B", "21 - no B", "23 - no B"), title = 'Normalized genomic coverage')

legend('topright', bty = 'n', col = c('black', 'red'), lty = 1, lwd = c(5, 3), c('whole gene', 'exons'), title = 'g2644')

# ?? TODO add RNA-seq mapping (so we don't have it saved)
dev.off()
```

### Getting expression per replicate

```bash
# /data/ross/mealybugs/analyses/B_viburni_2020/3_RNA_seq/4_genome_based/RSEM_results

echo -n "File:" > expression_of_the_funky_gene/g2644_expression.tsv
head -1 04F_1.genes.results >> expression_of_the_funky_gene/g2644_expression.tsv
grep "g2644" *.genes.results >> expression_of_the_funky_gene/g2644_expression.tsv
```

(note the location within directory has changed. I am not trying to reorganise the cluster, all commited files are the same tough)

```{R}
expression_tab <- read.table('data/4_cov_analysis/cov/funky_scaffolds/g2644_expression.tsv', header = T)

expressions <- expression_tab$TPM
names(expressions) <- sapply(strsplit(expression_tab[, 1], "[.]"), function(x){ x[1]} )

pal <- c(rep('dodgerblue3', 6), rep('limegreen', 7), rep('orangered3', 6), rep('deeppink', 7))
barplot(expressions, ylab = 'TPMs', main = 'Expression of g2644, the funky gene', col = pal)

legend('topleft', c("04 - BBB", "13 - B", "15 - no B", "21 - no B"), col = c('dodgerblue3', 'limegreen', 'orangered3', 'deeppink'), bty = 'n', pch = 20, cex = 1.3)
```

### More detailed expression data

```bash
conda activate /ceph/users/kjaron/../afilia/.conda/envs/afilia_trinity
mkdir RSEM_results_bams

qsub -o logs -e logs -cwd -N RSEM_ref -V -pe smp64 32 -b yes 'rsem-prepare-reference  --gff3 /data/ross/mealybugs/analyses/B_viburni_2020/1_pacbio_assembly/8_freeze_v0/p.viburni.freeze.v0.braker.gff3 --star -p 32 /data/ross/mealybugs/analyses/B_viburni_2020/1_pacbio_assembly/8_freeze_v0/p.viburni.freeze.v0.fa p.viburni.freeze.v0.'

qsub -o logs -e logs -cwd -N RSEM04F_1 -V -pe smp64 32 -b yes 'rsem-calculate-expression -p 32 --paired-end --star-gzipped-read-file --strandedness reverse --calc-ci --calc-pme --star ../0_reads/04F_1.trimmed_1.fastq.gz ../0_reads/04F_1.trimmed_2.fastq.gz p.viburni.freeze.v0. RSEM_results_bams/04F_1'
qsub -o logs -e logs -cwd -N RSEM04F_2 -V -pe smp64 32 -b yes 'rsem-calculate-expression -p 32 --paired-end --star-gzipped-read-file --strandedness reverse --calc-ci --calc-pme --star ../0_reads/04F_2.trimmed_1.fastq.gz ../0_reads/04F_2.trimmed_2.fastq.gz p.viburni.freeze.v0. RSEM_results_bams/04F_2'
qsub -o logs -e logs -cwd -N RSEM04M_2 -V -pe smp64 32 -b yes 'rsem-calculate-expression -p 32 --paired-end --star-gzipped-read-file --strandedness reverse --calc-ci --calc-pme --star ../0_reads/04M_2.trimmed_1.fastq.gz ../0_reads/04M_2.trimmed_2.fastq.gz p.viburni.freeze.v0. RSEM_results_bams/04M_2'

qsub -o logs -e logs -cwd -N RSEM15M_1 -V -pe smp64 8 -b yes 'mkdir -p /scratch/$USER/RSEM_15M_1 && rsem-calculate-expression -p 8 --paired-end --star-gzipped-read-file --strandedness reverse --calc-ci --calc-pme --star ../0_reads/15M_1.trimmed_1.fastq.gz ../0_reads/15M_1.trimmed_2.fastq.gz p.viburni.freeze.v0. /scratch/$USER/RSEM_15M_1/15M_1 && rm -r /scratch/$USER/RSEM_15M_1/15M_1.temp && rsync -av --remove-source-files /scratch/$USER/RSEM_15M_1 RSEM_results_bams/'
qsub -o logs -e logs -cwd -N RSEM15F_3 -V -pe smp64 8 -b yes 'mkdir -p /scratch/$USER/RSEM_15F_3 && rsem-calculate-expression -p 8 --paired-end --star-gzipped-read-file --strandedness reverse --calc-ci --calc-pme --star ../0_reads/15F_3.trimmed_1.fastq.gz ../0_reads/15F_3.trimmed_2.fastq.gz p.viburni.freeze.v0. /scratch/$USER/RSEM_15F_3/15F_3 && rm -r /scratch/$USER/RSEM_15F_3/15F_3.temp && rsync -av --remove-source-files /scratch/$USER/RSEM_15F_3 RSEM_results_bams/'
qsub -o logs -e logs -cwd -N RSEM15M_3 -V -pe smp64 8 -b yes 'mkdir -p /scratch/$USER/RSEM_15M_3 && rsem-calculate-expression -p 8 --paired-end --star-gzipped-read-file --strandedness reverse --calc-ci --calc-pme --star ../0_reads/15M_3.trimmed_1.fastq.gz ../0_reads/15M_3.trimmed_2.fastq.gz p.viburni.freeze.v0. /scratch/$USER/RSEM_15M_3/15M_3 && rm -r /scratch/$USER/RSEM_15M_3/15M_3.temp && rsync -av --remove-source-files /scratch/$USER/RSEM_15M_3 RSEM_results_bams/'

TODO: adjust the commands underneath
qsub -o logs -e logs -cwd -N RSEM21M_3 -V -pe smp64 8 -b yes 'rsem-calculate-expression -p 8 --paired-end --star-gzipped-read-file --strandedness reverse --calc-ci --calc-pme --star ../0_reads/21M_3.trimmed_1.fastq.gz ../0_reads/21M_3.trimmed_2.fastq.gz p.viburni.freeze.v0. RSEM_results_bams/21M_3'
qsub -o logs -e logs -cwd -N RSEM21M_4 -V -pe smp64 8 -b yes 'rsem-calculate-expression -p 8 --paired-end --star-gzipped-read-file --strandedness reverse --calc-ci --calc-pme --star ../0_reads/21M_4.trimmed_1.fastq.gz ../0_reads/21M_4.trimmed_2.fastq.gz p.viburni.freeze.v0. RSEM_results_bams/21M_4'
qsub -o logs -e logs -cwd -N RSEM21F_3 -V -pe smp64 8 -b yes 'rsem-calculate-expression -p 8 --paired-end --star-gzipped-read-file --strandedness reverse --calc-ci --calc-pme --star ../0_reads/21F_3.trimmed_1.fastq.gz ../0_reads/21F_3.trimmed_2.fastq.gz p.viburni.freeze.v0. RSEM_results_bams/21F_3'
```

---

```bash
cd RSEM_results_bams

qsub -o logs -e logs -cwd -N dep_s786 -V -pe smp64 3 -b yes "samtools sort 04F_1.transcript.bam > /scratch/kjaron/04F_1.transcript.sorted.bam && samtools index /scratch/kjaron/04F_1.transcript.sorted.bam && samtools view -h /scratch/kjaron/04F_1.transcript.sorted.bam 'g2644.t1' | samtools depth -aa - | grep '^g2644.t1' > funky_gene/04F_1_g2644.t1.depth"
qsub -o logs -e logs -cwd -N dep_04F_2 -V -pe smp64 3 -b yes "samtools sort 04F_2.transcript.bam > /scratch/kjaron/04F_2.transcript.sorted.bam && samtools index /scratch/kjaron/04F_2.transcript.sorted.bam && samtools view -h /scratch/kjaron/04F_2.transcript.sorted.bam 'g2644.t1' | samtools depth -aa - | grep '^g2644.t1' > funky_gene/04F_2_g2644.t1.depth"
qsub -o logs -e logs -cwd -N dep_04M_2 -V -pe smp64 3 -b yes "samtools sort 04M_2.transcript.bam > /scratch/kjaron/04M_2.transcript.sorted.bam && samtools index /scratch/kjaron/04M_2.transcript.sorted.bam && samtools view -h /scratch/kjaron/04M_2.transcript.sorted.bam 'g2644.t1' | samtools depth -aa - | grep '^g2644.t1' > funky_gene/04M_2_g2644.t1.depth"

qsub -o logs -e logs -cwd -N dep_s786 -V -pe smp64 3 -b yes "samtools sort RSEM_15M_1/15M_1.transcript.bam > /scratch/kjaron/15M_1.transcript.sorted.bam && samtools index /scratch/kjaron/15M_1.transcript.sorted.bam && samtools view -h /scratch/kjaron/15M_1.transcript.sorted.bam 'g2644.t1' | samtools depth -aa - | grep '^g2644.t1' > funky_gene/15M_1_g2644.t1.depth"
```

### Some notes about the gene function

The funky gene, `g2644`, is a homolog of `PGBD4_HUMAN`, which is a `PiggyBac transposable element-derived protein 4`, that contains Integrase zinc-binding domain - could be either contamination (unlikely in all replicates), or simply a bastard TE.
