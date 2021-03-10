# Origins and divergence of Bs

In this section we will try to ask:
 - How divergent are the B orthologs to the core genome? (approximated by sequence identity)
 - Whether the sequence identity to within-genome genes is higer than identity to orthologs in different species (solenopsis, longispinus)
 - Are the orthologs specific to a single core-genome chromosome?

### Input data

The root for these analses will be `/data/ross/mealybugs/analyses/B_viburni_2020/7_B_divergence_and_origin`, but I will need more data from all over the place. Here is the list relative to: `/data/ross/mealybugs/analyses/`.

**P. viburni**
  - Genome: `B_viburni_2020/1_pacbio_assembly/8_freeze_v0/p.viburni.freeze.v0.fa`
  - Annotation: `B_viburni_2020/1_pacbio_assembly/8_freeze_v0/p.viburni.freeze.v0.braker.gff3`

**P. solenopsis**

- Genome: `phenacoccus_solenopsis/HiC_genome_Phenacoccus_solenopsis.fasta`
- Annotation: `phenacoccus_solenopsis/HiC_genome.gff3`

Note, when quering the reference, suffixes must be removed, e.g. `>PSOL00002-TA` needs to be trimmed to `PSOL00002`.

**P. longispinus**

- Genome: `pseudococcus_longispinus/8_freeze_v1/pseudococcus_longispinus.v1.assembly.fa`
- Annotation: `pseudococcus_longispinus/8_freeze_v1/pseudococcus_longispinus.v1.braker.gff3`

To make it easier, I linked the files above to three directories: `viburni`, `longispinus`, and `solenopsis`. The links are callec `genome.fa`, and `annotation.gff3` respectivelly.

#### Data preprocessing

For `viburni` and `longispinus` I need to get transcripts first, so I do it for all three to have it done consistently using gffread.

```
qsub -o logs -e logs -cwd -N gffread -V -pe smp64 1 -b yes 'gffread viburni/annotation.gff3 -g viburni/genome.fa -w viburni/transcripts.fa -y viburni/proteins.faa'
qsub -o logs -e logs -cwd -N gffread -V -pe smp64 1 -b yes 'gffread longispinus/annotation.gff3 -g longispinus/genome.fa -w longispinus/transcripts.fa -y longispinus/proteins.faa'
qsub -o logs -e logs -cwd -N gffread -V -pe smp64 1 -b yes 'gffread solenopsis/annotation.gff3 -g solenopsis/genome.fa -w solenopsis/transcripts.fa -y solenopsis/proteins.faa'
```

I also generated `viburni/B_genes.tsv` file from [SUPPLFILE1.xlsx](https://github.com/RossLab/B_viburni/blob/master/manuscript/supplementary/SUPPLFILE1.xlsx). It's a simple list of gene names and B1/2/3. Non-B genes are not included.

I will use this list to filter B genes out of `viburni/transcripts.fa`. Following script

```python
from Bio import SeqIO
from collections import defaultdict

gene2asn = defaultdict(lambda: "A")
with open('viburni/B_genes.tsv', 'r') as gene_file:
    for line in gene_file:
        line = line.split('\t')
        gene = line[0]
        asn = line[1].rstrip('\n')
        gene2asn[gene] = asn

for seq in SeqIO.parse('viburni/transcripts.fa', "fasta"):
    gene,iso = seq.name.split('.')
    if gene2asn[gene] != 'A' and iso == 't1':
        print(">Pvib_" + gene + '_' + gene2asn[gene])
        print(seq.seq)
```

saved to `extracting_B_transcripts.py`, is then executed as foolows:

```
python3 extracting_B_transcripts.py > viburni/B_transcripts.fa
```

### Orthology

All the following analyses will be based on some sort of orthology inference.

#### BLAST

The simple approach would be reciprocal blast of

1. within _P. viburni_
2. _P. viburni B_ + _P. solenopsis_
3. _P. viburni B_ + _P. longispinus_

The co-linearity should be based on aa reciprocal blast, the divergences should be calculated on nt reciprocal blast. I will start with nt, the `aa`.

1. within _P. viburni_

```
GENES=viburni/transcripts.fa
makeblastdb -in $GENES -dbtype nucl

mkdir -p blastout
qsub -o logs -e logs -cwd -N selfblast -V -pe smp64 16 -b yes 'blastn -query '"$GENES"' -db '"$GENES"' -evalue 1e-10 -outfmt "6 std qlen slen" -num_threads 16 > blastout/viburni_all_vs_all.blast'
```

2. _P. viburni B_ + _P. solenopsis_

```
GENES=solenopsis/solenopsis_and_B_viburni_transcripts.fa
cat viburni/B_transcripts.fa > $GENES
cat solenopsis/transcripts.fa >> $GENES
makeblastdb -in $GENES -dbtype nucl
qsub -o logs -e logs -cwd -N selfblast -V -pe smp64 16 -b yes 'blastn -query '"$GENES"' -db '"$GENES"' -evalue 1e-10 -outfmt "6 std qlen slen" -num_threads 16 > blastout/solenopsis_and_B_viburni_all_vs_all.blast'
```

3. _P. viburni B_ + _P. longispinus_

```
GENES=longispinus/longispinus_and_B_viburni_transcripts.fa
cat viburni/B_transcripts.fa > $GENES
cat longispinus/transcripts.fa >> $GENES
makeblastdb -in $GENES -dbtype nucl
qsub -o logs -e logs -cwd -N selfblast -V -pe smp64 16 -b yes 'blastn -query '"$GENES"' -db '"$GENES"' -evalue 1e-10 -outfmt "6 std qlen slen" -num_threads 16 > blastout/longispinus_and_B_viburni_all_vs_all.blast'
```

And using `reciprocal_blast.py` (originally written for Hodson et al. 2021) get the reciprocal hits out of all-to-all blast.

```
python3 reciprocal_blast.py blastout/viburni_all_vs_all.blast blastout/Within_viburni_reciprocal
python3 reciprocal_blast.py blastout/solenopsis_and_B_viburni_all_vs_all.blast blastout/solenopsis_and_B_viburni_reciprocal
python3 reciprocal_blast.py blastout/longispinus_and_B_viburni_all_vs_all.blast blastout/longispinus_and_B_viburni_reciprocal
```

#### OrthoMCL

Perhaps a more formal Orthology search would not hurt. I will just try to run a simple OrthoMCL run on all three genomes simuntaneously.

### Divergence analysis

```R
viburni_B_genes_file <- 'viburni/B_genes.tsv'
within_viburni_file <- 'blastout/Within_viburni_reciprocal_OG_pairs.tsv'
longispinus_file <- 'blastout/longispinus_to_B_viburni_reciprocal_OG_pairs.tsv'
solenopsis_file <- 'blastout/solenopsis_to_B_viburni_reciprocal_OG_pairs.tsv'

B_genes_tab <- read.table(viburni_B_genes_file, col.names = c('gene', 'chr'))
viburni_orthologs <- read.table(within_viburni_file, header = T)

viburni_orthologs$gene1 <- sapply(strsplit(viburni_orthologs$gene1, '.t'), function(x){ x[1] })
viburni_orthologs$gene2 <- sapply(strsplit(viburni_orthologs$gene2, '.t'), function(x){ x[1] })

# filtering all self-hits of alternative transcripts
viburni_orthologs <- viburni_orthologs[viburni_orthologs$gene1 != viburni_orthologs$gene2, ]

row.names(B_genes_tab) <- B_genes_tab$gene

viburni_orthologs$ch1 <- B_genes_tab[viburni_orthologs$gene1, 'chr']
viburni_orthologs$ch2 <- B_genes_tab[viburni_orthologs$gene2, 'chr']

viburni_orthologs[is.na(viburni_orthologs$ch1), 'ch1'] <- 'A'
viburni_orthologs[is.na(viburni_orthologs$ch2), 'ch2'] <- 'A'

B_viburni <- viburni_orthologs[!c(viburni_orthologs$ch1 == 'A' & viburni_orthologs$ch2 == 'A'), ]
# removing alternative transcripts
B_viburni <- B_viburni[!duplicated(paste0(B_viburni$gene1, B_viburni$gene2)), ]

B1s <- c(B_viburni$ch1 == 'B1' | B_viburni$ch2 == 'B1')
B2s <- c(B_viburni$ch1 == 'B2' | B_viburni$ch2 == 'B2')
B3s <- c(B_viburni$ch1 == 'B3' | B_viburni$ch2 == 'B3')

png('ortholog_identity_B_to_core_viburni_genome.png')

hist(B_viburni$identity, breaks = 30, main = 'nucleotide identity of within viburni B orthologs', xlab = 'nt identity')
hist(B_viburni$identity[B3s], breaks = 30, add = T, col = 'green')
hist(B_viburni$identity[B2s], breaks = 15, add = T, col = 'blue')
hist(B_viburni$identity[B1s], breaks = 15, add = T, col = 'red')
legend('topleft', pch = 20, col = c('red', 'blue', 'green'), c('B1', 'B2', 'B3'))

dev.off()

solenopsis_orthologs <- read.table(solenopsis_file)
longispinus_orthologs <- read.table(longispinus_file)

colnames(solenopsis_orthologs) <- colnames(B_viburni)[1:7]
colnames(longispinus_orthologs) <- colnames(B_viburni)[1:7]

solenopsis_orthologs <- solenopsis_orthologs[!grepl('Pvib', solenopsis_orthologs$gene1), ]
longispinus_orthologs <- longispinus_orthologs[!grepl('Pvib', longispinus_orthologs$gene2), ]

png('ortholog_identity_B_to_other_species.png')

hist(B_viburni$identity, breaks = 30, main = 'nucleotide identity of B orthologs to three species', xlab = 'nt identity')
hist(longispinus_orthologs$identity, add = T, col = 'cyan', breaks = 30)
hist(solenopsis_orthologs$identity, add = T, col = 'magenta', breaks = 30)

legend('topleft', col = c('grey', 'cyan', 'magenta'), c('viburni', 'longispinus', 'solenopsis'), pch = 20)

dev.off()
```

### Exploration of the chromosomal origin

Solenopsis has a really sweet chr. lvl assembly, which allows to search for colinearity between B scaffolds and individual chromosomes

#### nt proxy

This is just a quick peek using nt reciprocal blast. I saved content of `solenopsis_orthologs$gene1` vector in a file `solenopsis/B_orthologs.tsv` and grepped it out of solenopsis annotation

```
grep "gene" solenopsis/annotation.gff3 | grep -f solenopsis/B_orthologs.tsv
```

And it turns out that out of 11 orthologs, 2 come from chr2, 2 from chr3, 3 from chr5, 4 from chr4.

This would be better done using colinear blocks and inidividual genes might be a bit missleading when it comes to chromosomal origin (let's face it, some of them jump a bit around). Yet, given this peek, I would not expect a clear single-chromosome origin of Bs.
