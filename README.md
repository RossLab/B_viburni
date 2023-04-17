# B_viburni
 Characterising the B chromosome of *Pseudococcus viburni*.


### Regenerating paper figures and tables

Preparing directory for all figures

```
mkdir -p manuscript/figures_revision
```

---

For panels of
 - Figure 3: Assignment of B chromosomes

```
Rscript R_scripts/B_scaffod_assignment_plots.R
```

---

For
- Figure 6: B+ vs B- gene expression analysis
- Figure S4: Distribution of library sizes and normalised log-CPMs and multi-dimensional scaling plots for the full differential expression model (sex and B presence/absence as factors)

```
Rscript R_scripts/B_differential_expression_plots.R
```

---

For
 - [Table of annotated genes located on B scaffolds](output/genes.in.Bs.anno.csv) -
 - Figure S3: Overall expression of B+ and B- lines

```
Rscript R_scripts/Exploring_AB_genes.R
```
---

LINK Figure 1
- Figure 1: A scheme of B chromosome tranmission
- Figure 2: Staining of B chromosomes


### Files

To be revised:
 - `output/rbbh_AB_transcripts.csv`


Revised files:
 - B assignment:
    - `output/B_scaffold_assignment_comeplete_window_coverage_table.tsv`: complete table with coverage information per 10k window
    - `output/scaffolds.final.assignment.tsv`: agregated table of all three sources of evidence into one table
    - `output/scaffolds.final.assignment.table.csv`: simplified table with assignments (scf; len; asn) - Supplementary file 1


 - Gene origin & function:
   - `output/B_genes.diamond_taxa_overview.tsv` - diamond hits of B-genes against ncbi nr
   - `output/genes.in.Bs.anno.csv` - functional annotation of B-linked genes

Preprint surviving files:
 - `output/pviburni.gene.GO`: GO terms for all genes
 - `output/freeze.v0.genes.anno.complete.csv`: agregated functional annotation for each gene
 - `annotation/p.viburni.freeze.v0.braker.transcripts.to.genes.txt`: gene <-> transcript table
 - `R_scripts/RSEM_digi.counts.matrix`: expression matrix


Preprint legacy files (_DO NOT USE!_):
 - `output/scaffolds.preprint.assignment.csv`: previous assignments
