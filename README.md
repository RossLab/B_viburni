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
