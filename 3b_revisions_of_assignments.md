### Revisions of assignments

This document shows how we re-done the coverage analysis, compared to the preprint assignments presented in [3_Coverage_analysis document](3_Coverage_analysis.md).

And finnal assignments use the original k-mer and assembly-based assginments, but we reduced the complexity of assignments by streamlining the rules:

All scaffolds well-supported by coverage differences between lines AND that have (k-mer support OR assembly support), are considered as reliably B-linked (`B`).

All scaffolds with a loose support by coverage differences between lines OR k-mer support OR assembly support are considered as "B candidates" (`Bc`). It is very likely that some of these will be false positives, and their interpretation should be treated with caution.

Also note that some of the B-linked scaffolds are expected to have homologous sequences in the core genome too.

### Coverage assignments

We filtered all non-perfectly-matching-reads when we generated the preprint assignments. However, during revisions we found out the divergences of the B+ and B- lineages to the reference genome vary and therefore such filteration created some unfortunate false positives (while most of the true positives were identified correctly).

Here is the analysis redone on just quality filtered mapping files and using median normalisation before calculating coverage ratios. We also deplyed the information about the B-copy number in the decision making process.

#### Extracting per-window coverage

Using `PV_18-??.initial.sorted.primary.only.bam` files (`/data/ross/mealybugs/analyses/B_viburni_2020/4_cov_analysis/cov` on cluster). I extracted-per 10k window (and per 1k window) coverages using `samtools depth` and I perl script. Like this (for all 4 strains)

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

I will load the tables in R and merge the four strains in a single table (`scf`, `window`, `cov_04`, ...). This would require `.bed` files that are NOT in this repository, however, the merged table is and the transformation was pretty streightforward (see the script)

```bash
Rscript R_scripts/B_scaffod_assignment_coverage_table.R
```

Now I have list with the coverage of all windows, all I need to do now is to merge them. They are constructed with the same reference, so there is a guarantee the order and number of windows is the same in all four files.

### Re-evaluating fractions of assigned kmers

So instead I will at least check for considency of individual assignments with the suboptimal kmer approach we already have implemented (`R_scripts/explore_assignments.R`).

All those plots seem quite alright. So, I will keep those.

### Final assignment

The final assignment table will be generated via (`R_scripts/B_scaffold_reassignments.R`) from the "preprint assignment" table and the table with updated coverage analysis generated above (and also included in this repository).

The two assignment approaches are as consistent as possible in logic, separating confidently B-linked and putative (candidate) B-linked scaffolds.

```
Rscript R_scripts/B_scaffold_reassignments.R
```
