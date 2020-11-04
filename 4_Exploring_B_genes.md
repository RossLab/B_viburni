
# Exploring B genes

	# working directory	
	/data/ross/mealybugs/analyses/B_viburni_2020/5_B_genes
	qlogin -pe smp64 32 -N bwa -l h=bigwig
    /ceph/software/utilities/sge/qlogin -pe smp64 32 -N bamfilter

We now have: 1) a genome assembly with an annotation (from blast, diamond, interproscan -- see 1_Genome_assembly), 2) a list of candidate B scaffolds (see 3_Coverage_analysis.md) and 3) lists of differentially expressed genes between B-carrying and non-carrying males and females. We can combine all this data and see what we can learn about the gene content of our putative B sequences.

## 1. Review our master annotation

Let's prepare a master file with annotated genes. For now, we are going to collapse transcripts into genes and add what we know from the different annotation sources [(R script here)](https://github.com/RossLab/B_viburni/blob/master/R_scripts/Gene_annotation.R). 8,914 genes have BLAST annotations, 10,915 have diamond annotations and 12,524 have a function assigned by interproscan (ignoring for now GO annotations, etc). In total, 13,515 genes are annotated, which is 57% of the predicted genes.

## 2. Genes on the B scaffolds

The R script is [here](https://github.com/RossLab/B_viburni/blob/master/R_scripts/Exploring_AB_genes.R). Let's start by exploring which genes fall in B candidate regions.