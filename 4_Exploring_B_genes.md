
# Exploring B genes

	# working directory	
	/data/ross/mealybugs/analyses/B_viburni_2020/5_B_genes
	qlogin -pe smp64 32 -N bwa -l h=bigwig
    /ceph/software/utilities/sge/qlogin -pe smp64 32 -N bamfilter

We now have: 1) a genome assembly with an annotation (from blast, diamond, interproscan -- see 1_Genome_assembly), 2) a list of candidate B scaffolds (see 3_Coverage_analysis.md) and 3) lists of differentially expressed genes between B-carrying and non-carrying males and females. We can combine all this data and see what we can learn about the gene content of our putative B sequences.

## 1. Review our master annotation

Let's prepare a master file with annotated genes. For now, we are going to collapse transcripts into genes and add what we know from the different annotation sources [(R script here)](https://github.com/RossLab/B_viburni/blob/master/R_scripts/Gene_annotation.R).