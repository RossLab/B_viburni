
# Coverage analysis

	# working directory	
	/data/ross/mealybugs/analyses/B_viburni_andres/2_short_read_DNA_seq/0_reads

## 1. Raw reads

	# Novaseq reads (/data/ross/mealybugs/analyses/B_viburni_andres/2_short_read_DNA_seq/0_reads)
	ln -s 180608_A00291_0042_BH3CC3DRXX_2_11372RL0001L01_1.fastq.gz -> /data/ross/mealybugs/raw/transfer.genomics.ed.ac.uk/11372_Ross_Laura/raw_data/all_reads/18_13_1B_550/180608_A00291_0042_BH3CC3DRXX_2_11372RL0001L01_1.fastq.gz
	ln -s 180608_A00291_0042_BH3CC3DRXX_2_11372RL0001L01_2.fastq.gz -> /data/ross/mealybugs/raw/transfer.genomics.ed.ac.uk/11372_Ross_Laura/raw_data/all_reads/18_13_1B_550/180608_A00291_0042_BH3CC3DRXX_2_11372RL0001L01_2.fastq.gz
	ln -s 180608_A00291_0042_BH3CC3DRXX_2_11372RL0002L01_1.fastq.gz -> /data/ross/mealybugs/raw/transfer.genomics.ed.ac.uk/11372_Ross_Laura/raw_data/all_reads/18_13_1B_350/180608_A00291_0042_BH3CC3DRXX_2_11372RL0002L01_1.fastq.gz
	ln -s 180608_A00291_0042_BH3CC3DRXX_2_11372RL0002L01_2.fastq.gz -> /data/ross/mealybugs/raw/transfer.genomics.ed.ac.uk/11372_Ross_Laura/raw_data/all_reads/18_13_1B_350/180608_A00291_0042_BH3CC3DRXX_2_11372RL0002L01_2.fastq.gz

18_21_0B/180608_A00291_0042_BH3CC3DRXX_2_11372RL0004L01_1.fastq.gz
18_21_0B/180608_A00291_0042_BH3CC3DRXX_2_11372RL0004L01_2.fastq.gz
18_23_0B/180608_A00291_0042_BH3CC3DRXX_2_11372RL0005L01_1.fastq.gz
18_23_0B/180608_A00291_0042_BH3CC3DRXX_2_11372RL0005L01_2.fastq.gz
18_4_2B/180608_A00291_0042_BH3CC3DRXX_2_11372RL0003L01_1.fastq.gz
18_4_2B/180608_A00291_0042_BH3CC3DRXX_2_11372RL0003L01_2.fastq.gz

## 2. Trim with fastp

	# fastp version 0.20.0 (conda env afilia)
	fastp -i 180608_A00291_0042_BH3CC3DRXX_2_11372RL0001L01_1.fastq.gz -I 180608_A00291_0042_BH3CC3DRXX_2_11372RL0001L01_2.fastq.gz -o /scratch/afilia/PV_18-13.Illumina.550.trimmed_1.fq.gz -O /scratch/afilia/PV_18-13.Illumina.550.trimmed_2.fq.gz --cut_by_quality5 --cut_by_quality3 --cut_window_size 4 --cut_mean_quality 20 --detect_adapter_for_pe --trim_poly_g
	fastp -i 180608_A00291_0042_BH3CC3DRXX_2_11372RL0002L01_1.fastq.gz -I 180608_A00291_0042_BH3CC3DRXX_2_11372RL0002L01_2.fastq.gz -o /scratch/afilia/PV_18-13.Illumina.350.trimmed_1.fq.gz -O /scratch/afilia/PV_18-13.Illumina.350.trimmed_2.fq.gz --cut_by_quality5 --cut_by_quality3 --cut_window_size 4 --cut_mean_quality 20 --detect_adapter_for_pe --trim_poly_g

Short read data			Nseq	Poor qual	Range	GC %	Dedup %	Avg length
PV_18-13.Illumina.350.trimmed_1.fq.gz	104477186.0	0.0	15-150	34.0	71.3	147.2
PV_18-13.Illumina.350.trimmed_2.fq.gz	104477186.0	0.0	15-150	34.0	71.7	147.2
PV_18-13.Illumina.550.trimmed_1.fq.gz	81156002.0	0.0	16-150	34.0	75.5	146.8
PV_18-13.Illumina.550.trimmed_2.fq.gz	81156002.0	0.0	15-150	34.0	77.3	146.7
