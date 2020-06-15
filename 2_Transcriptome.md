
# Transcriptome assembly and annotation

	# working directory	
	/data/ross/mealybugs/analyses/B_viburni_andres/2_short_read_DNA_seq/0_reads
	# raw reads
	/data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812

## 1. Raw RNAseq reads

	# Novaseq (/data/ross/mealybugs/analyses/B_viburni_andres/2_short_read_DNA_seq/0_reads)
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/13F_1/190805_A00291_0195_AHKLMLDMXX_2_11791RL0015L01_1.fastq.gz 13F_1_1.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/13F_1/190805_A00291_0195_AHKLMLDMXX_2_11791RL0015L01_2.fastq.gz 13F_1_2.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/13F_2/190805_A00291_0195_AHKLMLDMXX_2_11791RL0016L01_1.fastq.gz 13F_2_1.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/13F_2/190805_A00291_0195_AHKLMLDMXX_2_11791RL0016L01_2.fastq.gz 13F_2_2.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/13F_3/190805_A00291_0195_AHKLMLDMXX_2_11791RL0017L01_1.fastq.gz 13F_3_1.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/13F_3/190805_A00291_0195_AHKLMLDMXX_2_11791RL0017L01_2.fastq.gz 13F_3_2.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/13M_1/190805_A00291_0195_AHKLMLDMXX_2_11791RL0018L01_1.fastq.gz 13M_1_1.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/13M_1/190805_A00291_0195_AHKLMLDMXX_2_11791RL0018L01_2.fastq.gz 13M_1_2.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/13M_2/190805_A00291_0195_AHKLMLDMXX_2_11791RL0019L01_1.fastq.gz 13M_2_1.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/13M_2/190805_A00291_0195_AHKLMLDMXX_2_11791RL0019L01_2.fastq.gz 13M_2_2.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/13M_3/190805_A00291_0195_AHKLMLDMXX_2_11791RL0020L01_1.fastq.gz 13M_3_1.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/13M_3/190805_A00291_0195_AHKLMLDMXX_2_11791RL0020L01_2.fastq.gz 13M_3_2.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/13M_4/190805_A00291_0195_AHKLMLDMXX_2_11791RL0021L01_1.fastq.gz 13M_4_1.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/13M_4/190805_A00291_0195_AHKLMLDMXX_2_11791RL0021L01_2.fastq.gz 13M_4_2.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/15F_1/190805_A00291_0195_AHKLMLDMXX_2_11791RL0022L01_1.fastq.gz 15F_1_1.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/15F_1/190805_A00291_0195_AHKLMLDMXX_2_11791RL0022L01_2.fastq.gz 15F_1_2.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/15F_2/190805_A00291_0195_AHKLMLDMXX_2_11791RL0023L01_1.fastq.gz 15F_2_1.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/15F_2/190805_A00291_0195_AHKLMLDMXX_2_11791RL0023L01_2.fastq.gz 15F_2_2.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/15F_3/190805_A00291_0195_AHKLMLDMXX_2_11791RL0024L01_1.fastq.gz 15F_3_1.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/15F_3/190805_A00291_0195_AHKLMLDMXX_2_11791RL0024L01_2.fastq.gz 15F_3_2.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/15M_1/190805_A00291_0195_AHKLMLDMXX_2_11791RL0025L01_1.fastq.gz 15M_1_1.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/15M_1/190805_A00291_0195_AHKLMLDMXX_2_11791RL0025L01_2.fastq.gz 15M_1_2.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/15M_2/190805_A00291_0195_AHKLMLDMXX_2_11791RL0026L01_1.fastq.gz 15M_2_1.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/15M_2/190805_A00291_0195_AHKLMLDMXX_2_11791RL0026L01_2.fastq.gz 15M_2_2.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/15M_3/190805_A00291_0195_AHKLMLDMXX_2_11791RL0027L01_1.fastq.gz 15M_3_1.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/15M_3/190805_A00291_0195_AHKLMLDMXX_2_11791RL0027L01_2.fastq.gz 15M_3_2.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/21F_1/190805_A00291_0195_AHKLMLDMXX_2_11791RL0028L01_1.fastq.gz 21F_1_1.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/21F_1/190805_A00291_0195_AHKLMLDMXX_2_11791RL0028L01_2.fastq.gz 21F_1_2.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/21F_2/190805_A00291_0195_AHKLMLDMXX_2_11791RL0029L01_1.fastq.gz 21F_2_1.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/21F_2/190805_A00291_0195_AHKLMLDMXX_2_11791RL0029L01_2.fastq.gz 21F_2_2.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/21F_3/190805_A00291_0195_AHKLMLDMXX_2_11791RL0030L01_1.fastq.gz 21F_3_1.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/21F_3/190805_A00291_0195_AHKLMLDMXX_2_11791RL0030L01_2.fastq.gz 21F_3_2.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/21M_1/190805_A00291_0195_AHKLMLDMXX_2_11791RL0031L01_1.fastq.gz 21M_1_1.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/21M_1/190805_A00291_0195_AHKLMLDMXX_2_11791RL0031L01_2.fastq.gz 21M_1_2.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/21M_2/190805_A00291_0195_AHKLMLDMXX_2_11791RL0032L01_1.fastq.gz 21M_2_1.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/21M_2/190805_A00291_0195_AHKLMLDMXX_2_11791RL0032L01_2.fastq.gz 21M_2_2.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/21M_3/190805_A00291_0195_AHKLMLDMXX_2_11791RL0033L01_1.fastq.gz 21M_3_1.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/21M_3/190805_A00291_0195_AHKLMLDMXX_2_11791RL0033L01_2.fastq.gz 21M_3_2.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/21M_4/190805_A00291_0195_AHKLMLDMXX_2_11791RL0034L01_1.fastq.gz 21M_4_1.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/21M_4/190805_A00291_0195_AHKLMLDMXX_2_11791RL0034L01_2.fastq.gz 21M_4_2.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/4F_1/190805_A00291_0195_AHKLMLDMXX_2_11791RL0009L01_1.fastq.gz 04F_1_1.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/4F_1/190805_A00291_0195_AHKLMLDMXX_2_11791RL0009L01_2.fastq.gz 04F_1_2.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/4F_2/190805_A00291_0195_AHKLMLDMXX_2_11791RL0010L01_1.fastq.gz 04F_2_1.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/4F_2/190805_A00291_0195_AHKLMLDMXX_2_11791RL0010L01_2.fastq.gz 04F_2_2.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/4F_3/190805_A00291_0195_AHKLMLDMXX_2_11791RL0011L01_1.fastq.gz 04F_3_1.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/4F_3/190805_A00291_0195_AHKLMLDMXX_2_11791RL0011L01_2.fastq.gz 04F_3_2.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/4M_1/190805_A00291_0195_AHKLMLDMXX_2_11791RL0012L01_1.fastq.gz 04M_1_1.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/4M_1/190805_A00291_0195_AHKLMLDMXX_2_11791RL0012L01_2.fastq.gz 04M_1_2.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/4M_2/190805_A00291_0195_AHKLMLDMXX_2_11791RL0013L01_1.fastq.gz 04M_2_1.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/4M_2/190805_A00291_0195_AHKLMLDMXX_2_11791RL0013L01_2.fastq.gz 04M_2_2.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/4M_3/190805_A00291_0195_AHKLMLDMXX_2_11791RL0014L01_1.fastq.gz 04M_3_1.fastq.gz
	ln -s /data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/4M_3/190805_A00291_0195_AHKLMLDMXX_2_11791RL0014L01_2.fastq.gz 04M_3_2.fastq.gz

## 2. Trim reads

	# fastp version 0.20.0 (conda env afilia)
	fastp -i 13F_1_1.fastq.gz -I 13F_1_2.fastq.gz -o /scratch/afilia/13F_1.trimmed_1.fastq.gz -O /scratch/afilia/13F_1.trimmed_2.fastq.gz --cut_by_quality5 --cut_by_quality3 --cut_window_size 4 --cut_mean_quality 20 --detect_adapter_for_pe --trim_poly_g
	fastp -i 13F_2_1.fastq.gz -I 13F_2_2.fastq.gz -o /scratch/afilia/13F_2.trimmed_1.fastq.gz -O /scratch/afilia/13F_2.trimmed_2.fastq.gz --cut_by_quality5 --cut_by_quality3 --cut_window_size 4 --cut_mean_quality 20 --detect_adapter_for_pe 	--trim_poly_g
	fastp -i 13F_3_1.fastq.gz -I 13F_3_2.fastq.gz -o /scratch/afilia/13F_3.trimmed_1.fastq.gz -O /scratch/afilia/13F_3.trimmed_2.fastq.gz --cut_by_quality5 --cut_by_quality3 --cut_window_size 4 --cut_mean_quality 20 --detect_adapter_for_pe 	--trim_poly_g
	fastp -i 13M_1_1.fastq.gz -I 13M_1_2.fastq.gz -o /scratch/afilia/13M_1.trimmed_1.fastq.gz -O /scratch/afilia/13M_1.trimmed_2.fastq.gz --cut_by_quality5 --cut_by_quality3 --cut_window_size 4 --cut_mean_quality 20 --detect_adapter_for_pe 	--trim_poly_g
	fastp -i 13M_2_1.fastq.gz -I 13M_2_2.fastq.gz -o /scratch/afilia/13M_2.trimmed_1.fastq.gz -O /scratch/afilia/13M_2.trimmed_2.fastq.gz --cut_by_quality5 --cut_by_quality3 --cut_window_size 4 --cut_mean_quality 20 --detect_adapter_for_pe 	--trim_poly_g
	fastp -i 13M_3_1.fastq.gz -I 13M_3_2.fastq.gz -o /scratch/afilia/13M_3.trimmed_1.fastq.gz -O /scratch/afilia/13M_3.trimmed_2.fastq.gz --cut_by_quality5 --cut_by_quality3 --cut_window_size 4 --cut_mean_quality 20 --detect_adapter_for_pe 	--trim_poly_g
	fastp -i 13M_4_1.fastq.gz -I 13M_4_2.fastq.gz -o /scratch/afilia/13M_4.trimmed_1.fastq.gz -O /scratch/afilia/13M_4.trimmed_2.fastq.gz --cut_by_quality5 --cut_by_quality3 --cut_window_size 4 --cut_mean_quality 20 --detect_adapter_for_pe 	--trim_poly_g
	fastp -i 15F_1_1.fastq.gz -I 15F_1_2.fastq.gz -o /scratch/afilia/15F_1.trimmed_1.fastq.gz -O /scratch/afilia/15F_1.trimmed_2.fastq.gz --cut_by_quality5 --cut_by_quality3 --cut_window_size 4 --cut_mean_quality 20 --detect_adapter_for_pe 	--trim_poly_g
	fastp -i 15F_2_1.fastq.gz -I 15F_2_2.fastq.gz -o /scratch/afilia/15F_2.trimmed_1.fastq.gz -O /scratch/afilia/15F_2.trimmed_2.fastq.gz --cut_by_quality5 --cut_by_quality3 --cut_window_size 4 --cut_mean_quality 20 --detect_adapter_for_pe 	--trim_poly_g
	fastp -i 15F_3_1.fastq.gz -I 15F_3_2.fastq.gz -o /scratch/afilia/15F_3.trimmed_1.fastq.gz -O /scratch/afilia/15F_3.trimmed_2.fastq.gz --cut_by_quality5 --cut_by_quality3 --cut_window_size 4 --cut_mean_quality 20 --detect_adapter_for_pe 	--trim_poly_g
	fastp -i 15M_1_1.fastq.gz -I 15M_1_2.fastq.gz -o /scratch/afilia/15M_1.trimmed_1.fastq.gz -O /scratch/afilia/15M_1.trimmed_2.fastq.gz --cut_by_quality5 --cut_by_quality3 --cut_window_size 4 --cut_mean_quality 20 --detect_adapter_for_pe 	--trim_poly_g
	fastp -i 15M_2_1.fastq.gz -I 15M_2_2.fastq.gz -o /scratch/afilia/15M_2.trimmed_1.fastq.gz -O /scratch/afilia/15M_2.trimmed_2.fastq.gz --cut_by_quality5 --cut_by_quality3 --cut_window_size 4 --cut_mean_quality 20 --detect_adapter_for_pe 	--trim_poly_g
	fastp -i 15M_3_1.fastq.gz -I 15M_3_2.fastq.gz -o /scratch/afilia/15M_3.trimmed_1.fastq.gz -O /scratch/afilia/15M_3.trimmed_2.fastq.gz --cut_by_quality5 --cut_by_quality3 --cut_window_size 4 --cut_mean_quality 20 --detect_adapter_for_pe 	--trim_poly_g
	fastp -i 21F_1_1.fastq.gz -I 21F_1_2.fastq.gz -o /scratch/afilia/21F_1.trimmed_1.fastq.gz -O /scratch/afilia/21F_1.trimmed_2.fastq.gz --cut_by_quality5 --cut_by_quality3 --cut_window_size 4 --cut_mean_quality 20 --detect_adapter_for_pe 	--trim_poly_g
	fastp -i 21F_2_1.fastq.gz -I 21F_2_2.fastq.gz -o /scratch/afilia/21F_2.trimmed_1.fastq.gz -O /scratch/afilia/21F_2.trimmed_2.fastq.gz --cut_by_quality5 --cut_by_quality3 --cut_window_size 4 --cut_mean_quality 20 --detect_adapter_for_pe 	--trim_poly_g
	fastp -i 21F_3_1.fastq.gz -I 21F_3_2.fastq.gz -o /scratch/afilia/21F_3.trimmed_1.fastq.gz -O /scratch/afilia/21F_3.trimmed_2.fastq.gz --cut_by_quality5 --cut_by_quality3 --cut_window_size 4 --cut_mean_quality 20 --detect_adapter_for_pe 	--trim_poly_g
	fastp -i 21M_1_1.fastq.gz -I 21M_1_2.fastq.gz -o /scratch/afilia/21M_1.trimmed_1.fastq.gz -O /scratch/afilia/21M_1.trimmed_2.fastq.gz --cut_by_quality5 --cut_by_quality3 --cut_window_size 4 --cut_mean_quality 20 --detect_adapter_for_pe 	--trim_poly_g
	fastp -i 21M_2_1.fastq.gz -I 21M_2_2.fastq.gz -o /scratch/afilia/21M_2.trimmed_1.fastq.gz -O /scratch/afilia/21M_2.trimmed_2.fastq.gz --cut_by_quality5 --cut_by_quality3 --cut_window_size 4 --cut_mean_quality 20 --detect_adapter_for_pe 	--trim_poly_g
	fastp -i 21M_3_1.fastq.gz -I 21M_3_2.fastq.gz -o /scratch/afilia/21M_3.trimmed_1.fastq.gz -O /scratch/afilia/21M_3.trimmed_2.fastq.gz --cut_by_quality5 --cut_by_quality3 --cut_window_size 4 --cut_mean_quality 20 --detect_adapter_for_pe 	--trim_poly_g
	fastp -i 21M_4_1.fastq.gz -I 21M_4_2.fastq.gz -o /scratch/afilia/21M_4.trimmed_1.fastq.gz -O /scratch/afilia/21M_4.trimmed_2.fastq.gz --cut_by_quality5 --cut_by_quality3 --cut_window_size 4 --cut_mean_quality 20 --detect_adapter_for_pe 	--trim_poly_g
	fastp -i 04F_1_1.fastq.gz -I 04F_1_2.fastq.gz -o /scratch/afilia/04F_1.trimmed_1.fastq.gz -O /scratch/afilia/04F_1.trimmed_2.fastq.gz --cut_by_quality5 --cut_by_quality3 --cut_window_size 4 --cut_mean_quality 20 --detect_adapter_for_pe 	--trim_poly_g
	fastp -i 04F_2_1.fastq.gz -I 04F_2_2.fastq.gz -o /scratch/afilia/04F_2.trimmed_1.fastq.gz -O /scratch/afilia/04F_2.trimmed_2.fastq.gz --cut_by_quality5 --cut_by_quality3 --cut_window_size 4 --cut_mean_quality 20 --detect_adapter_for_pe 	--trim_poly_g
	fastp -i 04F_3_1.fastq.gz -I 04F_3_2.fastq.gz -o /scratch/afilia/04F_3.trimmed_1.fastq.gz -O /scratch/afilia/04F_3.trimmed_2.fastq.gz --cut_by_quality5 --cut_by_quality3 --cut_window_size 4 --cut_mean_quality 20 --detect_adapter_for_pe 	--trim_poly_g
	fastp -i 04M_1_1.fastq.gz -I 04M_1_2.fastq.gz -o /scratch/afilia/04M_1.trimmed_1.fastq.gz -O /scratch/afilia/04M_1.trimmed_2.fastq.gz --cut_by_quality5 --cut_by_quality3 --cut_window_size 4 --cut_mean_quality 20 --detect_adapter_for_pe 	--trim_poly_g
	fastp -i 04M_2_1.fastq.gz -I 04M_2_2.fastq.gz -o /scratch/afilia/04M_2.trimmed_1.fastq.gz -O /scratch/afilia/04M_2.trimmed_2.fastq.gz --cut_by_quality5 --cut_by_quality3 --cut_window_size 4 --cut_mean_quality 20 --detect_adapter_for_pe 	--trim_poly_g
	fastp -i 04M_3_1.fastq.gz -I 04M_3_2.fastq.gz -o /scratch/afilia/04M_3.trimmed_1.fastq.gz -O /scratch/afilia/04M_3.trimmed_2.fastq.gz --cut_by_quality5 --cut_by_quality3 --cut_window_size 4 --cut_mean_quality 20 --detect_adapter_for_pe 	--trim_poly_g

## 3. Trinity

	# Trinity-v2.8.5 - strand-specific
	Trinity --seqType fq --left 04F_1.trimmed_1.fastq.gz,04F_2.trimmed_1.fastq.gz,04F_3.trimmed_1.fastq.gz,04M_1.trimmed_1.fastq.gz,04M_2.trimmed_1.fastq.gz,04M_3.trimmed_1.fastq.gz,13F_1.trimmed_1.fastq.gz,13F_2.trimmed_1.fastq.gz,13F_3.trimmed_1.fastq.gz,13M_1.trimmed_1.fastq.gz,13M_2.trimmed_1.fastq.gz,13M_3.trimmed_1.fastq.gz,13M_4.trimmed_1.fastq.gz,15F_1.trimmed_1.fastq.gz,15F_2.trimmed_1.fastq.gz,15F_3.trimmed_1.fastq.gz,15M_1.trimmed_1.fastq.gz,15M_2.trimmed_1.fastq.gz,15M_3.trimmed_1.fastq.gz,21F_1.trimmed_1.fastq.gz,21F_2.trimmed_1.fastq.gz,21F_3.trimmed_1.fastq.gz,21M_1.trimmed_1.fastq.gz,21M_2.trimmed_1.fastq.gz,21M_3.trimmed_1.fastq.gz,21M_4.trimmed_1.fastq.gz --right 04F_1.trimmed_2.fastq.gz,04F_2.trimmed_2.fastq.gz,04F_3.trimmed_2.fastq.gz,04M_1.trimmed_2.fastq.gz,04M_2.trimmed_2.fastq.gz,04M_3.trimmed_2.fastq.gz,13F_1.trimmed_2.fastq.gz,13F_2.trimmed_2.fastq.gz,13F_3.trimmed_2.fastq.gz,13M_1.trimmed_2.fastq.gz,13M_2.trimmed_2.fastq.gz,13M_3.trimmed_2.fastq.gz,13M_4.trimmed_2.fastq.gz,15F_1.trimmed_2.fastq.gz,15F_2.trimmed_2.fastq.gz,15F_3.trimmed_2.fastq.gz,15M_1.trimmed_2.fastq.gz,15M_2.trimmed_2.fastq.gz,15M_3.trimmed_2.fastq.gz,21F_1.trimmed_2.fastq.gz,21F_2.trimmed_2.fastq.gz,21F_3.trimmed_2.fastq.gz,21M_1.trimmed_2.fastq.gz,21M_2.trimmed_2.fastq.gz,21M_3.trimmed_2.fastq.gz,21M_4.trimmed_2.fastq.gz --SS_lib_type RF --max_memory 100G --CPU 32 --full_cleanup

## 4. Annotate using the Trinonate pipeline

	conda install -c bioconda trinotate (3.2.0) # link to eggnoc doesn't work, had to update it manually using the 3.2.1 version (not in conda)
	conda install -c anaconda sqlite (3.31.1)
	# don't do signalP, tmhhm, rnammer for now

	# build sqlite database and prepare sequences
	Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate
	makeblastdb -in uniprot_sprot.pep -dbtype prot
	gunzip Pfam-A.hmm.gz
	hmmpress Pfam-A.hmm

	# use Transdecoder to predict ORFs and protein coding regions
	TransDecoder.LongOrfs -t ../viburni.trinity.fasta
	TransDecoder.Predict -t ../viburni.trinity.fasta (without homology options)

	# homology searches
	blastx -query ../viburni.trinity.fasta -db uniprot_sprot.pep -num_threads 16 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastx.outfmt6
	blastp -query viburni.trinity.fasta.transdecoder.pep -db uniprot_sprot.pep -num_threads 16 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastp.outfmt6
	hmmscan --cpu 16 --domtblout TrinotatePFAM.out Pfam-A.hmm viburni.trinity.fasta.transdecoder.pep > pfam.log

	# Load transcripts and coding regions
	Trinotate Trinotate.sqlite init --gene_trans_map ../viburni.trinity.fasta.gene_trans_map --transcript_fasta ../viburni.trinity.fasta --transdecoder_pep viburni.trinity.fasta.transdecoder.pep

	# Loading BLAST homologies
	Trinotate Trinotate.sqlite LOAD_swissprot_blastp blastp.outfmt6
	Trinotate Trinotate.sqlite LOAD_swissprot_blastx blastx.outfmt6
	Trinotate Trinotate.sqlite LOAD_pfam TrinotatePFAM.out

	# Generate report and extract GO terms
	Trinotate Trinotate.sqlite report -E 1e-3 > trinotate_annotation_report.xls
	extract_GO_assignments_from_Trinotate_xls.pl --Trinotate_xls trinotate_annotation_report.xls --trans --include_ancestral_terms > go_annotations_viburni_with_ancestral.txt

## 5. Explore the data

First, the distribution of samples looks okay -- no clear outlier

![](misc/tpm_distribution.jpeg)

We know from Scott's miniproject that the males cluster nicely according to B+/B-. Let's check that the samples cluster primarily by sex.

	/ceph/users/afilia/.conda/envs/afilia_trinity/bin/align_and_estimate_abundance.pl --transcripts ../1_trinity/viburni.trinity.fasta --seqType fq --samples_file trinity_sample_file.txt --est_method kallisto --SS_lib_type RF --thread_count 16 --trinity_mode --prep_reference --kallisto_add_opts "-b 100"
	# import tpm values to R, explore correlation
	tpm.cor.plot <- cor(tpm.merge[c(-1)], method="spearman", use = "complete.obs")
	corrplot(tpm.cor.plot, method = "number",number.cex = .7, cl.lim = c(0.3, 1),order = "hclust", addrect = 2)
	heatmap(tpm.cor.plot, scale = "row")

![](misc/rnaseq_heatmap.jpeg)

They cluster as expected (sex -> B status -> line). We see a clear clustering by sex in the PCA of estimated counts running sleuth, but no clear clustering by B status (note that the default threshold to filter out transcripts in sleuth has been lowered to at least 5 counts in 20% of samples)

![](misc/pca_sleuth.jpeg)

We see a similar thing when we look at the density plot of estimated counts.

![](misc/density_plot_kallisto.jpeg)

Let's inspect SPM distribition, which looks less extreme than in *P. citri* adults.

![](misc/spm_transcripts.jpeg)

Sleuth has estimated differential expression for 41,372 annotated transcripts (out of ca. 150,000 transcripts that passed the filter). Let's focus on annotated transcripts only for the time being. Estimating expression FC as log2((female TPM + 0.01)/(male TPM + 0.01)), the data looks like this:

![](misc/volcano_plot_sex.jpeg)

where NB are non-biased transcripts (q>0.05 and/or abs(FC) <= 1.5), FB and MB are sex-biased transcripts (q<0.05 and abs(FC) > 1.5) and FS and MS are sex-specific transcripts (sex-biased transcripts with FC > 2.5 and TPM in the other sex < 1).

- NB 35172
- FB  2644
- FS   747
- MB  1735
- MS  1074
- NB 35172

GO enrichment analysis for sex-biased transcripts

	find_enrichment.py --pval 0.05 --method fdr_bh --obo go-basic.obo --outfile results_sex/female.transcripts.GO.basic.tsv input_sex/female.biased.txt input_sex/backgound.pop.txt go_annotations_viburni_with_ancestral_for_GOATTOOLS.txt
	find_enrichment.py --pval 0.05 --method fdr_bh --obo go-basic.obo --outfile results_sex/male.transcripts.GO.basic.tsv input_sex/male.biased.txt input_sex/backgound.pop.txt go_annotations_viburni_with_ancestral_for_GOATTOOLS.txt

	find_enrichment.py --pval 0.05 --method fdr_bh --obo goslim_generic.obo --outfile results_sex/female.transcripts.GO.slim.tsv input_sex/female.biased.txt input_sex/backgound.pop.txt go_annotations_viburni_with_ancestral_for_GOATTOOLS.txt
	find_enrichment.py --pval 0.05 --method fdr_bh --obo goslim_generic.obo --outfile results_sex/male.transcripts.GO.slim.tsv input_sex/male.biased.txt input_sex/backgound.pop.txt go_annotations_viburni_with_ancestral_for_GOATTOOLS.txt

We have a lot of enriched/purified GO terms

 - 2536 female.transcripts.GO.basic.tsv
 - 118 female.transcripts.GO.slim.tsv
 - 837 male.transcripts.GO.basic.tsv
 - 52 male.transcripts.GO.slim.tsv

Let's narrow the datasets and look at sex-specific transcripts

	find_enrichment.py --pval 0.05 --method fdr_bh --obo go-basic.obo --outfile results_sex/female.sp.transcripts.GO.basic.tsv input_sex/female.sp.txt input_sex/backgound.pop.txt go_annotations_viburni_with_ancestral_for_GOATTOOLS.txt
	find_enrichment.py --pval 0.05 --method fdr_bh --obo go-basic.obo --outfile results_sex/male.sp.transcripts.GO.basic.tsv input_sex/male.sp.txt input_sex/backgound.pop.txt go_annotations_viburni_with_ancestral_for_GOATTOOLS.txt

	find_enrichment.py --pval 0.05 --method fdr_bh --obo goslim_generic.obo --outfile results_sex/female.sp.transcripts.GO.slim.tsv input_sex/female.sp.txt input_sex/backgound.pop.txt go_annotations_viburni_with_ancestral_for_GOATTOOLS.txt
	find_enrichment.py --pval 0.05 --method fdr_bh --obo goslim_generic.obo --outfile results_sex/male.sp.transcripts.GO.slim.tsv input_sex/male.sp.txt input_sex/backgound.pop.txt go_annotations_viburni_with_ancestral_for_GOATTOOLS.txt

## 6. Differential gene expression between B+ and B- samples (controlling for sex)

Sleuth has picked up 32,792 transcripts (23.1% out of 141,922 that passed the filters, and 16.3% of 200,993 transcripts) that are differentially expressed between B- and B+ samples. Of these, 9,518 transcrips are annotated (out of 41,372).

Let's categorise the transcripts according to the following criteria:

	sleuth.tpm.annotated$info <- "Equal expression"
	sleuth.tpm.annotated$info <- ifelse(sleuth.tpm.annotated$qval < 0.05 & sleuth.tpm.annotated$FC > 1.5,"Overexpressed in B+",sleuth.tpm.annotated$info)
	sleuth.tpm.annotated$info <- ifelse(sleuth.tpm.annotated$qval < 0.05 & sleuth.tpm.annotated$FC > 1.5 & sleuth.tpm.annotated$Bminus_F < 0.1 & 	sleuth.tpm.annotated$Bminus_M < 0.1,"Unique in B",sleuth.tpm.annotated$info)
	sleuth.tpm.annotated$info <- ifelse(sleuth.tpm.annotated$qval < 0.05 & sleuth.tpm.annotated$FC > 1.5 & sleuth.tpm.annotated$Bminus_F < 0.1 & 	sleuth.tpm.annotated$Bminus_M < 0.1 & sleuth.tpm.annotated$Bplus_F < 0.1 & sleuth.tpm.annotated$Bplus_M >= 0.5,"Unique in B males",sleuth.tpm.annotated$info)
	sleuth.tpm.annotated$info <- ifelse(sleuth.tpm.annotated$qval < 0.05 & sleuth.tpm.annotated$FC < -1.5,"Overexpressed in 0B",sleuth.tpm.annotated$info)

- Equal expression 		36192
- Overexpressed in 0B  	2113
- Overexpressed in B+   1909
- Unique in B  			1142
- Unique in B males  	16

![](misc/volcano_plot_B.jpeg)

Let's examine the 16 transcripts unique to B males

| transcript_id            | qval       | Bminus_F | Bminus_M | Bplus_F | Bplus_M | sprot_Top_BLASTX_hit | sprot_Top_BLASTP_hit                                                                  | FC    | info              |
|--------------------------|------------|----------|----------|---------|---------|----------------------|---------------------------------------------------------------------------------------|-------|-------------------|
| TRINITY_DN5504_c0_g1_i1  | 0.00517178 | 0        | 0.03     | 0.06    | 3.76    | TC3A_CAEEL           | Transposable element Tc3 transposase                                                  | 6.263 | Unique in B males |
| TRINITY_DN18538_c0_g3_i2 | 0.00064676 | 0        | 0.08     | 0       | 7.26    | CAD99_DROME          | Cadherin 99C                                                                          | 6.186 | Unique in B males |
| TRINITY_DN35228_c2_g1_i1 | 0.00077706 | 0        | 0.02     | 0       | 2.49    | PSH_DROME            | Serine protease persephone / Hayan                                                    | 5.972 | Unique in B males |
| TRINITY_DN16637_c2_g1_i1 | 0.00012846 | 0.01     | 0.04     | 0.03    | 3.78    | POLX_TOBAC           | Retrovirus-related Pol polyprotein from transposon TNT                                | 5.774 | Unique in B males |
| TRINITY_DN59503_c0_g2_i2 | 0.00073557 | 0.09     | 0.03     | 0.07    | 5.54    | SVOP_XENLA           | Synaptic vesicle 2-related protein                                                    | 5.33  | Unique in B males |
| TRINITY_DN446_c0_g1_i3   | 0.00803567 | 0.01     | 0.01     | 0.04    | 1.44    | PCBP3                | Poly(rC)-binding protein 3 (RNA binding proteins, translational activation/silencing) | 5.229 | Unique in B males |
| TRINITY_DN315_c0_g1_i2   | 0.00010257 | 0        | 0        | 0.07    | 0.64    | C6A13_DROME          | Probable cytochrome P450 6a13                                                         | 5.19  | Unique in B males |
| TRINITY_DN95495_c0_g1_i1 | 0.00029193 | 0        | 0.02     | 0.02    | 1.26    | ZIG8_CAEEL           | Zwei Ig domain protein zig-8                                                          | 5.022 | Unique in B males |
| TRINITY_DN6798_c0_g1_i1  | 0.00431661 | 0        | 0.01     | 0.09    | 0.85    | HIBCH_HUMAN          | Hydroxyisobutyryl-CoA hydrolase, mitochondrial                                        | 5     | Unique in B males |
| TRINITY_DN1761_c0_g1_i13 | 0.00589816 | 0        | 0.02     | 0.09    | 0.53    | PMS2_CHICK           | Mismatch repair endonuclease PMS2 (cell cycle/meiosis)                                | 4     | Unique in B males |
| TRINITY_DN2563_c1_g1_i8  | 0.00083441 | 0.04     | 0        | 0.04    | 0.78    | CP4G1_DROME          | Cytochrome P450 4g1                                                                   | 3.807 | Unique in B males |
| TRINITY_DN53103_c0_g1_i8 | 0.00053123 | 0        | 0.07     | 0.04    | 1.06    | FIG4_MOUSE           | Polyphosphoinositide phosphatase                                                      | 3.637 | Unique in B males |
| TRINITY_DN93_c0_g1_i2    | 0.00472272 | 0.01     | 0.06     | 0.01    | 1.09    | NA                   | Myb/SANT-like DNA-binding domain                                                      | 3.637 | Unique in B males |
| TRINITY_DN12933_c0_g3_i3 | 0.00943168 | 0        | 0.08     | 0       | 0.97    | LRFN4_MOUSE          | Leucine-rich repeat and fibronectin type-III domain-containing protein                | 3.307 | Unique in B males |
| TRINITY_DN6614_c0_g1_i14 | 0.02067184 | 0.01     | 0.04     | 0.03    | 0.61    | LAS1L_MOUSE          | Ribosomal biogenesis protein LAS1L                                                    | 3.237 | Unique in B males |
| TRINITY_DN28663_c0_g1_i1 | 0.03459472 | 0.03     | 0.08     | 0.04    | 0.59    | POL_SIVG             | Gag-Pol polyprotein                                                                   | 2.322 | Unique in B males |