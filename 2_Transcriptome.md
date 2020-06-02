
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

We know from Scott's miniproject that the males cluster nicely according to B+/B-. Let's check that the samples cluster primarily by sex.

	/ceph/users/afilia/.conda/envs/afilia_trinity/bin/align_and_estimate_abundance.pl --transcripts ../1_trinity/viburni.trinity.fasta --seqType fq --samples_file trinity_sample_file.txt --est_method kallisto --SS_lib_type RF --thread_count 16 --trinity_mode --prep_reference --kallisto_add_opts "-b 100"
	# import tpm values to R, explore correlation
	tpm.cor.plot <- cor(tpm.merge[c(-1)], method="spearman", use = "complete.obs")
	corrplot(tpm.cor.plot, method = "number",number.cex = .7, cl.lim = c(0.3, 1),order = "hclust", addrect = 2)
	heatmap(tpm.cor.plot, scale = "row")

![](misc/rnaseq_heatmap.jpeg)

They cluster as expected (sex -> B status -> line). We see a clear clustering by sex in the PCA of estimated counts running sleuth, but no clear clustering by B status (note that the default threshold to filter out transcripts in sleuth has been lowered to at least 5 counts in 20% of samples)

![](misc/pca_sleuth.jpeg)
