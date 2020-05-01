
# Scott Barlow's mini project

Looking for differentially expressed transcripts between B+ and B- males using the RNAseq data/

## 1. Transcript expression estimation with kallisto and sleuth

	# create new environment: conda create -n afilia_trinity
	# prepare sample file:

sample	condition
04M_1	B+
04M_2	B+
04M_3	B+
13M_1	B+
13M_2	B+
13M_3	B+
13M_4	B+
15M_1	B-
15M_2	B-
15M_3	B-
21M_1	B-
21M_2	B-
21M_3	B-
21M_4	B-

	/ceph/users/afilia/.conda/envs/afilia_trinity/bin/align_and_estimate_abundance.pl --transcripts viburni.trinity.fasta --seqType fq --samples_file trinity_sample_file.txt --est_method kallisto --SS_lib_type RF --thread_count 32 --trinity_mode --prep_reference --kallisto_add_opts "-b 100"

## 2. Pass on results to sleuth 

My version of the analysis here: scott_miniproject_andres.R
