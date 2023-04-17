
# Illumina coverage analysis
Current location as of December 2020

	# working directory
	/data/ross/mealybugs/analyses/B_viburni_2020/isabelle_old_2019/pseudococcus_viburni/
	# raw reads
	data/ross/mealybugs/analyses/B_viburni_2020/isabelle_old_2019/pseudococcus_viburni/0_reads/


# I. GENOME ASSEMBLY (copied from Quiver from 2018)

```
$ tmux attach-session -t genomics
```
## Made shortcut in 0_reads (26 June 2018)
```
cd 0_reads/

ln -s /data/ross/mealybugs/raw/transfer.genomics.ed.ac.uk/11372_Ross_Laura/raw_data/data_by_date/20180618/all_reads/18_13_1B_350/180608_A00291_0042_BH3CC3DRXX_2_11372RL0002L01_1.fastq.gz pviburni.1813.1B.350.r1.fastq.gz
ln -s /data/ross/mealybugs/raw/transfer.genomics.ed.ac.uk/11372_Ross_Laura/raw_data/data_by_date/20180618/all_reads/18_13_1B_350/180608_A00291_0042_BH3CC3DRXX_2_11372RL0002L01_2.fastq.gz pviburni.1813.1B.350.r2.fastq.gz

ln -s /data/ross/mealybugs/raw/transfer.genomics.ed.ac.uk/11372_Ross_Laura/raw_data/data_by_date/20180618/all_reads/18_23_0B/180608_A00291_0042_BH3CC3DRXX_2_11372RL0005L01_1.fastq.gz pviburni.1823.0B.350.r1.fastq.gz
ln -s /data/ross/mealybugs/raw/transfer.genomics.ed.ac.uk/11372_Ross_Laura/raw_data/data_by_date/20180618/all_reads/18_23_0B/180608_A00291_0042_BH3CC3DRXX_2_11372RL0005L01_2.fastq.gz pviburni.1823.0B.350.r2.fastq.gz

ln -s /data/ross/mealybugs/raw/transfer.genomics.ed.ac.uk/11372_Ross_Laura/raw_data/data_by_date/20180618/all_reads/18_13_1B_550/180608_A00291_0042_BH3CC3DRXX_2_11372RL0001L01_1.fastq.gz pviburni.1813.1B.550.r1.fastq.gz ln -s /data/ross/mealybugs/raw/transfer.genomics.ed.ac.uk/11372_Ross_Laura/raw_data/data_by_date/20180618/all_reads/18_13_1B_550/180608_A00291_0042_BH3CC3DRXX_2_11372RL0001L01_2.fastq.gz pviburni.1813.1B.550.r2.fastq.gz

ln -s /data/ross/mealybugs/raw/transfer.genomics.ed.ac.uk/11372_Ross_Laura/raw_data/data_by_date/20180618/all_reads/18_4_2B/180608_A00291_0042_BH3CC3DRXX_2_11372RL0003L01_1.fastq.gz pviburni.184.2B.350.r1.fastq.gz
ln -s /data/ross/mealybugs/raw/transfer.genomics.ed.ac.uk/11372_Ross_Laura/raw_data/data_by_date/20180618/all_reads/18_4_2B/180608_A00291_0042_BH3CC3DRXX_2_11372RL0003L01_2.fastq.gz pviburni.184.2B.350.r2.fastq.gz

ln -s /data/ross/mealybugs/raw/transfer.genomics.ed.ac.uk/11372_Ross_Laura/raw_data/data_by_date/20180618/all_reads/18_21_0B/180608_A00291_0042_BH3CC3DRXX_2_11372RL0004L01_1.fastq.gz pviburni.1821.0B.350.r1.fastq.gz
ln -s /data/ross/mealybugs/raw/transfer.genomics.ed.ac.uk/11372_Ross_Laura/raw_data/data_by_date/20180618/all_reads/18_21_0B/180608_A00291_0042_BH3CC3DRXX_2_11372RL0004L01_2.fastq.gz pviburni.1821.0B.350.r2.fastq.gz

```
OK

```
conda activate main_env
```
OK


## counting reads and bases
sub file: readbasecount.sub
```
/ceph/software/readbasecount/readbasecount.py 0_reads/*.fastq.gz > 0_reads/read_base_count.txt
```
submission command line
```
parallel -j1 'qsub -cwd -N readbase -V -pe **smp64** 16 -b yes {}' :::: 0_reads/readbasecount.sub
```

OUTPUT:
```
file=0_reads/pviburni.1813.1B.350.r1.fastq.gz count=107,793,344 bases=16,169,001,600

```

## fastqc
subfile: fastqc.sub
```
fastqc --nogroup --outdir 1_fastqc/ 0_reads/pviburni.1813.1B.350.r1.fastq.gz
fastqc --nogroup --outdir 1_fastqc/ 0_reads/pviburni.1813.1B.550.r1.fastq.gz
fastqc --nogroup --outdir 1_fastqc/ 0_reads/pviburni.1821.0B.350.r1.fastq.gz
fastqc --nogroup --outdir 1_fastqc/ 0_reads/pviburni.1823.0B.350.r1.fastq.gz
fastqc --nogroup --outdir 1_fastqc/ 0_reads/pviburni.184.2B.350.r1.fastq.gz
fastqc --nogroup --outdir 1_fastqc/ 0_reads/pviburni.1813.1B.350.r2.fastq.gz
fastqc --nogroup --outdir 1_fastqc/ 0_reads/pviburni.1813.1B.550.r2.fastq.gz
fastqc --nogroup --outdir 1_fastqc/ 0_reads/pviburni.1821.0B.350.r2.fastq.gz
fastqc --nogroup --outdir 1_fastqc/ 0_reads/pviburni.1823.0B.350.r2.fastq.gz
fastqc --nogroup --outdir 1_fastqc/ 0_reads/pviburni.184.2B.350.r2.fastq.gz
```
submission command
```
parallel -j1 'qsub -cwd -N fastqc -V -pe **smp64** 16 -b yes {}' :::: 1_fastqc/fastqc.sub
```
OK

## trimming with fastp
submission filename: fastp.sub
```
fastp -i 0_reads/pviburni.1813.1B.350.r1.fastq.gz -I 0_reads/pviburni.1813.1B.350.r2.fastq.gz -o 2_trim/pviburni.1813.1B.350.r1.trim.fastq.gz -O 2_trim/pviburni.1813.1B.350.r2.trim.fastq.gz --cut_by_quality5 --cut_by_quality3 --cut_window_size 4 --cut_mean_quality 20 --html 2_trim/2_trim/pviburni.1813.1B.350.html --thread 4
fastp -i 0_reads/pviburni.1813.1B.550.r1.fastq.gz -I 0_reads/pviburni.1813.1B.550.r2.fastq.gz -o 2_trim/pviburni.1813.1B.550.r1.trim.fastq.gz -O 2_trim/pviburni.1813.1B.550.r2.trim.fastq.gz --cut_by_quality5 --cut_by_quality3 --cut_window_size 4 --cut_mean_quality 20 --html 2_trim/2_trim/pviburni.1813.1B.550.html --thread 4
fastp -i 0_reads/pviburni.1821.0B.350.r1.fastq.gz -I 0_reads/pviburni.1821.0B.350.r2.fastq.gz -o 2_trim/pviburni.1821.0B.350.r1.trim.fastq.gz -O 2_trim/pviburni.1821.0B.350.r2.trim.fastq.gz --cut_by_quality5 --cut_by_quality3 --cut_window_size 4 --cut_mean_quality 20 --html 2_trim/2_trim/pviburni.1821.0B.350.html --thread 4
fastp -i 0_reads/pviburni.1823.0B.350.r1.fastq.gz -I 0_reads/pviburni.1823.0B.350.r2.fastq.gz -o 2_trim/pviburni.1823.0B.350.r1.trim.fastq.gz -O 2_trim/pviburni.1823.0B.350.r2.trim.fastq.gz --cut_by_quality5 --cut_by_quality3 --cut_window_size 4 --cut_mean_quality 20 --html 2_trim/2_trim/pviburni.1823.0B.350.html --thread 4
fastp -i 0_reads/pviburni.184.2B.350.r1.fastq.gz -I 0_reads/pviburni.184.2B.350.r2.fastq.gz -o 2_trim/pviburni.184.2B.350.r1.trim.fastq.gz -O 2_trim/pviburni.184.2B.350.r2.trim.fastq.gz --cut_by_quality5 --cut_by_quality3 --cut_window_size 4 --cut_mean_quality 20 --html 2_trim/2_trim/pviburni.184.2B.350.html --thread 4

```
submit to parallel
```
parallel -j1 'qsub -cwd -N fastp -V -pe **smp64** 16 -b yes {}' :::: 2_trim/fastp.sub

```


## fastqc after trimming
subfile: fastqc.trim.sub
```
fastqc --nogroup --outdir 3_fastqc/ 2_trim/pviburni.1813.1B.350.r1.trim.fastq.gz
fastqc --nogroup --outdir 3_fastqc/ 2_trim/pviburni.1813.1B.550.r1.trim.fastq.gz
fastqc --nogroup --outdir 3_fastqc/ 2_trim/pviburni.1821.0B.350.r1.trim.fastq.gz
fastqc --nogroup --outdir 3_fastqc/ 2_trim/pviburni.1823.0B.350.r1.trim.fastq.gz
fastqc --nogroup --outdir 3_fastqc/ 2_trim/pviburni.184.2B.350.r1.trim.fastq.gz
fastqc --nogroup --outdir 3_fastqc/ 2_trim/pviburni.1813.1B.350.r2.trim.fastq.gz
fastqc --nogroup --outdir 3_fastqc/ 2_trim/pviburni.1813.1B.550.r2.trim.fastq.gz
fastqc --nogroup --outdir 3_fastqc/ 2_trim/pviburni.1821.0B.350.r2.trim.fastq.gz
fastqc --nogroup --outdir 3_fastqc/ 2_trim/pviburni.1823.0B.350.r2.trim.fastq.gz
fastqc --nogroup --outdir 3_fastqc/ 2_trim/pviburni.184.2B.350.r2.trim.fastq.gz
```
submission command
```
parallel -j1 'qsub -cwd -N fastqc -V -pe smp64 32 -b yes {}' :::: 3_fastqc/fastqc.trim.sub
```




## CLC (only runs on bigfoot)

submission file: clc.sub
assemble all files in one assembly with single end option
```
/ceph/software/clc/clc-assembly-cell-5.0.0-linux_64/clc_assembler -q 2_trim/pviburni.1813.1B.350.r1.trim.fastq.gz 2_trim/pviburni.1813.1B.550.r1.trim.fastq.gz 2_trim/pviburni.1821.0B.350.r1.trim.fastq.gz 2_trim/pviburni.1823.0B.350.r1.trim.fastq.gz 2_trim/pviburni.184.2B.350.r1.trim.fastq.gz 2_trim/pviburni.1813.1B.350.r2.trim.fastq.gz 2_trim/pviburni.1813.1B.550.r2.trim.fastq.gz 2_trim/pviburni.1821.0B.350.r2.trim.fastq.gz 2_trim/pviburni.1823.0B.350.r2.trim.fastq.gz 2_trim/pviburni.184.2B.350.r2.trim.fastq.gz --cpus 24 -o 4_clc/pviburni.clc.se.fna

```

```
parallel -j1 'qsub -cwd -l h=bigfoot -N clc -V -pe **smp64** 32 -b yes {}' :::: 4_clc/clc.sub
```

clc assembly finished -> pviburni.clc.se.fna

## basic summary of clc assembly using QUAST (ran locally)

```
vpn2-083:quast-4.6.3 isabelle$ python quast.py -o ~/Dropbox/UniversityofEdinburgh/BNPGE/Deliverables/Genomics/ ~/Dropbox/UniversityofEdinburgh/BNPGE/Deliverables/Genomics/pviburni.clc.se.fna
/Users/isabelle/quast-4.6.3/quast.py -o /Users/isabelle/Dropbox/UniversityofEdinburgh/BNPGE/Deliverables/Genomics/ /Users/isabelle/Dropbox/UniversityofEdinburgh/BNPGE/Deliverables/Genomics/pviburni.clc.se.fna

```
```
Assembly                    pviburni.clc.se
# contigs (>= 0 bp)         740451         
# contigs (>= 1000 bp)      149617         
# contigs (>= 5000 bp)      4693           
# contigs (>= 10000 bp)     761            
# contigs (>= 25000 bp)     44             
# contigs (>= 50000 bp)     7              
Total length (>= 0 bp)      562345655      
Total length (>= 1000 bp)   296631269      
Total length (>= 5000 bp)   37593384       
Total length (>= 10000 bp)  11866022       
Total length (>= 25000 bp)  2162149        
Total length (>= 50000 bp)  1037753        
# contigs                   341580         
Largest contig              279562         
Total length                429348924      
GC (%)                      33.02          
N50                         1474           
N75                         876            
L50                         81892          
L75                         176765         
# N's per 100 kbp           0.00
```
# II. BLOBPLOTS (copied from Quiver 2018)

source activate blobtools_env
## Make index of the assembly with bwa
This was run in qmaster
```
bwa index pviburni.clc.se.fna
```
## get hits by blastn and blastx through diamond
Diamond is a program that is faster when running a blastx on genome data
subfile name: blast_commands.txt
```
blastn -query 4_clc/pviburni.clc.se.fna -num_threads 24 -db /ceph/software/blast_db/blast/ncbi_2018_01/nt -evalue 1e-25 -outfmt '6 qseqid staxids bitscore std' -out 4_clc/pviburni.clc.se.vs.nt.mts10.hsp1.out -max_target_seqs 10 -max_hsps 1

/ceph/software/diamond/diamond-v0.9.17/diamond blastx --threads 24 --query 4_clc/pviburni.clc.se.fna --db
/ceph/software/blast_db/diamond/uniprot_2017_12/uniprot_ref_proteomes.dmnd --sensitive --evalue 1e-25 --outfmt 6 --out 4_clc/pviburni.clc.se.vs.uniref.mts5.out --max-target-seqs 5

```

submitted to main.q (smp)
```
parallel -j1 'qsub -cwd -N blast -V -pe smp 24 -b yes {}' :::: blast_commands.txt
```

## mapping reads to the assembly using bwa
submission file: bwa_commands.txt
```
bwa mem -t 16 4_clc/pviburni.clc.se.fna 2_trim/pviburni.1813.1B.350.r1.trim.fastq.gz 2_trim/pviburni.1813.1B.350.r2.trim.fastq.gz | samtools view -b - > /scratch/ivea/pviburni.1813.1B.350.vs.pviburni.clc.se.bam && rsync /scratch/ivea/pviburni.1813.1B.350.vs.pviburni.clc.se.bam . && touch pviburni.1813.1B.350.vs.pviburni.clc.se.bam.done

bwa mem -t 16 4_clc/pviburni.clc.se.fna 2_trim/pviburni.1813.1B.550.r1.trim.fastq.gz 2_trim/pviburni.1813.1B.550.r2.trim.fastq.gz | samtools view -b - > /scratch/ivea/pviburni.1813.1B.550.vs.pviburni.clc.se.bam && rsync /scratch/ivea/pviburni.1813.1B.550.vs.pviburni.clc.se.bam . && touch pviburni.1813.1B.550.vs.pviburni.clc.se.bam.done
bwa mem -t 16 4_clc/pviburni.clc.se.fna 2_trim/pviburni.1821.0B.350.r1.trim.fastq.gz 2_trim/pviburni.1821.0B.350.r2.trim.fastq.gz | samtools view -b - > /scratch/ivea/pviburni.1821.0B.350.vs.pviburni.clc.se.bam && rsync /scratch/ivea/pviburni.1821.0B.350.vs.pviburni.clc.se.bam . && touch pviburni.1821.0B.350.vs.pviburni.clc.se.bam.done
bwa mem -t 16 4_clc/pviburni.clc.se.fna 2_trim/pviburni.1823.0B.350.r1.trim.fastq.gz 2_trim/pviburni.1823.0B.350.r2.trim.fastq.gz | samtools view -b - > /scratch/ivea/pviburni.1823.0B.350.vs.pviburni.clc.se.bam && rsync /scratch/ivea/pviburni.1823.0B.350.vs.pviburni.clc.se.bam . && touch pviburni.1823.0B.350.vs.pviburni.clc.se.bam.done
bwa mem -t 16 4_clc/pviburni.clc.se.fna 2_trim/pviburni.184.2B.350.r1.trim.fastq.gz 2_trim/pviburni.184.2B.350.r2.trim.fastq.gz | samtools view -b - > /scratch/ivea/pviburni.184.2B.350.vs.pviburni.clc.se.bam && rsync --remove-source-files /scratch/ivea/pviburni.184.2B.350.vs.pviburni.clc.se.bam . && touch pviburni.184.2B.350.vs.pviburni.clc.se.bam.done
```
run bwa for each library against the assembly and writes in /scratch, then rsync to qmaster, then creates an empty file.done to indicate that it is finished and the rsync was made

```
parallel -j1 'qsub -cwd -N bwa -V -pe smp 16 -b yes {}' :::: bwa_commands.txt
```

### Taxify on diamond output file
Allows to process the ID that is not recognisable from Diamond.
Run in qmaster directly => not done yet
```
/ceph/software/blobtools/blobtools taxify -f 4_clc/pviburni.clc.se.vs.uniref.mts5.out -m /ceph/software/blast_db/diamond/uniprot_2017_12/uniprot_ref_proteomes.taxids -s 0 -t 2
```


#5 July 2018
## Do blobtools in bigwig

ssh bigwig (from qm)

source activate blobtools

```
/ceph/software/blobtools/blobtools create -i 4_clc/pviburni.clc.se.fna -b pviburni.1813.1B.350.vs.pviburni.clc.se.bam -b pviburni.1813.1B.550.vs.pviburni.clc.se.bam -b pviburni.1821.0B.350.vs.pviburni.clc.se.bam -b pviburni.1823.0B.350.vs.pviburni.clc.se.bam -b pviburni.184.2B.350.vs.pviburni.clc.se.bam -t 4_clc/pviburni.clc.se.vs.nt.mts10.hsp1.out -t 4_clc/pviburni.clc.se.vs.uniref.mts5.taxified.out -x bestsumorder -o 5_blobtools/pviburni
```
running the command line directly in bigwig without submission

json db is here: /data/ross/mealybugs/analyses/pseudococcus_viburni/5_blobtools/pviburni.blobDB.json


```
/ceph/software/blobtools/blobtools plot -i pviburni.blobDB.json -x bestsumorder
```

```
/ceph/software/blobtools/blobtools view -i pviburni.blobDB.json -x bestsumorder --hits --rank all
```


# IV. We now need to map the contigs from the Illumina candidates to the PacBio assembly (From here, done in 2020 - 2021)
### December 2020
I chose List L for to map candidate contigs. Reminder, list L was determined by applying thresholds to the 5 libraries as follows:
- $5 2B >100
- $6 1B >25
- $7 1B >25
- $8 0B <1
- $9 0B <1
No taxa filter was applied.

1. Make a file of sequences containing the list of contigs from list L

List file: /data/ross/mealybugs/analyses/B_viburni_2020/isabelle_old_2019/pseudococcus_viburni/5_blobtools/1_blobplot1
Fasta file of contigs: /data/ross/mealybugs/analyses/B_viburni_2020/isabelle_old_2019/pseudococcus_viburni/4_clc
```
awk -F'>' 'NR==FNR{ids[$0]; next} NF>1{f=($2 in ids)} f' listL.txt /data/ross/mealybugs/analyses/B_viburni_2020/isabelle_old_2019/pseudococcus_viburni/4_clc/pviburni.clc.se.fna > pviburni.clc.se.listL.fna
```


2. Mapping the B contigs from list L to the assembly


### January 2021
=====January 15 2020=====
Working directory and environment
```
cd /data/ross/mealybugs/analyses/B_viburni_2020/isabelle_old_2019/pseudococcus_viburni/5_blobtools/1_blobplot1/

conda activate /ceph/users/afilia/.conda/envs/afilia
```
This directory contains the reads for each libraries and lists of B candidates.

- My contig files extract from list L: /data/ross/mealybugs/analyses/B_viburni_2020/isabelle_old_2019/pseudococcus_viburni/5_blobtools/1_blobplot1/pviburni.clc.se.listL.fna

- PacBio assembly: /data/ross/mealybugs/analyses/B_viburni_2020/1_pacbio_assembly/5_assembly_freeze_v0/p.viburni.freeze.v0.fa

I will retry bwa mapping using Christina submission command and use nucmer as in section 5 of 3.Coverage_analysis

```
qsub -o logs -e logs -cwd -N nucmer -V -pe smp64 1 -b yes 'nucmer -p pviburni.clc.se.listL.p.viburni.freeze.v0.fa -t 24 /data/ross/mealybugs/analyses/B_viburni_2020/1_pacbio_assembly/5_assembly_freeze_v0/p.viburni.freeze.v0.fa /data/ross/mealybugs/analyses/B_viburni_2020/isabelle_old_2019/pseudococcus_viburni/5_blobtools/1_blobplot1/pviburni.clc.se.listL.fna'
```
ok
```
show-coords -clrT  pviburni.clc.se.listL.p.viburni.freeze.v0.fa.delta > pviburni.clc.se.listL.p.viburni.freeze.v0.fa.delta.coords
```

```
awk '{ a[$12]++ } END { for (b in a) { print b } }' pviburni.clc.se.listL.p.viburni.freeze.v0.fa.delta.coords > pviburni.clc.se.listL.p.viburni.freeze.v0.scaffolds.list
```
The contigs mapped to 255 scaffolds of the PacBio assembly.


3. Chromosome status of scaffolds found in mapping analysis

Now I want to see how in which categories they are sorted. I used Andres script with the files that have chromosome status. The R script is here:

The table of scaffolds that correspond to the contigs from Illimina genome sequencing and were assigned as B sequences from the coverage analysis.

3.1. Table of B strict genes (based on scaffold number only)
There are 9 scaffolds where candidate genes identified in the coverage analysis are found.




Scaffold_497 has 3 contigs from Illumina. List of corresponding contigs can be found in /output/pviburni.clc.se.listL.p.viburni.freeze.v0.fa.delta.coord.B3

3.4. Scaffold assigned to A chromosomes.

194 scaffolds were assigned to A. I didn't compile the list of contigs but this can be done is necessary.

3.5. Going from the contigs that where assigned to Nematoda or Arthropoda in blobtools.

To triple check, I took list L again and looked at the sequences that had a hit for Arthropoda and Nematoda. There are only XX contigs and I checked which scaffolds and chromosome status they were.

| contig Illumina   assembly | PacBio_scaffold | chromosome status |
|---------------------------|-----------------|-------------------|
| contig_14807              | scaffold_148    | A                 |
| contig_33424              | scaffold_1148    | B1                 |
| contig_45379              | scaffold_28     | A                 |
| contig_50667              | scaffold_94     | A                 |
| contig_67571              | scaffold_359    | A                 |
| contig_143760             | scaffold_1269   | B1                |
| contig_145684             | scaffold_552    | A                 |
| contig_187923             | NA              |                   |
| contig_246315             | scaffold_467    | A                 |
| contig_260819             | scaffold_552    | B1                |
| contig_264982             | scaffold_91     | A                 |
| contig_273720             | scaffold_1542   | NA                |
| contig_321423             | scaffold_467    | A                 |
| contig_356463             | scaffold_9      | A                 |
| contig_385605             | scaffold_552    | B1                |
| contig_443093             | scaffold_193    | A                 |
| contig_518017             | scaffold_28     | A                 |

Very few contigs are B!

I wonder if we should try to blast the contigs against the PacBio assembly genes and see if there are more information on the annotation?
