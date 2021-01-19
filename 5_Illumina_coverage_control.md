
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


#########
##July 10 2018

###Examine blobplot and start building a list for only the contigs I want to keep
read only one column of interest
```
cut -f15 pviburni.blobDB.bestsumorder.table.txt
```
Showing all instances of one column - useful for taxa filtering
```
cut -f15 pviburni.blobDB.bestsumorder.table.txt | sort | uniq -c
```

subsets the table by whatever factors 
```
awk '$15=="Arthropoda"' pviburni.blobDB.bestsumorder.table.txt 
```

####A - contigs top left side
GC proportion <0.2 and = GC ==> $3
coverage > 8000 and ===> cov_sum   = $10
Arthropods = phylum.t.15 ==> $15
```
awk '$3 < 0.2 && $10 > 10000 && $15=="Arthropoda" ' pviburni.blobDB.bestsumorder.table.txt >  pviburni.blobDB.bestsumorder.table.A.txt
```




# III. Lists of candidates by coverage
Using: pviburni.blobDB.bestsumorder.table.txt

Column numbers are:
- $5 2B
- $6 1B
- $7 1B
- $8 0B
- $9 0B

I used different threshold of coverage for each library and filtered by taxa or not.
The lists that are relevant:

## List D - Potential candidates of B sequences in the Arthropod hits
Very conservative criteria: need to be Arthropod and coverage as follow
- $5 2B >100
- $6 1B >100
- $7 1B >100
- $8 0B <10
- $9 0B <10
```
awk '$15 =="Arthropoda" && $5 > 100 && $6 >100  && $7 >100 && $8 <10 && $9 < 10 >'  pviburni.blobDB.bestsumorder.table.txt >  pviburni.blobDB.bestsumorder.table.D.txt 
```
## List L
No Arthropod filtering and coverage as follows
- $5 2B >100
- $6 1B >25
- $7 1B >25
- $8 0B <1
- $9 0B <1

```
awk '$5 > 100 && $6 >25  && $7 >25 && $8 <1 && $9 <1'  pviburni.blobDB.bestsumorder.table.txt >  pviburni.blobDB.bestsumorder.table.L.txt 
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

/ceph/users/afilia/.conda/envs/afilia

1.1. using bwa
In submission file: 5_blobtools/1_blobplot1/bwa_listL_p.viburni.freeze.v0.sub
```
bwa ##mem -R '@RG\tID:WYE3_pacbio\tSM:WYE3_pacbio'## -x nanoseq -t 32 /data/ross/mealybugs/analyses/B_viburni_2020/1_pacbio_assembly/5_freeze_v0/p.viburni.freeze.v0.fa 5_blobtools/1_blobplot1/pviburni.clc.se.listL.fna | samtools view -bS - > 5_blobtools/1_blobplot1/pviburni.clc.se.listL.fna.vs.p.viburni.freeze.v0.fa
```

parallel -j1 'qsub -cwd -N bwa -V -pe **smp64** 32 -b yes {}' :::: 5_blobtools/1_blobplot1/bwa_listL_p.viburni.freeze.v0.sub

===> stopped here because parallel not working. 

Christina sent this command to submit to the cluster

qsub -o logs -e logs -cwd -N td.lorf -V -pe smp64 1 -b yes 'TransDecoder.LongOrfs -t transcriptomes/cech.transcriptome.longiso.fasta --output_dir 2_transdecoder/cech/'

```
qsub -o logfiles -e logfiles -cwd -N td.lorf -V -pe smp64 1 -b yes 'bwa ##mem -R '@RG\tID:WYE3_pacbio\tSM:WYE3_pacbio'## -x nanoseq -t 32 /data/ross/mealybugs/analyses/B_viburni_2020/1_pacbio_assembly/5_freeze_v0/p.viburni.freeze.v0.fa 5_blobtools/1_blobplot1/pviburni.clc.se.listL.fna | samtools view -bS - > 5_blobtools/1_blobplot1/pviburni.clc.se.listL.fna.vs.p.viburni.freeze.v0.fa'
```

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

There are 255 scaffolds that mapped to the PacBio genome

3. Chromosome status of scaffolds found in mapping analysis

Now I want to see how in which categories they are sorted. I used Andres script with the files that have chromosome status. The R script is here: 

The table of scaffolds that correspond to the contigs from Illimina genome sequencing and were assigned as B sequences from the coverage analysis.

3.1. Table of B strict genes (based on scaffold number only)
There are 9 scaffolds where candidate genes identified in the coverage analysis are found.


| seq           | length | b.status.final | cov.04v13  | b.status | b.status.asn | b.status.kmer | gene   | gene_len | blast       | diamond    | pfam_acc        | pfam_descr                                                                                 | ipr_acc             | ipr_descr                                                          | GO                                          | anno |
|---------------|--------|----------------|------------|----------|--------------|---------------|--------|----------|-------------|------------|-----------------|--------------------------------------------------------------------------------------------|---------------------|--------------------------------------------------------------------|---------------------------------------------|------|
| scaffold_1124 | 20190  | B1             | 1.711547   | B.strict | A            | B             | g19189 | 2414     | NA          | NA         | NA              | NA                                                                                         | NA                  | NA                                                                 | NA                                          | N    |
| scaffold_1124 | 20190  | B1             | 1.711547   | B.strict | A            | B             | g19188 | 285      | NA          | NA         | NA              | NA                                                                                         | NA                  | NA                                                                 | NA                                          | N    |
| scaffold_1124 | 20190  | B1             | 1.711547   | B.strict | A            | B             | g19190 | 1792     | NA          | NA         | NA              | NA                                                                                         | NA                  | NA                                                                 | NA                                          | N    |
| scaffold_1148 | 19549  | B1             | 1.88146    | B.strict | A            | B             | g16323 | 725      | KIF23_HUMAN | K7J331     | PF00225         | Kinesin motor domain                                                                       | IPR001752           | Kinesin motor domain                                               | GO:0003777;GO:0005524;GO:0007018;GO:0008017 | Y    |
| scaffold_1148 | 19549  | B1             | 1.88146    | B.strict | A            | B             | g16324 | 797      | KIF23_HUMAN | K7J331     | PF00225         | Kinesin motor domain                                                                       | IPR001752           | Kinesin motor domain                                               | GO:0003777;GO:0005524;GO:0007018;GO:0008017 | Y    |
| scaffold_1269 | 17142  | B1             | 1.907414   | B.strict | A            | B             | g12159 | 1076     | NA          | A0A226CTZ9 | NA              | NA                                                                                         | NA                  | NA                                                                 | NA                                          | Y    |
| scaffold_1269 | 17142  | B1             | 1.907414   | B.strict | A            | B             | g12160 | 1038     | NA          | A0A226DDY8 | PF05699         | hAT family C-terminal dimerisation region                                                  | IPR008906           | HAT, C-terminal dimerisation domain                                | GO:0046983                                  | Y    |
| scaffold_1269 | 17142  | B1             | 1.907414   | B.strict | A            | B             | g12161 | 369      | NA          | NA         | NA              | NA                                                                                         | NA                  | NA                                                                 | NA                                          | N    |
| scaffold_1271 | 17114  | B1             | 2.013767   | B.strict | A            | B             | g4788  | 222      | NA          | NA         | NA              | NA                                                                                         | NA                  | NA                                                                 | NA                                          | N    |
| scaffold_1768 | 9394   | B1             | 0.4527596  | B.strict | B            | B             | g23391 | 2182     | NA          | J9LQI3     | NA              | NA                                                                                         | NA                  | NA                                                                 | NA                                          | Y    |
| scaffold_552  | 189883 | B1             | 1.94663467 | B.strict | A            | B             | g20107 | 596      | NA          | NA         | NA              | NA                                                                                         | NA                  | NA                                                                 | NA                                          | N    |
| scaffold_552  | 189883 | B1             | 1.94663467 | B.strict | A            | B             | g20086 | 742      | NA          | NA         | NA              | NA                                                                                         | NA                  | NA                                                                 | NA                                          | N    |
| scaffold_552  | 189883 | B1             | 1.94663467 | B.strict | A            | B             | g20097 | 552      | NA          | T1HLH9     | NA              | NA                                                                                         | NA                  | NA                                                                 | NA                                          | Y    |
| scaffold_552  | 189883 | B1             | 1.94663467 | B.strict | A            | B             | g20084 | 1170     | NA          | A0A0J7N8Y6 | PF17921         | Integrase zinc binding domain                                                              | IPR041588           | Integrase zinc-binding domain                                      | NA                                          | Y    |
| scaffold_552  | 189883 | B1             | 1.94663467 | B.strict | A            | B             | g20087 | 737      | NA          | NA         | NA              | NA                                                                                         | NA                  | NA                                                                 | NA                                          | N    |
| scaffold_552  | 189883 | B1             | 1.94663467 | B.strict | A            | B             | g20095 | 928      | NA          | NA         | NA              | NA                                                                                         | NA                  | NA                                                                 | NA                                          | N    |
| scaffold_552  | 189883 | B1             | 1.94663467 | B.strict | A            | B             | g20111 | 965      | NA          | NA         | PF07679         | Immunoglobulin I-set domain                                                                | IPR013098           | Immunoglobulin I-set                                               | NA                                          | Y    |
| scaffold_552  | 189883 | B1             | 1.94663467 | B.strict | A            | B             | g20106 | 1316     | YMD2_CAEEL  | A0A151IDY8 | NA              | NA                                                                                         | NA                  | NA                                                                 | NA                                          | Y    |
| scaffold_552  | 189883 | B1             | 1.94663467 | B.strict | A            | B             | g20089 | 1599     | NA          | NA         | NA              | NA                                                                                         | NA                  | NA                                                                 | NA                                          | N    |
| scaffold_552  | 189883 | B1             | 1.94663467 | B.strict | A            | B             | g20114 | 569      | NA          | NA         | NA              | NA                                                                                         | NA                  | NA                                                                 | NA                                          | N    |
| scaffold_552  | 189883 | B1             | 1.94663467 | B.strict | A            | B             | g20096 | 526      | NA          | NA         | NA              | NA                                                                                         | NA                  | NA                                                                 | NA                                          | N    |
| scaffold_552  | 189883 | B1             | 1.94663467 | B.strict | A            | B             | g20113 | 227      | NA          | NA         | NA              | NA                                                                                         | NA                  | NA                                                                 | NA                                          | N    |
| scaffold_552  | 189883 | B1             | 1.94663467 | B.strict | A            | B             | g20103 | 842      | NA          | J9L5G6     | NA              | NA                                                                                         | NA                  | NA                                                                 | NA                                          | Y    |
| scaffold_552  | 189883 | B1             | 1.94663467 | B.strict | A            | B             | g20083 | 563      | NA          | NA         | NA              | NA                                                                                         | NA                  | NA                                                                 | NA                                          | N    |
| scaffold_552  | 189883 | B1             | 1.94663467 | B.strict | A            | B             | g20088 | 677      | NA          | NA         | NA              | NA                                                                                         | NA                  | NA                                                                 | NA                                          | N    |
| scaffold_552  | 189883 | B1             | 1.94663467 | B.strict | A            | B             | g20104 | 1025     | NA          | NA         | NA              | NA                                                                                         | NA                  | NA                                                                 | NA                                          | N    |
| scaffold_552  | 189883 | B1             | 1.94663467 | B.strict | A            | B             | g20099 | 1158     | NA          | NA         | NA              | NA                                                                                         | NA                  | NA                                                                 | NA                                          | N    |
| scaffold_552  | 189883 | B1             | 1.94663467 | B.strict | A            | B             | g20090 | 686      | NA          | NA         | PF03564         | Protein of unknown function (DUF1759)                                                      | IPR005312           | Protein of unknown function DUF1759                                | NA                                          | Y    |
| scaffold_552  | 189883 | B1             | 1.94663467 | B.strict | A            | B             | g20091 | 1279     | NA          | NA         | NA              | NA                                                                                         | NA                  | NA                                                                 | NA                                          | N    |
| scaffold_552  | 189883 | B1             | 1.94663467 | B.strict | A            | B             | g20112 | 971      | NA          | NA         | PF00047         | Immunoglobulin domain                                                                      | IPR013151           | Immunoglobulin                                                     | NA                                          | Y    |
| scaffold_552  | 189883 | B1             | 1.94663467 | B.strict | A            | B             | g20108 | 2547     | NA          | A0A0R3QAL5 | PF00665         | Integrase core domain                                                                      | IPR001584           | Integrase, catalytic core                                          | GO:0015074                                  | Y    |
| scaffold_552  | 189883 | B1             | 1.94663467 | B.strict | A            | B             | g20085 | 1087     | NA          | NA         | NA              | NA                                                                                         | NA                  | NA                                                                 | NA                                          | N    |
| scaffold_552  | 189883 | B1             | 1.94663467 | B.strict | A            | B             | g20092 | 2589     | NA          | A0A0R3QBG2 | PF14529;PF00078 | Endonuclease-reverse transcriptase; Reverse transcriptase (RNA-dependent   DNA polymerase) | IPR005135;IPR000477 | Endonuclease/exonuclease/phosphatase; Reverse transcriptase domain | NA                                          | Y    |
| scaffold_552  | 189883 | B1             | 1.94663467 | B.strict | A            | B             | g20105 | 734      | NA          | NA         | NA              | NA                                                                                         | NA                  | NA                                                                 | NA                                          | N    |
| scaffold_552  | 189883 | B1             | 1.94663467 | B.strict | A            | B             | g20093 | 1501     | NA          | NA         | NA              | NA                                                                                         | NA                  | NA                                                                 | NA                                          | N    |
| scaffold_552  | 189883 | B1             | 1.94663467 | B.strict | A            | B             | g20110 | 641      | NA          | D6WQK2     | NA              | NA                                                                                         | NA                  | NA                                                                 | NA                                          | Y    |
| scaffold_552  | 189883 | B1             | 1.94663467 | B.strict | A            | B             | g20098 | 1001     | NA          | NA         | NA              | NA                                                                                         | NA                  | NA                                                                 | NA                                          | N    |
| scaffold_552  | 189883 | B1             | 1.94663467 | B.strict | A            | B             | g20102 | 310      | NA          | NA         | NA              | NA                                                                                         | NA                  | NA                                                                 | NA                                          | N    |
| scaffold_552  | 189883 | B1             | 1.94663467 | B.strict | A            | B             | g20094 | 405      | NA          | NA         | NA              | NA                                                                                         | NA                  | NA                                                                 | NA                                          | N    |
| scaffold_552  | 189883 | B1             | 1.94663467 | B.strict | A            | B             | g20100 | 3042     | NA          | A0A087T926 | PF07679;PF10551 | Immunoglobulin I-set domain; MULE transposase domain                                       | IPR013098;IPR018289 | Immunoglobulin I-set; MULE transposase domain                      | NA                                          | Y    |
| scaffold_552  | 189883 | B1             | 1.94663467 | B.strict | A            | B             | g20101 | 475      | NA          | NA         | NA              | NA                                                                                         | NA                  | NA                                                                 | NA                                          | N    |
| scaffold_552  | 189883 | B1             | 1.94663467 | B.strict | A            | B             | g20109 | 341      | NA          | NA         | NA              | NA                                                                                         | NA                  | NA                                                                 | NA                                          | N    |
| scaffold_552  | 189883 | B1             | 1.94663467 | B.strict | A            | B             | g20115 | 482      | NA          | J9KK61     | NA              | NA                                                                                         | NA                  | NA                                                                 | NA                                          | Y    |
| scaffold_786  | 48770  | B1             | 1.38485149 | B.strict | A            | B             | g2644  | 4614     | PGBD4_HUMAN | X1X046     | PF13843         | Transposase IS4                                                                            | IPR029526           | PiggyBac transposable element-derived protein                      | NA                                          | Y    |
| scaffold_902  | 30666  | B1             | 1.82721635 | B.strict | B            | B             | g19061 | 1061     | NA          | K7JWY9     | PF00078         | Reverse transcriptase (RNA-dependent DNA polymerase)                                       | IPR000477           | Reverse transcriptase domain                                       | NA                                          | Y    |
| scaffold_932  | 28592  | B1             | 2.02197511 | B.strict | A            | B             | g23135 | 122      | NA          | NA         | NA              | NA                                                                                         | NA                  | NA                                                                 | NA                                          | N    |


The list of contig numbers (Illumina assembly) that mapped to these scaffolds are in output/pviburni.clc.se.listL.p.viburni.freeze.v0.fa.delta.coord.Bstrict

Scaffolds that had the most contigs are scaffold_1148, scaffold_1271, scaffold_552 and scaffold_902.
NB: contig_260819 (that contains a chromo domain) is in scaffold_552 and contig_33424 is in scaffold_1148. 

3.2. List of scaffolds with B2 status

| seq           | length | b.status.final | cov.04v13  | b.status | b.status.asn | b.status.kmer | gene   | gene_len | blast       | diamond    | pfam_acc                                                | pfam_descr                                                                                                                                                                                                                                                                                   | ipr_acc                                                               | ipr_descr                                                                                                                                                                                                                                                                                                             | GO                                                                                      | anno |
|---------------|--------|----------------|------------|----------|--------------|---------------|--------|----------|-------------|------------|---------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------|------|
| scaffold_1330 | 16017  | B2             | 0.3300829  | B.loose  | A            | A             | g22793 | 330      | NA          | NA         | NA                                                      | NA                                                                                                                                                                                                                                                                                           | NA                                                                    | NA                                                                                                                                                                                                                                                                                                                    | NA                                                                                      | N    |
| scaffold_1330 | 16017  | B2             | 0.3300829  | B.loose  | A            | A             | g22792 | 103      | NA          | NA         | NA                                                      | NA                                                                                                                                                                                                                                                                                           | NA                                                                    | NA                                                                                                                                                                                                                                                                                                                    | NA                                                                                      | N    |
| scaffold_1637 | 11086  | B2             | 0.08185333 | B.loose  | A            | A             | g18895 | 224      | NA          | NA         | NA                                                      | NA                                                                                                                                                                                                                                                                                           | NA                                                                    | NA                                                                                                                                                                                                                                                                                                                    | NA                                                                                      | N    |
| scaffold_2062 | 6708   | B2             | 0.2535727  | B.loose  | A            | A             | g12685 | 847      | NA          | NA         | NA                                                      | NA                                                                                                                                                                                                                                                                                           | NA                                                                    | NA                                                                                                                                                                                                                                                                                                                    | NA                                                                                      | N    |
| scaffold_295  | 464503 | B2             | 0.44526067 | B.loose  | A            | A             | g17787 | 700      | NA          | NA         | NA                                                      | NA                                                                                                                                                                                                                                                                                           | NA                                                                    | NA                                                                                                                                                                                                                                                                                                                    | NA                                                                                      | N    |
| scaffold_295  | 464503 | B2             | 0.44526067 | B.loose  | A            | A             | g17791 | 1933     | NA          | A0A0R3Q549 | NA                                                      | NA                                                                                                                                                                                                                                                                                           | NA                                                                    | NA                                                                                                                                                                                                                                                                                                                    | NA                                                                                      | Y    |
| scaffold_295  | 464503 | B2             | 0.44526067 | B.loose  | A            | A             | g17807 | 1413     | NA          | NA         | NA                                                      | NA                                                                                                                                                                                                                                                                                           | NA                                                                    | NA                                                                                                                                                                                                                                                                                                                    | NA                                                                                      | N    |
| scaffold_295  | 464503 | B2             | 0.44526067 | B.loose  | A            | A             | g17804 | 391      | LSM4_HUMAN  | D6WJA0     | PF01423                                                 | LSM domain                                                                                                                                                                                                                                                                                   | IPR001163                                                             | LSM domain, eukaryotic/archaea-type                                                                                                                                                                                                                                                                                   | NA                                                                                      | Y    |
| scaffold_295  | 464503 | B2             | 0.44526067 | B.loose  | A            | A             | g17802 | 433      | CNI_DROME   | A0A482X1C0 | PF03311                                                 | Cornichon protein                                                                                                                                                                                                                                                                            | IPR003377                                                             | Cornichon                                                                                                                                                                                                                                                                                                             | GO:0016192                                                                              | Y    |
| scaffold_295  | 464503 | B2             | 0.44526067 | B.loose  | A            | A             | g17788 | 1377     | NA          | NA         | NA                                                      | NA                                                                                                                                                                                                                                                                                           | NA                                                                    | NA                                                                                                                                                                                                                                                                                                                    | NA                                                                                      | N    |
| scaffold_295  | 464503 | B2             | 0.44526067 | B.loose  | A            | A             | g17796 | 779      | NA          | NA         | NA                                                      | NA                                                                                                                                                                                                                                                                                           | NA                                                                    | NA                                                                                                                                                                                                                                                                                                                    | NA                                                                                      | N    |
| scaffold_295  | 464503 | B2             | 0.44526067 | B.loose  | A            | A             | g17810 | 196      | NA          | NA         | NA                                                      | NA                                                                                                                                                                                                                                                                                           | NA                                                                    | NA                                                                                                                                                                                                                                                                                                                    | NA                                                                                      | N    |
| scaffold_295  | 464503 | B2             | 0.44526067 | B.loose  | A            | A             | g17795 | 629      | NA          | NA         | NA                                                      | NA                                                                                                                                                                                                                                                                                           | NA                                                                    | NA                                                                                                                                                                                                                                                                                                                    | NA                                                                                      | N    |
| scaffold_295  | 464503 | B2             | 0.44526067 | B.loose  | A            | A             | g17792 | 1241     | NA          | A0A0R3Q943 | NA                                                      | NA                                                                                                                                                                                                                                                                                           | NA                                                                    | NA                                                                                                                                                                                                                                                                                                                    | NA                                                                                      | Y    |
| scaffold_295  | 464503 | B2             | 0.44526067 | B.loose  | A            | A             | g17811 | 244      | NA          | NA         | NA                                                      | NA                                                                                                                                                                                                                                                                                           | NA                                                                    | NA                                                                                                                                                                                                                                                                                                                    | NA                                                                                      | N    |
| scaffold_295  | 464503 | B2             | 0.44526067 | B.loose  | A            | A             | g17813 | 1666     | ARHGC_MOUSE | A0A151JPR7 | PF00621;PF17838                                         | RhoGEF domain; PH domain                                                                                                                                                                                                                                                                     | IPR000219;IPR041020                                                   | Dbl homology (DH) domain; ARHGEF1-like, PH domain                                                                                                                                                                                                                                                                     | GO:0005089;GO:0035023                                                                   | Y    |
| scaffold_295  | 464503 | B2             | 0.44526067 | B.loose  | A            | A             | g17789 | 467      | NA          | NA         | PF13920                                                 | Zinc finger, C3HC4 type (RING finger)                                                                                                                                                                                                                                                        | NA                                                                    | NA                                                                                                                                                                                                                                                                                                                    | NA                                                                                      | Y    |
| scaffold_295  | 464503 | B2             | 0.44526067 | B.loose  | A            | A             | g17794 | 247      | NA          | NA         | NA                                                      | NA                                                                                                                                                                                                                                                                                           | NA                                                                    | NA                                                                                                                                                                                                                                                                                                                    | NA                                                                                      | N    |
| scaffold_295  | 464503 | B2             | 0.44526067 | B.loose  | A            | A             | g17814 | 215      | NA          | NA         | NA                                                      | NA                                                                                                                                                                                                                                                                                           | NA                                                                    | NA                                                                                                                                                                                                                                                                                                                    | NA                                                                                      | N    |
| scaffold_295  | 464503 | B2             | 0.44526067 | B.loose  | A            | A             | g17803 | 2832     | MCM4_DROME  | A0A067R0Z3 | PF00493;PF17207;PF17855;PF14551                         | MCM P-loop domain; MCM OB domain; MCM AAA-lid domain; MCM N-terminal   domain                                                                                                                                                                                                                | IPR001208;IPR033762;IPR041562;IPR027925                               | MCM domain; MCM OB domain; MCM, AAA-lid domain; MCM N-terminal domain                                                                                                                                                                                                                                                 | GO:0003677;GO:0005524;GO:0006270                                                        | Y    |
| scaffold_295  | 464503 | B2             | 0.44526067 | B.loose  | A            | A             | g17815 | 509      | NA          | NA         | NA                                                      | NA                                                                                                                                                                                                                                                                                           | NA                                                                    | NA                                                                                                                                                                                                                                                                                                                    | NA                                                                                      | N    |
| scaffold_295  | 464503 | B2             | 0.44526067 | B.loose  | A            | A             | g17799 | 318      | NA          | NA         | NA                                                      | NA                                                                                                                                                                                                                                                                                           | NA                                                                    | NA                                                                                                                                                                                                                                                                                                                    | NA                                                                                      | N    |
| scaffold_295  | 464503 | B2             | 0.44526067 | B.loose  | A            | A             | g17801 | 4870     | DYH7_HUMAN  | A0A482WYZ6 | PF12777;PF18198;PF03028;PF00104;PF18199;PF12781;PF00105 | Microtubule-binding stalk of dynein motor; Dynein heavy chain AAA lid   domain; Dynein heavy chain region D6 P-loop domain; Ligand-binding domain of   nuclear hormone receptor; Dynein heavy chain C-terminal domain; ATP-binding   dynein motor region; Zinc finger, C4 type (two domains) | IPR024743;IPR041658;IPR004273;IPR000536;IPR041228;IPR035706;IPR001628 | Dynein heavy chain, coiled coil stalk; Dynein heavy chain AAA lid domain;   Dynein heavy chain region D6 P-loop domain; Nuclear hormone receptor,   ligand-binding domain; Dynein heavy chain, C-terminal domain; Dynein heavy   chain, ATP-binding dynein motor region; Zinc finger, nuclear hormone   receptor-type | GO:0003777;GO:0007018;GO:0030286;GO:0003700;GO:0005634;GO:0006355;GO:0008270;GO:0043565 | Y    |
| scaffold_295  | 464503 | B2             | 0.44526067 | B.loose  | A            | A             | g17817 | 1903     | NA          | A0A2J7QCN0 | NA                                                      | NA                                                                                                                                                                                                                                                                                           | NA                                                                    | NA                                                                                                                                                                                                                                                                                                                    | NA                                                                                      | Y    |
| scaffold_295  | 464503 | B2             | 0.44526067 | B.loose  | A            | A             | g17816 | 376      | NA          | NA         | NA                                                      | NA                                                                                                                                                                                                                                                                                           | NA                                                                    | NA                                                                                                                                                                                                                                                                                                                    | NA                                                                                      | N    |
| scaffold_295  | 464503 | B2             | 0.44526067 | B.loose  | A            | A             | g17793 | 1052     | NA          | NA         | NA                                                      | NA                                                                                                                                                                                                                                                                                           | NA                                                                    | NA                                                                                                                                                                                                                                                                                                                    | NA                                                                                      | N    |
| scaffold_295  | 464503 | B2             | 0.44526067 | B.loose  | A            | A             | g17805 | 584      | NA          | NA         | NA                                                      | NA                                                                                                                                                                                                                                                                                           | NA                                                                    | NA                                                                                                                                                                                                                                                                                                                    | NA                                                                                      | N    |
| scaffold_295  | 464503 | B2             | 0.44526067 | B.loose  | A            | A             | g17812 | 509      | NA          | NA         | NA                                                      | NA                                                                                                                                                                                                                                                                                           | NA                                                                    | NA                                                                                                                                                                                                                                                                                                                    | NA                                                                                      | N    |
| scaffold_295  | 464503 | B2             | 0.44526067 | B.loose  | A            | A             | g17818 | 1395     | NXF1_RAT    | A0A2J7R106 | PF09162;PF03943                                         | Tap, RNA-binding; TAP C-terminal domain                                                                                                                                                                                                                                                      | IPR015245;IPR005637                                                   | Nuclear RNA export factor Tap, RNA-binding domain; TAP C-terminal (TAP-C)   domain                                                                                                                                                                                                                                    | GO:0003723;GO:0005634;GO:0005737;GO:0006406;GO:0051028                                  | Y    |
| scaffold_295  | 464503 | B2             | 0.44526067 | B.loose  | A            | A             | g17809 | 293      | NA          | NA         | NA                                                      | NA                                                                                                                                                                                                                                                                                           | NA                                                                    | NA                                                                                                                                                                                                                                                                                                                    | NA                                                                                      | N    |
| scaffold_295  | 464503 | B2             | 0.44526067 | B.loose  | A            | A             | g17806 | 645      | BHE22_XENTR | A0A1S3D3L0 | PF00010                                                 | Helix-loop-helix DNA-binding domain                                                                                                                                                                                                                                                          | IPR011598                                                             | Myc-type, basic helix-loop-helix (bHLH) domain                                                                                                                                                                                                                                                                        | GO:0046983                                                                              | Y    |
| scaffold_295  | 464503 | B2             | 0.44526067 | B.loose  | A            | A             | g17800 | 2289     | NA          | A0A151JC66 | PF03175                                                 | DNA polymerase type B, organellar and viral                                                                                                                                                                                                                                                  | IPR004868                                                             | DNA-directed DNA polymerase, family B, mitochondria/virus                                                                                                                                                                                                                                                             | GO:0000166;GO:0003677;GO:0003887;GO:0006260                                             | Y    |
| scaffold_295  | 464503 | B2             | 0.44526067 | B.loose  | A            | A             | g17797 | 638      | NA          | NA         | NA                                                      | NA                                                                                                                                                                                                                                                                                           | NA                                                                    | NA                                                                                                                                                                                                                                                                                                                    | NA                                                                                      | N    |
| scaffold_295  | 464503 | B2             | 0.44526067 | B.loose  | A            | A             | g17808 | 2723     | NA          | A0A0L7RH29 | PF09128;PF00130;PF00595                                 | Regulator of G protein signalling-like domain; Phorbol   esters/diacylglycerol binding domain (C1 domain); PDZ domain                                                                                                                                                                        | IPR015212;IPR002219;IPR001478                                         | Regulator of G protein signalling-like domain; Protein kinase C-like,   phorbol ester/diacylglycerol-binding domain; PDZ domain                                                                                                                                                                                       | GO:0005089;GO:0005737;GO:0035556;GO:0005515                                             | Y    |
| scaffold_295  | 464503 | B2             | 0.44526067 | B.loose  | A            | A             | g17790 | 988      | NA          | NA         | NA                                                      | NA                                                                                                                                                                                                                                                                                           | NA                                                                    | NA                                                                                                                                                                                                                                                                                                                    | NA                                                                                      | N    |
| scaffold_295  | 464503 | B2             | 0.44526067 | B.loose  | A            | A             | g17798 | 1807     | YMD2_CAEEL  | A0A151IDY8 | NA                                                      | NA                                                                                                                                                                                                                                                                                           | NA                                                                    | NA                                                                                                                                                                                                                                                                                                                    | NA                                                                                      | Y    |
| scaffold_607  | 138163 | B2             | 0.53066585 | B.loose  | A            | A             | g14106 | 2281     | BAB2_DROME  | A0A482WRV5 | PF00651;PF00096                                         | BTB/POZ domain; Zinc finger, C2H2 type                                                                                                                                                                                                                                                       | IPR000210;IPR013087                                                   | BTB/POZ domain; Zinc finger C2H2-type                                                                                                                                                                                                                                                                                 | GO:0005515;GO:0003676                                                                   | Y    |
| scaffold_607  | 138163 | B2             | 0.53066585 | B.loose  | A            | A             | g14107 | 470      | NA          | NA         | NA                                                      | NA                                                                                                                                                                                                                                                                                           | NA                                                                    | NA                                                                                                                                                                                                                                                                                                                    | NA                                                                                      | N    |
