
# Illumina coverage analysis

	# working directory	
	/data/ross/mealybugs/analyses/B_viburni_2020/isabelle_old_2019/pseudococcus_viburni/
	# raw reads
	data/ross/mealybugs/analyses/B_viburni_2020/isabelle_old_2019/pseudococcus_viburni/0_reads/ 


# I. GENOME ASSEMBLY => to copy from Quiver



$ tmux attach-session -t genomics

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

# II. BLOBPLOTS ==> to copy from Quiver

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




# III. LISTs by coverage
Column numbers are
$5 2B
$6 1B
$7 1B
$8 0B
$9 0B

I used different threshold of coverage for each library and filtered by taxa or not.

## List D - Potential candidates of B sequences in the Arthropod hits
Very conservative criteria: need to be Arthropod and coverage as follow
$5 2B >100
$6 1B >100
$7 1B >100
$8 0B <10
$9 0B <10
```
awk '$15 =="Arthropoda" && $5 > 100 && $6 >100  && $7 >100 && $8 <10 && $9 < 10 >'  pviburni.blobDB.bestsumorder.table.txt >  pviburni.blobDB.bestsumorder.table.D.txt 
```
## List L
No Arthropod filtering and coverage as follows
$5 2B >100
$6 1B >25
$7 1B >25
$8 0B <1
$9 0B <1

```
awk '$5 > 100 && $6 >25  && $7 >25 && $8 <1 && $9 <1'  pviburni.blobDB.bestsumorder.table.txt >  pviburni.blobDB.bestsumorder.table.L.txt 
```

Using List L

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

1.2. using minimap

 minimap2 --secondary=no -ax map-pb -t 32 pseudococcus_viburni.hypo3.fa ../../0_reads/PV_18-13.1.subreads.fasta.gz ../../0_reads/PV_18-13.2.subreads.fasta.gz ../../0_reads/PV_18-13.3.subreads.fasta.gz | samtools view -hF 0x900 - | samtools sort -@32 -O BAM -o /scratch/afilia/pseudococcus_viburni.hypo3.sorted.bam - && rsync -av /scratch/afilia/pseudococcus_viburni.hypo3.sorted.bam .


