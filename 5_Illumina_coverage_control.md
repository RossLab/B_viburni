
# Illumina coverage analysis

	# working directory	
	/data/ross/mealybugs/analyses/B_viburni_2020/isabelle_old_2019/pseudococcus_viburni/
	# raw reads
	/data/ross/mealybugs/raw/XXXXX

# I. GENOME ASSEMBLY => to copy from Quiver

## 1. Raw RNAseq reads

## 1. Raw RNAseq reads


## 1. Raw RNAseq reads


## 1. Raw RNAseq reads

# II. BLOBPLOTS ==> to copy from Quiver




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


