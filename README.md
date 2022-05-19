# Viral_Metagenomics_Ants
## Dependencies 
Java
Perl
Bowtie2
samtools
Spades 3.14.0
trimmomatic-0.35 

## Workflow
This workflow starts with raw paired-end HiSeq data in demultiplexed FASTQ formated assumed to be located within a folder called raw_seq

## Concatenates the reads from multiple lanes of sequencing together  
```sh
cat /scratch/midway2/pflynn/15/15A_TTGCCTAG-TAAGTGGT-AHNFFJBBXX_L004_R1.fastq.gz /scratch/midway2/pflynn/15/15A_TTGCCTAG-TAAGTGGT-AHNFFJBBXX_L005_R1.fastq.gz /scratch/midway2/pflynn/15/15A_TTGCCTAG-TAAGTGGT-AHWYVLBBXX_L005_R1.fastq.gz > /scratch/midway2/pflynn/15/15A_concatenated_R1.fastq.gz
```
## Run fastqc for manual inspection of the quality of the sequences 
```sh
mkdir fastqc_out
fastqc -t 4 raw_data/* -o fastqc_out/
```
##trims the contactenated reads for adapters using Trimmomatic
```sh
java -jar /project2/mlcoleman/src/trimmomatic/Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 \
/scratch/midway2/pflynn/15/15A_concatenated_R1.fastq.gz /scratch/midway2/pflynn/15/15A_concatenated_R2.fastq.gz  \
/scratch/midway2/pflynn/15/15_pe1.fq /scratch/midway2/pflynn/15/15_unp1.fq /scratch/midway2/pflynn/15/15_pe2.fq /scratch/midway2/pflynn/15/15_unp2.fq \
ILLUMINACLIP:/project2/mlcoleman/src/trimmomatic/Trimmomatic-0.35/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```
##builds ant species specific genome database in bowtie2
```sh
bowtie2-build /scratch/midway2/pflynn/genomes/D_quadriceps.fna /scratch/midway2/pflynn/genomes/Dquad

```
##builds ant species specific genome database in bowtie2
```sh
bowtie2 -t -x /scratch/midway2/pflynn/genomes/Dquad -1 /scratch/midway2/pflynn/15/15_pe1.fq  -2 /scratch/midway2/pflynn/15/15_pe2.fq --un-conc /scratch/midway2/pflynn/15/15_conc_unmapped.fastq --al-conc /scratch/midway2/pflynn/15/15_conc_mapped.sam

```

