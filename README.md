# Ant Associated Viral Metagenomics Workflow
## Dependencies 
Java <br>
Perl <br>
Bowtie2 <br>
samtools <br>
trimmomatic-0.35 <br>
Spades 3.14.0 <br>
CheckV <br>

## Bioinformatics Workflow
This workflow starts with raw paired-end HiSeq data in demultiplexed FASTQ formated assumed to be located within a folder called raw_seq
Example for one sample
1. Concatenates the reads from multiple lanes of sequencing together  
```sh
cat /scratch/midway2/pflynn/15/15A_TTGCCTAG-TAAGTGGT-AHNFFJBBXX_L004_R1.fastq.gz /scratch/midway2/pflynn/15/15A_TTGCCTAG-TAAGTGGT-AHNFFJBBXX_L005_R1.fastq.gz /scratch/midway2/pflynn/15/15A_TTGCCTAG-TAAGTGGT-AHWYVLBBXX_L005_R1.fastq.gz > /scratch/midway2/pflynn/15/15A_concatenated_R1.fastq.gz
```
2. Run fastqc for manual inspection of the quality of the sequences 
```sh
mkdir fastqc_out
fastqc -t 4 raw_data/* -o fastqc_out/
```
3. trims the contactenated reads for adapters using Trimmomatic
```sh
java -jar /project2/mlcoleman/src/trimmomatic/Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 \
/scratch/midway2/pflynn/15/15A_concatenated_R1.fastq.gz /scratch/midway2/pflynn/15/15A_concatenated_R2.fastq.gz  \
/scratch/midway2/pflynn/15/15_pe1.fq /scratch/midway2/pflynn/15/15_unp1.fq /scratch/midway2/pflynn/15/15_pe2.fq /scratch/midway2/pflynn/15/15_unp2.fq \
ILLUMINACLIP:/project2/mlcoleman/src/trimmomatic/Trimmomatic-0.35/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```
 4. Builds ant species specific genome database in bowtie2
```sh
bowtie2-build /scratch/midway2/pflynn/genomes/D_quadriceps.fna /scratch/midway2/pflynn/genomes/Dquad

```
5. map ant genome to paired reads, unmapped reads are carried through to the assembly
```sh
bowtie2 -t -x /scratch/midway2/pflynn/genomes/Dquad -1 /scratch/midway2/pflynn/15/15_pe1.fq  -2 /scratch/midway2/pflynn/15/15_pe2.fq --un-conc /scratch/midway2/pflynn/15/15_conc_unmapped.fastq --al-conc /scratch/midway2/pflynn/15/15_conc_mapped.sam

```
6. spades assembly with the unmapped reads
```sh
/project2/mlcoleman/src/SPAdes-3.14.0-Linux/bin/spades.py --sc --pe1-1 /scratch/midway2/pflynn/15/15_conc_unmapped.1.fastq --pe1-2 /scratch/midway2/pflynn/15/15_conc_unmapped.2.fastq -k 21,33,55,77,99,127  -o /scratch/midway2/pflynn/15/15_spades
```

7. only contigs 300 bp or larger
```sh
perl /scratch/midway2/pflynn/removesmalls.pl 300 /scratch/midway2/pflynn/15/15_spades/scaffolds.fasta > /scratch/midway2/pflynn/Scaffolds/15_scaffolds_300.fasta
```

8. map contigs back to reads for that sample to assess read coverage
```sh
bowtie2-build /scratch/midway2/pflynn/15/15_scaffolds_300.fasta /scratch/midway2/pflynn/15/15_scaffolds_300.fasta

bowtie2 -p 12 -x /scratch/midway2/pflynn/15/15_scaffolds_300.fasta -1 /scratch/midway2/pflynn/15/15_pe1.fq  -2 /scratch/midway2/pflynn/15/15_pe2.fq -S /scratch/midway2/pflynn/15/15_reads.map.sam

samtools faidx /scratch/midway2/pflynn/15/15_scaffolds_300.fasta
```
## Scaffolds workflow 
Combine all scaffolds together 

9. CD hit for scaffolds i.e. 1_scaffolds_300.fasta
```sh
for file in /Scaffolds_CDHIT/*.fasta;
do
echo "$file";
/cd-hit-v4.8.1-2019-0228/cd-hit -i "$file" -o "${file//_merged.fasta}"  -aS 0.95 -c 0.95 -n 5 -d 0;
done
```