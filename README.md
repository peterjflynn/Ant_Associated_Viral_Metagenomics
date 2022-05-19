# Viral_Metagenomics_Ants
#Workflow
This workflow starts with raw paired-end HiSeq data in demultiplexed FASTQ formated assumed to be located within a folder called raw_seq

## Concatenates the reads from different runs  
```sh
cat /scratch/midway2/pflynn/15/15A_TTGCCTAG-TAAGTGGT-AHNFFJBBXX_L004_R1.fastq.gz /scratch/midway2/pflynn/15/15A_TTGCCTAG-TAAGTGGT-AHNFFJBBXX_L005_R1.fastq.gz /scratch/midway2/pflynn/15/15A_TTGCCTAG-TAAGTGGT-AHWYVLBBXX_L005_R1.fastq.gz > /scratch/midway2/pflynn/15/15A_concatenated_R1.fastq.gz
```
# this is next
```sh
summary(cars)
```