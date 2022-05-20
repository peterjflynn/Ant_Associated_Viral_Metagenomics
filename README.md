# Ant Associated Viral Metagenomics Workflow
## Dependencies 
Java <br>
Perl <br>
Bowtie2 <br>
samtools <br>
trimmomatic-0.35 <br>
Spades 3.14.0 <br>
CheckV <br>
CD Hit <br>

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

7. only contigs 500 bp or larger
```sh
perl /scratch/midway2/pflynn/removesmalls.pl 500 /scratch/midway2/pflynn/15/15_spades/scaffolds.fasta > /scratch/midway2/pflynn/Scaffolds/15_scaffolds_300.fasta
```

8. map contigs back to reads for that sample to assess read coverage
```sh
bowtie2-build /scratch/midway2/pflynn/15/15_scaffolds_300.fasta /scratch/midway2/pflynn/15/15_scaffolds_300.fasta

bowtie2 -p 12 -x /scratch/midway2/pflynn/15/15_scaffolds_300.fasta -1 /scratch/midway2/pflynn/15/15_pe1.fq  -2 /scratch/midway2/pflynn/15/15_pe2.fq -S /scratch/midway2/pflynn/15/15_reads.map.sam

samtools faidx /scratch/midway2/pflynn/15/15_scaffolds_300.fasta
```
## Scaffolds workflow 

9. CD hit for scaffolds i.e. 1_scaffolds_300.fasta
```sh
work_dir="./Scaffolds_CDHIT/"
read_dir="../output"

cd "${work_dir}"

for i in *.fasta
do

../cd-hit-v4.8.1-2019-0228/cd-hit -i "${i}" -o "${read_dir}/${i//}_cdhit.fasta" -aS 0.95 -c 0.95 -n 5 -d 0

done

```

10. concatenate all hits together
```sh
cat * > all_contigs_cdhit_decon.fasta
```
11. sort contigs by length
```sh
~/bbmap_2/sortbyname.sh in=all_contigs_cdhit.fasta out=all_contigs_cdhit_decon_sorted.fasta length descending
```
12. make single line
```sh
perl -pe '/^>/ ? print "\n" : chomp’  /all_contigs_cdhit_decon_sorted.fasta >  /all_contigs_cdhit_decon_sorted_single.fasta
```
# Virsorter2  
13. Use VirSorter2 version 2.1.0 to identify further viral contigs from your cross-assembled samples. I found that the CyVerse version of VirSorter2 worked better than the command line versions. Did not change any parameters Output file is: virsorter_contigs.fa

#Decontaminate
database with taxonomy for decontamination, these files are several gb, but you can download them to your server from NCBI. This tutorial is helpful: https://andreirozanski.com/2020/01/03/building-a-diamond-db-using-refseq-protein/

```sh
~/diamond/diamond makedb --in ~/23_scaffolds_300.fasta_cdhit.fasta --db ~/decontamination_db --taxonmap ~/nr/prot.accession2taxid.gz --taxonnodes ~/nr/nodes.dmp --taxonnames  ~/nr/names.dmp --threads 20 &
```
decontamination
```sh
~/diamond/diamond makedb --in ~/decontamination/23_scaffolds_300.fasta_cdhit.fasta -d ~/decontamination/decontamination_db1
makeblastdb -in ~/decontamination/23_scaffolds_300.fasta_cdhit.fasta -out ~/Decontamination/Decon -dbtype nucl -input_type fasta
```

blastn evalue 1e-5 outfmt 6 to decontaminate contigs 300 with 300
```sh
blastn -num_threads "40" -db ~/Decontamination/Decon -outfmt "6" -max_target_seqs "1" -evalue "1e-5" -max_hsps 1  -out ~/decontamination/contaminated_contigs_300.out -query ~/decontamination/all_contigs_cdhit_decon_sorted.fasta &
```

single line fasta
```sh
perl -pe '$. > 1 and /^>/ ? print "\n" : chomp' all_contigs_cdhit_decon_sorted.fasta > all_contigs_cdhit_decon_sorted_single.fasta
```
delete contaminated sequences from contig file
```sh
awk 'BEGIN{while((getline<"contam.txt")>0)l[">"$1]=1}/^>/{f=!l[$1]}f' all_contigs_cdhit_decon_sorted_single.fasta > all_contigs_300_decontam_cdhit_single.fasta
```
blastx with contig list on RefSeq database
```sh
~/diamond/diamond blastx -p "40" -d ~/refseq_protein_vir_diamond.dmnd -f "6" qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms skingdoms sphylums stitle qtitle qstrand -k "1" --evalue "1e-3" --max-hsps 1 --sensitive -o ~/300_blast/RefSeq_blastx_contigs_300.out -q ~/300_blast/all_contigs_300_decontam_cdhit_single.fasta &
```
extract viral fasta sequences from  blastx refseq results (to make txt file need to do this in text wrangler)
```sh
awk -F'>' 'NR==FNR{ids[$0]; next} NF>1{f=($2 in ids)} f' ~/300_blast/RefSeq_viral_contigs.txt ~/300_blast/all_contigs_300_decontam_cdhit_single.fasta > ~/300_blast/viral_contigs_refseq_300.fasta
```
delete anything after | in fasta header
```sh
cut -d'|' -f1 virsorter_contigs.fa > virsorter_contigs_1.fa
```
find in common headers
```sh
awk '/^>/{if (a[$1]>=1){print $1}a[$1]++}' virsorter_contigs_1.fa viral_contigs_refseq_300.fasta > common_viral.txt
```
takes away > from each line
```sh
perl -pe 's,.*>,,' common_viral.txt > common_viral2.txt
```

delete same sequences from virsorter contig file
```sh
awk 'BEGIN{while((getline<"common_viral2.txt")>0)l[">"$1]=1}/^>/{f=!l[$1]}f' virsorter_contigs_1.fa > virsorter_unique_viruses.fa
```
combine viral RefSeq contigs and Virsorter2 contigs
```sh
cat virsorter_unique_viruses.fa viral_contigs_refseq_300.fasta > Refseq_virsorter_contigs.fasta
```

blastx viral RefSeq and VirSorter contigs on nr
```sh
~/diamond/diamond blastx -p "50" -d ~/nr/nr_diamond.dmnd -f "6" qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms skingdoms sphylums stitle qtitle qstrand -k "1" --evalue "1e-3" --max-hsps 1 --sensitive -o ~/300_blast/NR_blastx_contigs_300.out -q ~/300_blast/Refseq_virsorter_contigs.fasta &
```
extract viral fasta sequences from  blastx NR results (to make txt file need to do this in text wrangler)
```sh
awk -F'>' 'NR==FNR{ids[$0]; next} NF>1{f=($2 in ids)} f' ~/300_blast/viruses_NR_300.txt ~/300_blast/Refseq_virsorter_contigs.fasta > ~/300_blast/final_NR_viral_contigs.fasta
```
CHECKV to look for proviral contamination (on home desktop)
```sh
export CHECKVDB=~/checkv-db-v1.0
checkv contamination ~/final_NR_viral_contigs.fasta  ~/checkv_output2
```

extract viral fasta sequences from  blastx NR results with proviruses and retroviruses and endogenous viruses removed (to make txt file need to do this in text wrangler) and 500 bp
```sh
awk -F'>' 'NR==FNR{ids[$0]; next} NF>1{f=($2 in ids)} f' ~/final_viruses/final_viruses.txt ~/300_blast/final_NR_viral_contigs.fasta > ~/final_viruses/final_viruses_aftertaxonomy.fasta
```


#Phylogenetics Workflow 

Geneious for manual alignment 

![Happy Christmas](atta.png)
